###############################################################################
#  RegressionExtensions.R
#  Price Elasticity of Demand for Marijuana in Mexico
#  -----------------------------------------------------------------------
#  Additional specifications that strengthen identification and credibility.
#  Meant to run AFTER RegressionAnalysis.R and IV_SeizureGravity.R.
#
#  Extensions implemented:
#    A. Municipality fixed effects (finer geography)
#    B. Hausman-type leave-one-out mean price instrument
#    C. Reduced-form estimates (regress quantity directly on instruments)
#    D. Wu-Hausman endogeneity test (OLS vs IV)
#    E. Anderson-Rubin weak-IV robust inference
#    F. Nonlinear price effects (quadratic in ln_price)
#    G. Combined IV table: OSM driving-time + seizure-gravity side by side
#
#  All standard errors clustered at state level unless noted.
###############################################################################

# ---- 0. Packages ---------------------------------------------------------
project_lib <- file.path(getwd(), ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(fixest)
library(modelsummary)
library(lubridate)
library(stringr)
library(sf)

fixest::setFixest_estimation(fixef.rm = "infinite_coef")
options(
  modelsummary_factory_latex = "kableExtra",
  modelsummary_format_numeric_latex = "plain"
)
source(here("Code", "OSMTravelTime.R"))

insert_label_after_caption <- function(path, label) {
  lines <- readLines(path, warn = FALSE)
  cap_idx <- grep("^\\\\caption\\{", lines)[1]

  if (is.na(cap_idx)) {
    stop("Could not find caption in ", path)
  }

  if (any(grepl("^\\\\label\\{", lines))) {
    return(invisible(NULL))
  }

  lines <- append(lines, paste0("\\label{", label, "}"), after = cap_idx)
  writeLines(lines, path)
}

force_table_here <- function(path) {
  lines <- readLines(path, warn = FALSE)
  table_idx <- grep("^\\\\begin\\{table\\}", lines)[1]

  if (is.na(table_idx)) {
    stop("Could not find table environment in ", path)
  }

  lines[table_idx] <- "\\begin{table}[H]"
  writeLines(lines, path)
}

wrap_tabular_in_resizebox <- function(path) {
  lines <- readLines(path, warn = FALSE)
  begin_idx <- grep("^\\\\begin\\{tabular", lines)[1]
  end_idx <- grep("^\\\\end\\{tabular\\}", lines)[1]

  if (is.na(begin_idx) || is.na(end_idx)) {
    stop("Could not find tabular environment in ", path)
  }

  lines <- append(lines, "\\resizebox{\\linewidth}{!}{%", after = begin_idx - 1)
  end_idx <- end_idx + 1
  lines <- append(lines, "}", after = end_idx)
  writeLines(lines, path)
}

educ_coef_map <- c(
  "educ_fMiddle school" = "Education: Middle school",
  "educ_fHigh school" = "Education: High school",
  "educ_fBachelor" = "Education: Bachelor",
  "educ_fGraduate" = "Education: Graduate"
)

gender_coef_map <- c(
  "gender_fOther" = "Gender: Other",
  "gender_fWomen" = "Gender: Women"
)

# ---- 1. Load and prepare data (mirrors RegressionAnalysis.R) -------------
df <- read_parquet(here("Data", "Marijuana_Prices_in_Mexico_clean.parquet"))

df <- df |>
  mutate(
    ln_quantity  = log(quantity_g),
    ln_price     = log(price_mxn_g),
    quality_good = as.integer(quality == "good"),
    cvegeo       = str_pad(as.character(muni_code), width = 5, pad = "0"),
    state_code   = str_sub(cvegeo, 1, 2),
    purchase_date = as.Date(purchase_date),
    month         = floor_date(purchase_date, "month"),
    month_str     = format(month, "%Y-%m"),
    month_num     = month(purchase_date),
    month_id      = year(purchase_date) * 12L + month(purchase_date),
    state_month   = paste0(state_code, "_", month_str),
    educ_f        = factor(education,
                           levels = c("Elementary school", "Middle school",
                                      "High school", "Bachelor", "Graduate")),
    gender_f      = factor(gender)
  )

# Trimming flags
p01 <- quantile(df$price_mxn_g, 0.01)
p99 <- quantile(df$price_mxn_g, 0.99)
q01 <- quantile(df$quantity_g, 0.01)
q99 <- quantile(df$quantity_g, 0.99)

df <- df |>
  mutate(
    in_trim = price_mxn_g >= p01 & price_mxn_g <= p99 &
              quantity_g  >= q01 & quantity_g  <= q99
  )

cat("Analytic sample:", nrow(df), "observations\n")
cat("States:", n_distinct(df$state_code),
    "| Municipalities:", n_distinct(df$cvegeo), "\n")

# ---- 2. OSM driving-time IV (from RegressionAnalysis.R logic) -------------
mun_centroids <- municipality_centroids_osm(
  here("Data/Shapefiles/conjunto_de_datos", "00mun.shp")
)

hub_state_codes <- c("25", "08", "12", "10", "16")

hub_points <- hub_points_from_states(
  here("Data/Shapefiles/conjunto_de_datos", "00ent.shp"),
  hub_state_codes
)

df_cvegeo <- df |> distinct(cvegeo)

df_hub_time <- tryCatch(
  df_cvegeo |>
    left_join(mun_centroids, by = "cvegeo") |>
    nearest_hub_travel_time(
      hub_points = hub_points,
      cache_path = here("Data", "Derived", "osrm_statehub_travel_time_min.rds")
    ),
  error = function(e) {
    warning(
      "OSM production-hub driving-time IV unavailable. ",
      conditionMessage(e),
      call. = FALSE
    )
    tibble(
      cvegeo = df_cvegeo$cvegeo,
      travel_time_nearest_hub_min = NA_real_,
      nearest_hub_state = NA_character_,
      ln_travel_time_hub = NA_real_
    )
  }
)

df <- df |>
  left_join(df_hub_time, by = "cvegeo")

has_osm_iv <- any(is.finite(df$ln_travel_time_hub))

# ---- 3. Seizure-gravity IV (from IV_SeizureGravity.R logic) -------------
#  Load pre-built gravity IV if available; otherwise skip seizure IV models.
mucd_dir <- here("Data", "MUCD")
has_mucd <- dir.exists(mucd_dir) && length(list.files(mucd_dir, "\\.csv$")) > 0

if (has_mucd) {
  cat("Loading MUCD seizure data for gravity IV...\n")
  source_env <- new.env()

  # --- inline gravity construction (self-contained) ---
  mucd_files <- list.files(mucd_dir, pattern = "\\.csv$", full.names = TRUE)

  read_mucd <- function(fp) {
    readr::read_csv(
      fp,
      skip = 3,
      col_names = c(
        "CodEntidad",
        "nombreEntidad",
        "nombreMunicipio",
        "idMunicipio",
        "Anio",
        "Mes",
        "Destrucciondecultivos_marihuana",
        "Aseguramientode_marihuana",
        "Aseguramientode_semillademarihuana",
        "trailing_empty"
      ),
      col_types = readr::cols(.default = readr::col_character()),
      locale = readr::locale(encoding = "latin1"),
      show_col_types = FALSE
    )
  }
  seizure_raw <- bind_rows(lapply(mucd_files, read_mucd))

  seizures <- seizure_raw |>
    transmute(
      cvegeo_seizure = str_pad(as.character(idMunicipio), 5, pad = "0"),
      anio = as.integer(Anio),
      mes = as.integer(Mes),
      month_id = as.integer(Anio) * 12L + as.integer(Mes),
      destruccion_kg = readr::parse_number(
        as.character(Destrucciondecultivos_marihuana),
        na = c("", "NA")
      ),
      seizure_kg = readr::parse_number(
        as.character(Aseguramientode_marihuana),
        na = c("", "NA")
      ),
      semilla_kg = readr::parse_number(
        as.character(Aseguramientode_semillademarihuana),
        na = c("", "NA")
      )
    ) |>
    mutate(
      across(c(destruccion_kg, seizure_kg, semilla_kg), ~ replace_na(.x, 0)),
      total_mj_kg = destruccion_kg + seizure_kg + semilla_kg
    )

  seizure_mm <- seizures |>
    filter(!is.na(cvegeo_seizure), nchar(cvegeo_seizure) == 5) |>
    group_by(cvegeo_seizure, month_id) |>
    summarise(
      seizure_kg = sum(seizure_kg, na.rm = TRUE),
      total_mj_kg = sum(total_mj_kg, na.rm = TRUE),
      .groups = "drop"
    )

  # Centroids and OSRM driving-time matrix for gravity
  purchase_munis <- unique(df$cvegeo)
  seizure_cvegeos <- seizure_mm |>
    filter(seizure_kg > 0 | total_mj_kg > 0) |>
    distinct(cvegeo_seizure) |>
    pull(cvegeo_seizure)
  cent_p <- mun_centroids |> filter(cvegeo %in% purchase_munis)
  cent_s <- mun_centroids |> filter(cvegeo %in% seizure_cvegeos)

  time_sp <- osrm_duration_matrix(
    src = cent_s,
    dst = cent_p,
    cache_path = here("Data", "Derived", "osrm_seizure_purchase_time_min.rds"),
    src_id = "cvegeo",
    dst_id = "cvegeo"
  )
  time_ps <- t(time_sp)

  purchase_months <- sort(unique(df$month_id))
  gravity_iv <- gravity_from_travel_time(
    time_mat = time_ps,
    seizure_muni_month = seizure_mm,
    purchase_months = purchase_months,
    kg_col = "seizure_kg",
    out_col = "Z_gravity_lag1",
    lag_months = 1,
    missing_lag = "na",
    month_col = "month_id"
  )

  df <- df |>
    left_join(gravity_iv, by = c("cvegeo", "month_id")) |>
    mutate(
      ln_Z_gravity_lag1 = log(1 + Z_gravity_lag1)
    )
  cat("Lagged gravity IV merged. Non-missing:",
      sum(!is.na(df$Z_gravity_lag1)), "\n")
} else {
  cat("MUCD data not found; seizure-gravity IV models will be skipped.\n")
  df$Z_gravity_lag1 <- NA_real_
  df$ln_Z_gravity_lag1 <- NA_real_
}

# ============================================================================
# A. MUNICIPALITY FIXED EFFECTS
# ============================================================================
cat("\n\n========= A. MUNICIPALITY FIXED EFFECTS =========\n")

# Drop municipalities with fewer than 2 observations (singleton FE)
muni_counts <- df |> count(cvegeo)
df <- df |> left_join(muni_counts, by = "cvegeo", suffix = c("", "_muni"))

# State-FE comparison
mA_state <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
                  | state_code,
                  data = df, cluster = ~state_code)

# State-by-month FE benchmark
mA_state_month <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
                        | state_month,
                        data = df, cluster = ~state_code)

# Municipality FE added to the state-by-month benchmark.
mA_muni_state_month <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
                             | state_month + cvegeo,
                             data = filter(df, n >= 2),
                             cluster = ~state_code)

cat("\n  State FE:        beta =", round(coef(mA_state)["ln_price"], 4),
    " SE =", round(se(mA_state)["ln_price"], 4), "\n")
cat("  State x Month FE: beta =", round(coef(mA_state_month)["ln_price"], 4),
    " SE =", round(se(mA_state_month)["ln_price"], 4), "\n")
cat("  State x Month + Municipality FE: beta =",
    round(coef(mA_muni_state_month)["ln_price"], 4),
    " SE =", round(se(mA_muni_state_month)["ln_price"], 4), "\n")

# ============================================================================
# B. HAUSMAN-TYPE LEAVE-ONE-OUT MEAN PRICE INSTRUMENT
# ============================================================================
cat("\n\n========= B. HAUSMAN-TYPE LEAVE-ONE-OUT PRICE IV =========\n")
#
# Instrument: for each observation i in state s, quality q, month t:
#   Z_i = mean(ln_price) among all OTHER observations in the same
#         (state, quality, month) cell, excluding observation i.
#
# Identifying assumption (à la Hausman 1996, Nevo 2001):
#   RELEVANCE: common supply-cost shocks within a state-quality-month
#     cell drive shared price variation → strong first stage.
#   EXCLUSION: after controlling for individual covariates and FE,
#     the cell-level mean price (excluding own) reflects supply-side
#     factors (shared costs, local market conditions) rather than
#     individual demand shocks.  Valid if individual demand shocks are
#     independent within cells.

df <- df |>
  group_by(state_code, quality_good, month_str) |>
  mutate(
    cell_n       = n(),
    cell_sum_lnp = sum(ln_price),
    # Leave-one-out mean
    loo_mean_lnp = (cell_sum_lnp - ln_price) / (cell_n - 1)
  ) |>
  ungroup()

# Drop cells with only 1 observation (no leave-one-out possible)
df_hausman <- df |> filter(cell_n >= 2)

cat("Hausman IV sample:", nrow(df_hausman),
    "(dropped", nrow(df) - nrow(df_hausman), "singleton cells)\n")

# IV regression: Hausman instrument (state-by-month FE for comparability with m6)
mB_iv <- feols(ln_quantity ~ quality_good + age + gender_f + educ_f
               | state_month
               | ln_price ~ loo_mean_lnp,
               data = df_hausman, cluster = ~state_code)

# First stage
mB_fs <- feols(ln_price ~ loo_mean_lnp + quality_good + age + gender_f + educ_f
               | state_month,
               data = df_hausman, cluster = ~state_code)

cat("  Hausman IV:  beta =", round(coef(mB_iv)["fit_ln_price"], 4),
    " SE =", round(se(mB_iv)["fit_ln_price"], 4), "\n")
cat("  First stage: coef =", round(coef(mB_fs)["loo_mean_lnp"], 4),
    " t =", round(coef(mB_fs)["loo_mean_lnp"] / se(mB_fs)["loo_mean_lnp"], 2), "\n")
cat("  First-stage F:\n  ")
print(fitstat(mB_iv, "ivf"))

# ============================================================================
# C. REDUCED-FORM ESTIMATES
# ============================================================================
cat("\n\n========= C. REDUCED-FORM ESTIMATES =========\n")
#
# Regress ln_quantity directly on each instrument.
# Under the IV assumptions, the reduced-form coefficient = beta * pi,
# where beta is the elasticity and pi is the first-stage coefficient.
# A significant reduced form is necessary for the IV to be informative.

# C1: OSM driving-time instrument
if (has_osm_iv) {
  mC_rf_dist <- feols(ln_quantity ~ ln_travel_time_hub + quality_good + age + gender_f + educ_f
                      | state_month,
                      data = df, cluster = ~state_code)

  cat("  Reduced form (OSM driving time):",
      "coef =", round(coef(mC_rf_dist)["ln_travel_time_hub"], 5),
      "SE =", round(se(mC_rf_dist)["ln_travel_time_hub"], 5),
      "t =", round(coef(mC_rf_dist)["ln_travel_time_hub"] /
                    se(mC_rf_dist)["ln_travel_time_hub"], 2), "\n")
} else {
  mC_rf_dist <- NULL
  cat("  Reduced form (OSM driving time): skipped; no cached state-hub OSRM matrix.\n")
}

# C2: seizure-gravity instrument
if (has_mucd) {
  mC_rf_grav <- feols(ln_quantity ~ ln_Z_gravity_lag1 + quality_good + age + gender_f + educ_f
                      | state_month,
                      data = df, cluster = ~state_code)
  cat("  Reduced form (gravity):",
      "coef =", round(coef(mC_rf_grav)["ln_Z_gravity_lag1"], 5),
      "SE =", round(se(mC_rf_grav)["ln_Z_gravity_lag1"], 5),
      "t =", round(coef(mC_rf_grav)["ln_Z_gravity_lag1"] /
                    se(mC_rf_grav)["ln_Z_gravity_lag1"], 2), "\n")
}

# C3: Hausman instrument
mC_rf_haus <- feols(ln_quantity ~ loo_mean_lnp + quality_good + age + gender_f + educ_f
                    | state_month,
                    data = df_hausman, cluster = ~state_code)
cat("  Reduced form (Hausman):",
    "coef =", round(coef(mC_rf_haus)["loo_mean_lnp"], 5),
    "SE =", round(se(mC_rf_haus)["loo_mean_lnp"], 5),
    "t =", round(coef(mC_rf_haus)["loo_mean_lnp"] /
                  se(mC_rf_haus)["loo_mean_lnp"], 2), "\n")

# ============================================================================
# D. WU-HAUSMAN ENDOGENEITY TEST
# ============================================================================
cat("\n\n========= D. WU-HAUSMAN ENDOGENEITY TEST =========\n")
#
# Test H0: OLS is consistent (price is exogenous)
# If rejected, the IV estimates are preferred.
#
# Implementation: add the first-stage residuals to the OLS equation.
# If the coefficient on the residuals is significant → endogeneity.

# Using OSM driving-time IV
if (has_osm_iv) {
  fs_dist <- feols(ln_price ~ ln_travel_time_hub + quality_good + age + gender_f + educ_f
                   | state_month,
                   data = df, cluster = ~state_code)

  df$resid_fs_dist <- residuals(fs_dist)

  mD_wu <- feols(ln_quantity ~ ln_price + resid_fs_dist + quality_good + age +
                   gender_f + educ_f
                 | state_month,
                 data = df, cluster = ~state_code)

  cat("  Wu-Hausman (OSM driving-time IV):\n")
  cat("    Coef on FS residual:", round(coef(mD_wu)["resid_fs_dist"], 4),
      " SE:", round(se(mD_wu)["resid_fs_dist"], 4),
      " t:", round(coef(mD_wu)["resid_fs_dist"] /
                    se(mD_wu)["resid_fs_dist"], 2), "\n")
  cat("    p-value:", round(2 * pt(-abs(coef(mD_wu)["resid_fs_dist"] /
                                         se(mD_wu)["resid_fs_dist"]),
                                   df = mD_wu$nobs - length(coef(mD_wu))), 4), "\n")
} else {
  mD_wu <- NULL
  cat("  Wu-Hausman (OSM driving-time IV): skipped; no cached state-hub OSRM matrix.\n")
}

# Using Hausman instrument
fs_haus <- feols(ln_price ~ loo_mean_lnp + quality_good + age + gender_f + educ_f
                 | state_month,
                 data = df_hausman, cluster = ~state_code)

df_hausman$resid_fs_haus <- residuals(fs_haus)

mD_wu_haus <- feols(ln_quantity ~ ln_price + resid_fs_haus + quality_good + age +
                      gender_f + educ_f
                    | state_month,
                    data = df_hausman, cluster = ~state_code)

cat("  Wu-Hausman (Hausman IV):\n")
cat("    Coef on FS residual:", round(coef(mD_wu_haus)["resid_fs_haus"], 4),
    " SE:", round(se(mD_wu_haus)["resid_fs_haus"], 4),
    " t:", round(coef(mD_wu_haus)["resid_fs_haus"] /
                  se(mD_wu_haus)["resid_fs_haus"], 2), "\n")
wu_p_haus <- 2 * pt(-abs(coef(mD_wu_haus)["resid_fs_haus"] /
                           se(mD_wu_haus)["resid_fs_haus"]),
                    df = mD_wu_haus$nobs - length(coef(mD_wu_haus)))
cat("    p-value:", round(wu_p_haus, 4), "\n")

if (has_mucd) {
  df_grav_diag <- df |>
    filter(!is.na(ln_Z_gravity_lag1))

  fs_grav_diag <- feols(ln_price ~ ln_Z_gravity_lag1 + quality_good + age + gender_f + educ_f
                        | state_month,
                        data = df_grav_diag, cluster = ~state_code)
  df_grav_diag$resid_fs_grav <- residuals(fs_grav_diag)

  mD_wu_grav <- feols(ln_quantity ~ ln_price + resid_fs_grav + quality_good + age +
                        gender_f + educ_f
                      | state_month,
                      data = df_grav_diag, cluster = ~state_code)

  cat("  Wu-Hausman (gravity IV):\n")
  cat("    Coef on FS residual:", round(coef(mD_wu_grav)["resid_fs_grav"], 4),
      " SE:", round(se(mD_wu_grav)["resid_fs_grav"], 4),
      " t:", round(coef(mD_wu_grav)["resid_fs_grav"] /
                    se(mD_wu_grav)["resid_fs_grav"], 2), "\n")
  wu_p_grav <- 2 * pt(-abs(coef(mD_wu_grav)["resid_fs_grav"] /
                             se(mD_wu_grav)["resid_fs_grav"]),
                      df = mD_wu_grav$nobs - length(coef(mD_wu_grav)))
  cat("    p-value:", round(wu_p_grav, 4), "\n")
}

# ============================================================================
# E. ANDERSON-RUBIN WEAK-IV ROBUST INFERENCE
# ============================================================================
cat("\n\n========= E. ANDERSON-RUBIN CONFIDENCE SETS =========\n")
#
# The AR test inverts the hypothesis test for beta_0:
#   H0: beta = beta_0
# by regressing (Y - beta_0 * X) on Z and testing joint significance.
# The 95% confidence set is all beta_0 values not rejected at 5%.
#
# This is robust to weak instruments.

ar_grid <- seq(-2.5, 0.5, by = 0.01)

# AR test using OSM driving-time IV
if (has_osm_iv) {
  ar_pvals_dist <- sapply(ar_grid, function(b0) {
    df$y_tilde <- df$ln_quantity - b0 * df$ln_price
    m_ar <- feols(y_tilde ~ ln_travel_time_hub + quality_good + age + gender_f + educ_f
                  | state_month,
                  data = df, cluster = ~state_code)
    # Wald test on ln_travel_time_hub coefficient
    t_stat <- coef(m_ar)["ln_travel_time_hub"] / se(m_ar)["ln_travel_time_hub"]
    2 * pt(-abs(t_stat), df = m_ar$nobs - length(coef(m_ar)))
  })

  ar_ci_dist <- range(ar_grid[ar_pvals_dist >= 0.05])
  cat("  AR 95% CI (OSM driving-time IV): [", round(ar_ci_dist[1], 3), ",",
      round(ar_ci_dist[2], 3), "]\n")
} else {
  ar_ci_dist <- c(NA_real_, NA_real_)
  cat("  AR 95% CI (OSM driving-time IV): skipped; no cached state-hub OSRM matrix.\n")
}

# AR test using Hausman IV
ar_pvals_haus <- sapply(ar_grid, function(b0) {
  df_hausman$y_tilde <- df_hausman$ln_quantity - b0 * df_hausman$ln_price
  m_ar <- feols(y_tilde ~ loo_mean_lnp + quality_good + age + gender_f + educ_f
                | state_month,
                data = df_hausman, cluster = ~state_code)
  t_stat <- coef(m_ar)["loo_mean_lnp"] / se(m_ar)["loo_mean_lnp"]
  2 * pt(-abs(t_stat), df = m_ar$nobs - length(coef(m_ar)))
})

ar_ci_haus <- range(ar_grid[ar_pvals_haus >= 0.05])
cat("  AR 95% CI (Hausman IV):  [", round(ar_ci_haus[1], 3), ",",
    round(ar_ci_haus[2], 3), "]\n")

# AR test using seizure-gravity IV
if (has_mucd) {
  df_grav_ar <- df |>
    filter(!is.na(ln_Z_gravity_lag1))

  ar_pvals_grav <- sapply(ar_grid, function(b0) {
    df_grav_ar$y_tilde <- df_grav_ar$ln_quantity - b0 * df_grav_ar$ln_price
    m_ar <- feols(y_tilde ~ ln_Z_gravity_lag1 + quality_good + age + gender_f + educ_f
                  | state_month,
                  data = df_grav_ar, cluster = ~state_code)
    t_stat <- coef(m_ar)["ln_Z_gravity_lag1"] / se(m_ar)["ln_Z_gravity_lag1"]
    2 * pt(-abs(t_stat), df = m_ar$nobs - length(coef(m_ar)))
  })
  ar_ci_grav <- range(ar_grid[ar_pvals_grav >= 0.05])
  cat("  AR 95% CI (gravity IV):  [", round(ar_ci_grav[1], 3), ",",
      round(ar_ci_grav[2], 3), "]\n")
}

# ============================================================================
# F. NONLINEAR PRICE EFFECTS
# ============================================================================
cat("\n\n========= F. NONLINEAR PRICE EFFECTS =========\n")
#
# Add (ln_price)^2 to test for curvature in the demand relationship.
# If coefficient on the squared term is significant, the constant-elasticity
# specification may be misspecified.

df <- df |> mutate(ln_price_sq = ln_price^2)

mG_quad <- feols(ln_quantity ~ ln_price + ln_price_sq + quality_good + age +
                   gender_f + educ_f
                 | state_month,
                 data = df, cluster = ~state_code)

cat("  Quadratic specification:\n")
cat("    ln_price:    coef =", round(coef(mG_quad)["ln_price"], 4),
    " SE =", round(se(mG_quad)["ln_price"], 4), "\n")
cat("    ln_price^2:  coef =", round(coef(mG_quad)["ln_price_sq"], 5),
    " SE =", round(se(mG_quad)["ln_price_sq"], 5), "\n")
cat("    Implied elasticity at median ln_price (",
    round(median(df$ln_price), 2), "):",
    round(coef(mG_quad)["ln_price"] + 2 * coef(mG_quad)["ln_price_sq"] *
            median(df$ln_price), 4), "\n")

# ============================================================================
# G. COMBINED IV SUMMARY TABLE
# ============================================================================
cat("\n\n========= G. COMBINED IV SUMMARY TABLE =========\n")

# OLS baseline
mH_ols <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
                | state_month,
                data = df, cluster = ~state_code)

# Hausman IV
mH_iv_haus <- feols(ln_quantity ~ quality_good + age + gender_f + educ_f
                    | state_month
                    | ln_price ~ loo_mean_lnp,
                    data = df_hausman, cluster = ~state_code)

tab_combined <- list(
  "(1) OLS"         = mH_ols
)

if (has_osm_iv) {
  mH_iv_dist <- feols(ln_quantity ~ quality_good + age + gender_f + educ_f
                      | state_month
                      | ln_price ~ ln_travel_time_hub,
                      data = df, cluster = ~state_code)
  tab_combined[["(2) IV: OSM driving time"]] <- mH_iv_dist
} else {
  mH_iv_dist <- NULL
}

tab_combined[[paste0("(", length(tab_combined) + 1L, ") IV: Hausman")]] <- mH_iv_haus

# Add seizure-gravity if available
if (has_mucd) {
  mH_iv_grav <- feols(ln_quantity ~ quality_good + age + gender_f + educ_f
                      | state_month
                      | ln_price ~ ln_Z_gravity_lag1,
                      data = df, cluster = ~state_code)
  tab_combined[[paste0("(", length(tab_combined) + 1L, ") IV: gravity")]] <- mH_iv_grav
}

msummary(tab_combined, output = "markdown",
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("ln_price" = "ln(Price)",
                      "fit_ln_price" = "ln(Price) [IV]",
                      "quality_good" = "Good quality",
                      "age" = "Age",
                      gender_coef_map,
                      educ_coef_map),
         gof_map = c("nobs", "r.squared", "adj.r.squared"))

# ---- Export combined table to LaTeX ----
output_dir <- here("LaTeX", "tables")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

combined_notes <- c(
  "Dependent variable: ln(quantity in grams).",
  "All columns include state-by-month FE and demographic controls.",
  if (has_osm_iv) {
    "OSM driving-time IV: instrument = ln(OSM driving time from nearest constructed production-hub point)."
  } else {
    "OSM hub driving-time IV omitted; regenerated state-centroid OSRM cache required."
  },
  "Hausman IV: instrument = leave-one-out mean ln(price) within state-quality-month cell.",
  "Gravity IV: instrument = ln(1 + one-month-lag seizure-gravity index).",
  "For weak geography-based IVs, inference should rely on Anderson-Rubin confidence sets rather than conventional significance stars.",
  "SEs clustered at state level.",
  "* p $<$ 0.1, ** p $<$ 0.05, *** p $<$ 0.01."
)

msummary(tab_combined,
         output = file.path(output_dir, "table_combined_iv.tex"),
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("ln_price" = "ln(Price)",
                      "fit_ln_price" = "ln(Price) [IV]",
                      "quality_good" = "Good quality",
                      "age" = "Age",
                      gender_coef_map,
                      educ_coef_map),
         gof_map = c("nobs", "r.squared", "adj.r.squared"),
         title = "OLS and IV estimates: alternative identification strategies",
         notes = combined_notes)
insert_label_after_caption(
  file.path(output_dir, "table_combined_iv.tex"),
  "tab:iv-combined"
)
force_table_here(file.path(output_dir, "table_combined_iv.tex"))
wrap_tabular_in_resizebox(file.path(output_dir, "table_combined_iv.tex"))

# ---- Export IV diagnostic table ----
extract_ivf <- function(model) {
  txt <- paste(capture.output(fitstat(model, "ivf")), collapse = " ")
  stat <- sub(".*stat = ([^,]+),.*", "\\1", txt)
  as.numeric(stat)
}

fmt_num <- function(x, digits = 3) {
  if (is.na(x) || !is.finite(x)) return("--")
  formatC(x, format = "f", digits = digits)
}

fmt_p <- function(x) {
  if (is.na(x) || !is.finite(x)) return("--")
  formatC(x, format = "f", digits = 3)
}

fmt_coef_se <- function(coef_value, se_value) {
  paste0(fmt_num(coef_value), " (", fmt_num(se_value), ")")
}

fmt_iv_coef_se <- function(model) {
  nm <- intersect(c("fit_ln_price", "ln_price"), names(coef(model)))[1]
  if (is.na(nm)) return("--")
  fmt_coef_se(coef(model)[nm], se(model)[nm])
}

fmt_ci <- function(ci) {
  if (length(ci) != 2 || any(!is.finite(ci))) return("--")
  paste0("[", fmt_num(ci[1], 2), ", ", fmt_num(ci[2], 2), "]")
}

wu_p_dist <- if (has_osm_iv) {
  2 * pt(-abs(coef(mD_wu)["resid_fs_dist"] /
                se(mD_wu)["resid_fs_dist"]),
         df = mD_wu$nobs - length(coef(mD_wu)))
} else {
  NA_real_
}

iv_diag_path <- file.path(output_dir, "table_iv_diagnostics.tex")
writeLines(
  c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Instrumental-variables estimates and diagnostics}",
    "\\label{tab:iv-diagnostics}",
    "\\centering",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}[t]{lccc}",
    "\\toprule",
    "Diagnostic & OSM driving-time IV & Hausman IV & Seizure-gravity IV\\\\",
    "\\midrule",
    paste0("First-stage F & ", if (has_osm_iv) fmt_num(extract_ivf(mH_iv_dist), 2) else "--", " & ",
           fmt_num(extract_ivf(mH_iv_haus), 2), " & ",
           if (has_mucd) fmt_num(extract_ivf(mH_iv_grav), 2) else "--", "\\\\"),
    paste0("IV estimate & ",
           if (has_osm_iv) fmt_iv_coef_se(mH_iv_dist) else "--", " & ",
           fmt_iv_coef_se(mH_iv_haus), " & ",
           if (has_mucd) fmt_iv_coef_se(mH_iv_grav) else "--", "\\\\"),
    paste0("Wu--Hausman p-value & ", fmt_p(wu_p_dist), " & ",
           fmt_p(wu_p_haus), " & ",
           if (has_mucd) fmt_p(wu_p_grav) else "--", "\\\\"),
    paste0("Anderson--Rubin 95\\% CI & ", fmt_ci(ar_ci_dist), " & ",
           fmt_ci(ar_ci_haus), " & ",
           if (has_mucd) fmt_ci(ar_ci_grav) else "--", "\\\\"),
    "\\bottomrule",
    "\\multicolumn{4}{l}{\\rule{0pt}{1em}All diagnostics use the state-by-month FE specifications in Table~\\ref{tab:iv-combined}.}\\\\",
    "\\multicolumn{4}{l}{\\rule{0pt}{1em}IV estimate entries reproduce the coefficient on ln(Price) [IV] and clustered SE from Table~\\ref{tab:iv-combined}.}\\\\",
    "\\multicolumn{4}{l}{\\rule{0pt}{1em}Seizure-gravity diagnostics use one-month-lagged MUCD marijuana seizures from December 2018--August 2019.}\\\\",
    "\\multicolumn{4}{l}{\\rule{0pt}{1em}For weak geography-based IVs, inference should rely on Anderson-Rubin confidence sets rather than conventional significance stars.}\\\\",
    "\\multicolumn{4}{l}{\\rule{0pt}{1em}OSM and seizure-gravity have first-stage F-statistics below 10; Hausman has a strong first stage but relies on a leave-one-out exclusion restriction.}\\\\",
    "\\end{tabular}",
    "}",
    "\\end{table}"
  ),
  iv_diag_path
)

# ---- Export municipality FE table ----
tab_muni_fe <- list(
  "(1) St×Mo FE"        = mA_state_month,
  "(2) St×Mo + Muni FE" = mA_muni_state_month
)

msummary(tab_muni_fe,
         output = file.path(output_dir, "table_muni_fe.tex"),
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("ln_price" = "ln(Price)",
                      "quality_good" = "Good quality",
                      "age" = "Age",
                      gender_coef_map,
                      educ_coef_map),
         gof_map = c("nobs", "r.squared", "adj.r.squared"),
         title = "Demand estimates with alternative fixed-effects structures",
         notes = c(
           "Dependent variable: ln(quantity in grams).",
           "All columns include state-by-month FE and demographic controls.",
           "Col 2 additionally includes municipality FE and restricts to municipalities with at least 2 observations.",
           "SEs clustered at state level.",
           "* p $<$ 0.1, ** p $<$ 0.05, *** p $<$ 0.01."
         ))
insert_label_after_caption(
  file.path(output_dir, "table_muni_fe.tex"),
  "tab:muni-fe"
)
force_table_here(file.path(output_dir, "table_muni_fe.tex"))
wrap_tabular_in_resizebox(file.path(output_dir, "table_muni_fe.tex"))

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n\n==================== SUMMARY ====================\n")
cat("A. Municipality FE: elasticity",
    round(coef(mA_muni_state_month)["ln_price"], 3), "(vs St x Mo benchmark:",
    round(coef(mA_state_month)["ln_price"], 3), ")\n")
cat("B. Hausman IV: elasticity",
    round(coef(mB_iv)["fit_ln_price"], 3), "\n")
cat("C. Reduced forms reported above\n")
cat("D. Wu-Hausman tests reported above\n")
cat("E. Anderson-Rubin CIs reported above\n")
cat("F. Quadratic: ln_price^2 coef =",
    round(coef(mG_quad)["ln_price_sq"], 5), "\n")
cat("G. Combined IV table exported\n")
cat("\nAll extension tables saved to:", output_dir, "\n")
cat("Done.\n")
