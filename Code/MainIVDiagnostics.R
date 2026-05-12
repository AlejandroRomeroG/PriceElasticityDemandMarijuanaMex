###############################################################################
#  MainIVDiagnostics.R
#  Main-paper IV specifications and diagnostics
#  -----------------------------------------------------------------------
#  Generates:
#    - LaTeX/tables/table_combined_iv.tex
#    - LaTeX/tables/table_iv_diagnostics.tex
#    - Data/Derived/main_iv_diagnostics.csv
###############################################################################

project_lib <- file.path(getwd(), ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(arrow)
library(dplyr)
library(fixest)
library(here)
library(lubridate)
library(modelsummary)
library(readr)
library(sf)
library(stringr)
library(tidyr)

fixest::setFixest_estimation(fixef.rm = "infinite_coef")
options(
  modelsummary_factory_latex = "kableExtra",
  modelsummary_format_numeric_latex = "plain"
)

source(here("Code", "OSMTravelTime.R"))

output_dir <- here("LaTeX", "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(here("Data", "Derived"), recursive = TRUE, showWarnings = FALSE)

insert_label_after_caption <- function(path, label) {
  lines <- readLines(path, warn = FALSE)
  cap_idx <- grep("^\\\\caption\\{", lines)[1]
  if (is.na(cap_idx)) stop("Could not find caption in ", path)
  if (any(grepl("^\\\\label\\{", lines))) return(invisible(NULL))
  lines <- append(lines, paste0("\\label{", label, "}"), after = cap_idx)
  writeLines(lines, path)
}

force_table_here <- function(path) {
  lines <- readLines(path, warn = FALSE)
  table_idx <- grep("^\\\\begin\\{table\\}", lines)[1]
  if (is.na(table_idx)) stop("Could not find table environment in ", path)
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

fmt_num <- function(x, digits = 3) {
  if (is.na(x) || !is.finite(x)) return("--")
  formatC(x, format = "f", digits = digits)
}

fmt_p <- function(x) {
  if (is.na(x) || !is.finite(x)) return("--")
  formatC(x, format = "f", digits = 3)
}

fmt_ci <- function(ci) {
  if (length(ci) != 2 || any(!is.finite(ci))) return("--")
  paste0("[", fmt_num(ci[1], 2), ", ", fmt_num(ci[2], 2), "]")
}

fmt_coef_se <- function(coef_value, se_value) {
  paste0(fmt_num(coef_value), " (", fmt_num(se_value), ")")
}

fmt_num_vec <- function(x, digits = 3) {
  vapply(x, fmt_num, character(1), digits = digits)
}

fmt_p_vec <- function(x) {
  vapply(x, fmt_p, character(1))
}

safe_coef <- function(model, term) {
  out <- coef(model)
  if (!term %in% names(out)) return(NA_real_)
  unname(out[[term]])
}

safe_se <- function(model, term) {
  out <- se(model)
  if (!term %in% names(out)) return(NA_real_)
  unname(out[[term]])
}

extract_ivf <- function(model) {
  out <- tryCatch(fitstat(model, "ivf"), error = function(e) NULL)
  if (is.null(out) || length(out) == 0) return(NA_real_)
  as.numeric(unname(out[[1]]$stat))
}

iv_coef_name <- function(model) {
  intersect(c("fit_ln_price", "ln_price"), names(coef(model)))[1]
}

fmt_iv_coef_se <- function(model) {
  nm <- iv_coef_name(model)
  if (is.na(nm)) return("--")
  fmt_coef_se(coef(model)[nm], se(model)[nm])
}

make_region <- function(state_code) {
  case_when(
    state_code %in% c("02", "03", "08", "10", "25", "26") ~ "Northwest",
    state_code %in% c("05", "19", "28") ~ "Northeast",
    state_code %in% c("01", "06", "11", "14", "16", "18", "22", "24", "32") ~ "WestBajio",
    state_code %in% c("09", "13", "15", "17", "21", "29") ~ "Central",
    state_code %in% c("07", "12", "20", "30") ~ "South",
    state_code %in% c("04", "23", "27", "31") ~ "Southeast",
    TRUE ~ "Other"
  )
}

read_mucd <- function(filepath) {
  read_csv(
    filepath,
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
    col_types = cols(.default = col_character()),
    locale = locale(encoding = "latin1"),
    show_col_types = FALSE
  )
}

# ---- Analytic sample ---------------------------------------------------------
df <- read_parquet(here("Data", "Marijuana_Prices_in_Mexico_clean.parquet")) |>
  mutate(
    ln_quantity = log(quantity_g),
    ln_price = log(price_mxn_g),
    quality_good = as.integer(quality == "good"),
    cvegeo = str_pad(as.character(muni_code), width = 5, pad = "0"),
    state_code = str_sub(cvegeo, 1, 2),
    purchase_date = as.Date(purchase_date),
    month = floor_date(purchase_date, "month"),
    month_str = format(month, "%Y-%m"),
    month_id = year(purchase_date) * 12L + month(purchase_date),
    state_month = paste0(state_code, "_", month_str),
    quality_month = paste0(quality_good, "_", month_str),
    region = make_region(state_code),
    educ_f = factor(
      education,
      levels = c(
        "Elementary school", "Middle school", "High school",
        "Bachelor", "Graduate"
      )
    ),
    gender_f = factor(gender)
  )

hub_state_codes <- c("25", "08", "12", "10", "16")

# ---- OSM driving time to production hubs ------------------------------------
mun_centroids <- municipality_centroids_osm(
  here("Data/Shapefiles/conjunto_de_datos", "00mun.shp")
)
hub_points <- hub_points_from_states(
  here("Data/Shapefiles/conjunto_de_datos", "00ent.shp"),
  hub_state_codes
)

df_hub_time <- tryCatch(
  df |>
    distinct(cvegeo) |>
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
      cvegeo = unique(df$cvegeo),
      travel_time_nearest_hub_min = NA_real_,
      nearest_hub_state = NA_character_,
      ln_travel_time_hub = NA_real_
    )
  }
)

df <- df |>
  left_join(df_hub_time, by = "cvegeo")

has_osm_iv <- any(is.finite(df$ln_travel_time_hub))

# ---- External Hausman instruments ------------------------------------------
qm_totals <- df |>
  group_by(quality_good, month_id) |>
  summarise(qm_n = n(), qm_sum_lnp = sum(ln_price), .groups = "drop")

sqm_totals <- df |>
  group_by(state_code, quality_good, month_id) |>
  summarise(sqm_n = n(), sqm_sum_lnp = sum(ln_price), .groups = "drop")

rqm_totals <- df |>
  group_by(region, quality_good, month_id) |>
  summarise(rqm_n = n(), rqm_sum_lnp = sum(ln_price), .groups = "drop")

df <- df |>
  left_join(qm_totals, by = c("quality_good", "month_id")) |>
  left_join(sqm_totals, by = c("state_code", "quality_good", "month_id")) |>
  left_join(rqm_totals, by = c("region", "quality_good", "month_id")) |>
  mutate(
    z_other_states_qm = if_else(
      qm_n > sqm_n,
      (qm_sum_lnp - sqm_sum_lnp) / (qm_n - sqm_n),
      NA_real_
    ),
    z_other_regions_qm = if_else(
      qm_n > rqm_n,
      (qm_sum_lnp - rqm_sum_lnp) / (qm_n - rqm_n),
      NA_real_
    )
  )

# ---- Seizure-gravity instruments -------------------------------------------
mucd_dir <- here("Data", "MUCD")
has_mucd <- dir.exists(mucd_dir) && length(list.files(mucd_dir, "\\.csv$")) > 0

if (has_mucd) {
  mucd_files <- list.files(mucd_dir, pattern = "\\.csv$", full.names = TRUE)
  seizure_muni_month <- bind_rows(lapply(mucd_files, read_mucd)) |>
    transmute(
      cvegeo_seizure = str_pad(as.character(idMunicipio), 5, pad = "0"),
      seizure_state = str_sub(cvegeo_seizure, 1, 2),
      seizure_region = make_region(seizure_state),
      month_id = as.integer(Anio) * 12L + as.integer(Mes),
      seizure_kg = parse_number(
        as.character(Aseguramientode_marihuana),
        na = c("", "NA")
      ),
      seed_kg = parse_number(
        as.character(Aseguramientode_semillademarihuana),
        na = c("", "NA")
      ),
      crop_destroy = parse_number(
        as.character(Destrucciondecultivos_marihuana),
        na = c("", "NA")
      )
    ) |>
    mutate(
      across(c(seizure_kg, seed_kg, crop_destroy), ~ replace_na(.x, 0)),
      enforcement_total = seizure_kg + seed_kg + 100 * crop_destroy
    ) |>
    filter(!is.na(cvegeo_seizure), nchar(cvegeo_seizure) == 5) |>
    group_by(cvegeo_seizure, seizure_state, seizure_region, month_id) |>
    summarise(
      seizure_kg = sum(seizure_kg, na.rm = TRUE),
      enforcement_total = sum(enforcement_total, na.rm = TRUE),
      .groups = "drop"
    )

  time_sp <- readRDS(here("Data", "Derived", "osrm_seizure_purchase_time_min.rds"))
  time_ps <- t(time_sp)

  purchase_cvegeo <- rownames(time_ps)
  seizure_cvegeo <- colnames(time_ps)
  seizure_state <- str_sub(seizure_cvegeo, 1, 2)
  purchase_region <- make_region(str_sub(purchase_cvegeo, 1, 2))
  seizure_region <- make_region(seizure_state)

  make_weight_matrix <- function(kind = c("hub_only", "outside_region")) {
    kind <- match.arg(kind)
    w <- 1 / (1 + time_ps)
    if (kind == "hub_only") {
      w[, !seizure_state %in% hub_state_codes] <- 0
    }
    if (kind == "outside_region") {
      same_region <- outer(purchase_region, seizure_region, FUN = "==")
      w[same_region] <- 0
    }
    w
  }

  gravity_for_month <- function(weight_mat, target_month_id, value_col) {
    sz <- seizure_muni_month |>
      filter(
        .data$month_id == .env$target_month_id - 1L,
        cvegeo_seizure %in% seizure_cvegeo
      )
    kg <- rep(0, length(seizure_cvegeo))
    names(kg) <- seizure_cvegeo
    if (nrow(sz) > 0) {
      kg[sz$cvegeo_seizure] <- sz[[value_col]]
    }
    as.numeric(weight_mat %*% kg)
  }

  purchase_months <- sort(unique(df$month_id))
  w_hub <- make_weight_matrix("hub_only")
  w_strict <- make_weight_matrix("outside_region")

  gravity_panel <- bind_rows(lapply(purchase_months, function(mm) {
    tibble(
      cvegeo = purchase_cvegeo,
      month_id = mm,
      Z_hub_seizure_lag1 = gravity_for_month(w_hub, mm, "seizure_kg"),
      Z_strict_upstream_lag1 = gravity_for_month(w_strict, mm, "seizure_kg")
    )
  })) |>
    mutate(
      ln_Z_hub_seizure_lag1 = log1p(Z_hub_seizure_lag1),
      ln_Z_strict_upstream_lag1 = log1p(Z_strict_upstream_lag1)
    )

  df <- df |>
    left_join(gravity_panel, by = c("cvegeo", "month_id"))
} else {
  df <- df |>
    mutate(
      ln_Z_hub_seizure_lag1 = NA_real_,
      ln_Z_strict_upstream_lag1 = NA_real_
    )
}

has_hub_gravity <- any(is.finite(df$ln_Z_hub_seizure_lag1))
has_strict_gravity <- any(is.finite(df$ln_Z_strict_upstream_lag1))

# ---- Estimation helpers ------------------------------------------------------
rhs_controls_state <- "quality_good + age + gender_f + educ_f"
rhs_controls_qm <- "age + gender_f + educ_f"
fe_state_month <- "state_month"
fe_quality_month <- "state_month + quality_month"

estimate_ols <- function(data = df, fe = fe_state_month,
                         controls = rhs_controls_state) {
  feols(
    as.formula(paste0("ln_quantity ~ ln_price + ", controls, " | ", fe)),
    data = data,
    cluster = ~state_code
  )
}

estimate_iv <- function(z, data = df, fe = fe_state_month,
                        controls = rhs_controls_state) {
  feols(
    as.formula(paste0("ln_quantity ~ ", controls, " | ", fe, " | ln_price ~ ", z)),
    data = data,
    cluster = ~state_code
  )
}

wu_hausman_p <- function(z, data = df, fe = fe_state_month,
                         controls = rhs_controls_state) {
  data_z <- data |>
    filter(!is.na(.data[[z]]), is.finite(.data[[z]]))
  fs <- feols(
    as.formula(paste0("ln_price ~ ", z, " + ", controls, " | ", fe)),
    data = data_z,
    cluster = ~state_code
  )
  data_z$resid_fs <- residuals(fs)
  wu <- feols(
    as.formula(paste0("ln_quantity ~ ln_price + resid_fs + ", controls, " | ", fe)),
    data = data_z,
    cluster = ~state_code
  )
  t_stat <- safe_coef(wu, "resid_fs") / safe_se(wu, "resid_fs")
  2 * pt(-abs(t_stat), df = wu$nobs - length(coef(wu)))
}

ar_ci <- function(z, data = df, fe = fe_state_month,
                  controls = rhs_controls_state,
                  grid = seq(-2.5, 0.5, by = 0.01)) {
  data_z <- data |>
    filter(!is.na(.data[[z]]), is.finite(.data[[z]]))
  pvals <- sapply(grid, function(b0) {
    data_z$y_tilde <- data_z$ln_quantity - b0 * data_z$ln_price
    model <- feols(
      as.formula(paste0("y_tilde ~ ", z, " + ", controls, " | ", fe)),
      data = data_z,
      cluster = ~state_code
    )
    t_stat <- safe_coef(model, z) / safe_se(model, z)
    2 * pt(-abs(t_stat), df = model$nobs - length(coef(model)))
  })
  keep <- grid[pvals >= 0.05]
  if (length(keep) == 0) return(c(NA_real_, NA_real_))
  range(keep)
}

model_specs <- list()
diagnostics <- list()

add_model <- function(key, label, instrument = NA_character_, model, role,
                      fe = fe_state_month, controls = rhs_controls_state) {
  model_specs[[label]] <<- model
  z <- instrument
  coef_nm <- iv_coef_name(model)
  fstat <- if (!is.na(z)) extract_ivf(model) else NA_real_
  diagnostics[[key]] <<- tibble(
    key = key,
    label = label,
    role = role,
    instrument = z,
    n = nobs(model),
    first_stage_F = fstat,
    estimate = if (!is.na(coef_nm)) unname(coef(model)[coef_nm]) else NA_real_,
    se = if (!is.na(coef_nm)) unname(se(model)[coef_nm]) else NA_real_,
    wu_hausman_p = if (!is.na(z) && is.finite(fstat) && fstat >= 10) {
      wu_hausman_p(z, fe = fe, controls = controls)
    } else {
      NA_real_
    },
    ar_ci_low = if (!is.na(z)) ar_ci(z, fe = fe, controls = controls)[1] else NA_real_,
    ar_ci_high = if (!is.na(z)) ar_ci(z, fe = fe, controls = controls)[2] else NA_real_
  )
}

m_ols <- estimate_ols()
add_model(
  key = "ols",
  label = "(1) OLS",
  model = m_ols,
  role = "Comparable OLS benchmark with state-by-month FE and demographics"
)

if (has_osm_iv) {
  m_osm <- estimate_iv("ln_travel_time_hub")
  add_model(
    key = "osm",
    label = "(2) OSM time",
    instrument = "ln_travel_time_hub",
    model = m_osm,
    role = "Geographic production-hub driving-time diagnostic"
  )
}

if (has_hub_gravity) {
  m_hub_gravity <- estimate_iv("ln_Z_hub_seizure_lag1")
  add_model(
    key = "hub_gravity",
    label = "(3) Hub seizures",
    instrument = "ln_Z_hub_seizure_lag1",
    model = m_hub_gravity,
    role = "Lagged hub-only marijuana-seizure gravity"
  )
}

if (has_strict_gravity) {
  m_strict_gravity <- estimate_iv("ln_Z_strict_upstream_lag1")
  add_model(
    key = "strict_gravity",
    label = "(4) Strict upstream",
    instrument = "ln_Z_strict_upstream_lag1",
    model = m_strict_gravity,
    role = "Lagged outside-region upstream marijuana-seizure gravity"
  )
}

m_haus_state <- estimate_iv(
  "z_other_states_qm",
  fe = fe_quality_month,
  controls = rhs_controls_qm
)
add_model(
  key = "haus_state",
  label = "(5) Hausman: other states",
  instrument = "z_other_states_qm",
  model = m_haus_state,
  role = "External Hausman price instrument; excludes own state",
  fe = fe_quality_month,
  controls = rhs_controls_qm
)

m_haus_region <- estimate_iv(
  "z_other_regions_qm",
  fe = fe_quality_month,
  controls = rhs_controls_qm
)
add_model(
  key = "haus_region",
  label = "(6) Hausman: other regions",
  instrument = "z_other_regions_qm",
  model = m_haus_region,
  role = "External Hausman price instrument; excludes own macroregion",
  fe = fe_quality_month,
  controls = rhs_controls_qm
)

diagnostic_table <- bind_rows(diagnostics)
write_csv(diagnostic_table, here("Data", "Derived", "main_iv_diagnostics.csv"))

# ---- Export combined IV table -----------------------------------------------
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

combined_notes <- c(
  "Dependent variable: ln(quantity in grams).",
  "All columns include state-by-month fixed effects and demographic controls.",
  "OSM time: instrument = ln(OSM/OSRM source-to-market driving time from nearest production-hub state centroid).",
  "Hub seizures: instrument = ln(1 + one-month-lag hub-only marijuana-seizure gravity index).",
  "Strict upstream: instrument = ln(1 + one-month-lag outside-region marijuana-seizure gravity index).",
  "Hausman instruments use other-market mean ln(price) within the same quality-by-month cell, excluding the buyer's own state or macroregion.",
  "The external Hausman columns also include quality-by-month fixed effects, so the standalone quality indicator is absorbed.",
  "For weak IV columns, rely on Anderson--Rubin diagnostics rather than conventional stars.",
  "SEs clustered at state level."
)

combined_path <- file.path(output_dir, "table_combined_iv.tex")
msummary(
  model_specs,
  output = combined_path,
  stars = c('*' = .1, '**' = .05, '***' = .01),
  coef_map = c(
    "ln_price" = "ln(Price)",
    "fit_ln_price" = "ln(Price) [IV]",
    "quality_good" = "Good quality",
    "age" = "Age",
    gender_coef_map,
    educ_coef_map
  ),
  gof_map = c("nobs", "r.squared", "adj.r.squared"),
  title = "OLS and IV estimates: alternative identification strategies",
  notes = combined_notes
)
insert_label_after_caption(combined_path, "tab:iv-combined")
force_table_here(combined_path)
wrap_tabular_in_resizebox(combined_path)

# ---- Export diagnostics table ------------------------------------------------
diag_cols <- diagnostic_table |>
  filter(key != "ols")

diag_header <- paste0(
  "Diagnostic & ",
  paste(diag_cols$label, collapse = " & "),
  "\\\\"
)

diag_line <- function(name, values) {
  paste0(name, " & ", paste(values, collapse = " & "), "\\\\")
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
    paste0("\\begin{tabular}[t]{l", paste(rep("c", nrow(diag_cols)), collapse = ""), "}"),
    "\\toprule",
    diag_header,
    "\\midrule",
    diag_line("First-stage F", fmt_num_vec(diag_cols$first_stage_F, 2)),
    diag_line("IV estimate", mapply(fmt_coef_se, diag_cols$estimate, diag_cols$se)),
    diag_line("Wu--Hausman p-value", fmt_p_vec(diag_cols$wu_hausman_p)),
    diag_line(
      "Anderson--Rubin 95\\% CI",
      mapply(function(lo, hi) fmt_ci(c(lo, hi)), diag_cols$ar_ci_low, diag_cols$ar_ci_high)
    ),
    "\\bottomrule",
    paste0("\\multicolumn{", nrow(diag_cols) + 1L, "}{l}{\\rule{0pt}{1em}All diagnostics correspond to Table~\\ref{tab:iv-combined}.}\\\\"),
    paste0("\\multicolumn{", nrow(diag_cols) + 1L, "}{l}{\\rule{0pt}{1em}Wu--Hausman p-values are reported only for specifications with first-stage F $\\geq$ 10.}\\\\"),
    paste0("\\multicolumn{", nrow(diag_cols) + 1L, "}{l}{\\rule{0pt}{1em}The OSM column is retained to compare with distance-based instruments in prior cannabis studies.}\\\\"),
    paste0("\\multicolumn{", nrow(diag_cols) + 1L, "}{l}{\\rule{0pt}{1em}Hub-only seizure gravity is statistically stronger; the stricter outside-region seizure instrument reports the exclusion-strength tradeoff.}\\\\"),
    "\\end{tabular}",
    "}",
    "\\end{table}"
  ),
  iv_diag_path
)

cat("\nMain IV diagnostics written to:\n")
cat("  ", combined_path, "\n", sep = "")
cat("  ", iv_diag_path, "\n", sep = "")
print(diagnostic_table)
