###############################################################################
#  IV_SeizureGravity.R
#  Gravity-style IV using MUCD drug-seizure data
#  -----------------------------------------------------------------------
#  Instrument idea (Bartik/shift-share flavour):
#    For each purchase observation i in municipality m, month t:
#
#    Z_{m,t-1} = sum_j [ kg_marihuana_{j,t-1} / (1 + time(j,m)) ]
#
#    where j indexes municipalities with seizures in month t - 1,
#    kg_marihuana is the quantity seized, and time(j,m) is the OSRM/OpenStreetMap
#    driving time in minutes from seizure municipality j to purchase
#    municipality m.
#
#  Identifying assumptions:
#    RELEVANCE: seizures near municipality m disrupt local supply chains,
#      raising retail prices.  The driving-time decay weighting captures the
#      idea that seizures that are closer in the road network have a stronger
#      supply-disruption effect.
#    EXCLUSION: conditional on state-by-month FE and buyer demographics, the
#      geographic pattern of enforcement actions (driven by federal/military
#      interdiction campaigns) affects buyer quantity only through the
#      price channel.  This is plausible if enforcement targets production
#      and trafficking routes rather than retail demand locations.
#
#  Data source:
#    MUCD — Datos Abiertos sobre Acciones Antidrogas
#    https://datosabiertosdrogas.mucd.org.mx/
#    Monthly CSV files covering December 2018--December 2019, saved in Data/MUCD/
#    MUCD compiles official responses to public-information requests sent to
#    federal security and justice institutions, then cleans and harmonizes the
#    heterogeneous source files into municipality-month data.
#
#  File format (actual):
#    Row 1: "DataDASAACSV" (junk header)
#    Row 2: blank
#    Row 3: CodEntidad, nombreEntidad, nombreMunicipio, idMunicipio,
#            Anio, Mes, Destrucciondecultivos_marihuana,
#            Aseguramientode_marihuana, Aseguramientode_semillademarihuana
#    Data rows contain a trailing empty column.
#    idMunicipio = INEGI 5-digit cvegeo (e.g., "09003")
#    Aseguramientode_marihuana = kg seized (our key variable)
#    Encoding: latin1
###############################################################################

project_lib <- file.path(getwd(), ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(sf)
library(fixest)
library(modelsummary)

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

# ============================================================================
# 1. LOAD AND STACK MONTHLY MUCD FILES
# ============================================================================
mucd_dir <- here("Data", "MUCD")
mucd_files <- list.files(mucd_dir, pattern = "\\.csv$", full.names = TRUE)
cat("Found", length(mucd_files), "MUCD CSV files\n")

read_mucd <- function(filepath) {
  readr::read_csv(
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
    col_types = readr::cols(.default = readr::col_character()),
    locale = readr::locale(encoding = "latin1"),
    show_col_types = FALSE
  )
}

seizure_raw <- bind_rows(lapply(mucd_files, read_mucd))

cat("\nTotal seizure records loaded:", nrow(seizure_raw), "\n")
cat("Columns:", paste(names(seizure_raw), collapse = ", "), "\n")
cat("Year-months present:",
    paste(sort(unique(paste(seizure_raw$Anio, seizure_raw$Mes, sep = "-"))),
          collapse = ", "), "\n")

# ---- Rename to canonical names ----
seizures <- seizure_raw |>
  transmute(
    cvegeo_seizure = str_pad(as.character(idMunicipio), 5, pad = "0"),
    state_code_s   = str_sub(str_pad(as.character(idMunicipio), 5, pad = "0"), 1, 2),
    state_name     = nombreEntidad,
    muni_name      = nombreMunicipio,
    anio           = as.integer(Anio),
    mes            = as.integer(Mes),
    month_id       = as.integer(Anio) * 12L + as.integer(Mes),
    month_date     = as.Date(sprintf("%04d-%02d-01", as.integer(Anio), as.integer(Mes))),
    destruccion_kg = readr::parse_number(
      as.character(Destrucciondecultivos_marihuana),
      na = c("", "NA")
    ),
    aseguramiento_kg = readr::parse_number(
      as.character(Aseguramientode_marihuana),
      na = c("", "NA")
    ),
    semilla_kg = readr::parse_number(
      as.character(Aseguramientode_semillademarihuana),
      na = c("", "NA")
    )
  ) |>
  # Total marijuana enforcement = seizures + crop destruction + seeds
  mutate(
    across(c(destruccion_kg, aseguramiento_kg, semilla_kg), ~ replace_na(.x, 0)),
    total_mj_kg = destruccion_kg + aseguramiento_kg + semilla_kg,
    # Primary IV variable: seizures only (most relevant for supply disruption)
    seizure_kg  = aseguramiento_kg
  )

cat("\nRecords with positive marijuana seizure (aseguramiento):",
    sum(seizures$seizure_kg > 0, na.rm = TRUE), "\n")
cat("Records with positive total enforcement:",
    sum(seizures$total_mj_kg > 0, na.rm = TRUE), "\n")

# Quick summary of seizure magnitudes
cat("\nSeizure kg summary (conditional on > 0):\n")
print(summary(seizures$seizure_kg[seizures$seizure_kg > 0]))

# ---- Aggregate: total kg seized per municipality-month ----
# (Data is already at muni-month level, but some munis may appear
#  in multiple files if there were corrections — aggregate to be safe.)
seizure_muni_month <- seizures |>
  filter(!is.na(cvegeo_seizure), nchar(cvegeo_seizure) == 5) |>
  group_by(cvegeo_seizure, month_id, month_date) |>
  summarise(
    seizure_kg  = sum(seizure_kg, na.rm = TRUE),
    total_mj_kg = sum(total_mj_kg, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nMunicipality-months with data:", nrow(seizure_muni_month), "\n")
cat("  with positive seizure_kg:", sum(seizure_muni_month$seizure_kg > 0), "\n")
cat("Unique seizure municipalities:", n_distinct(seizure_muni_month$cvegeo_seizure), "\n")

# ============================================================================
# 2. COMPUTE MUNICIPALITY CENTROIDS (reuse shapefile)
# ============================================================================
mun_centroids <- municipality_centroids_osm(
  here("Data/Shapefiles/conjunto_de_datos", "00mun.shp")
)

# ============================================================================
# 3. BUILD GRAVITY IV: Z_{m,t}
# ============================================================================
#  Z_{m,t-1} = sum_j [ kg_{j,t-1} / (1 + time(j,m)) ]
#  where time is source-to-market OSRM/OpenStreetMap driving time in minutes.

cat("\nBuilding gravity instrument...\n")

# --- Load purchase data ---
purchase_data <- read_parquet(here("Data", "Marijuana_Prices_in_Mexico_clean.parquet")) |>
  mutate(
    cvegeo    = str_pad(as.character(muni_code), 5, pad = "0"),
    purchase_date = as.Date(purchase_date),
    month_num = month(purchase_date),
    month_id = year(purchase_date) * 12L + month(purchase_date)
  )

purchase_munis  <- unique(purchase_data$cvegeo)
purchase_months <- sort(unique(purchase_data$month_id))

cat("Purchase municipalities:", length(purchase_munis), "\n")
cat("Purchase months:",
    paste(format(as.Date(sprintf("%04d-%02d-01",
                                 (purchase_months - 1L) %/% 12L,
                                 (purchase_months - 1L) %% 12L + 1L)),
                 "%Y-%m"),
          collapse = ", "), "\n")

# --- Centroids ---
cent_purchase <- mun_centroids |> filter(cvegeo %in% purchase_munis)

seizure_cvegeos <- seizure_muni_month |>
  filter(seizure_kg > 0 | total_mj_kg > 0) |>
  distinct(cvegeo_seizure) |>
  pull(cvegeo_seizure)
cent_seizure  <- mun_centroids |> filter(cvegeo %in% seizure_cvegeos)

cat("Purchase municipalities with centroids:", nrow(cent_purchase),
    "of", length(purchase_munis), "\n")
cat("Seizure municipalities with centroids:", nrow(cent_seizure),
    "of", length(seizure_cvegeos), "\n")

# --- Pairwise OSRM driving-time matrix (seizure × purchase), in minutes ---
time_sp <- osrm_duration_matrix(
  src = cent_seizure,
  dst = cent_purchase,
  cache_path = here("Data", "Derived", "osrm_seizure_purchase_time_min.rds"),
  src_id = "cvegeo",
  dst_id = "cvegeo"
)

n_s_raw <- nrow(time_sp)
n_p_raw <- ncol(time_sp)
cat("OSRM source-to-market driving-time matrix:", n_s_raw, "×", n_p_raw, "\n")

# Gravity code expects purchase rows and seizure columns; transpose while
# preserving source-to-market durations.
time_ps <- t(time_sp)

n_p <- nrow(time_ps)
n_s <- ncol(time_ps)

cat("Gravity weighting matrix:", n_p, "×", n_s, "\n")

# --- Compute gravity with one-month lag ---
gravity_iv <- gravity_from_travel_time(
  time_mat = time_ps,
  seizure_muni_month = seizure_muni_month,
  purchase_months = purchase_months,
  kg_col = "seizure_kg",
  out_col = "Z_gravity_lag1",
  lag_months = 1,
  missing_lag = "na",
  month_col = "month_id"
)

gravity_iv_total <- gravity_from_travel_time(
  time_mat = time_ps,
  seizure_muni_month = seizure_muni_month,
  purchase_months = purchase_months,
  kg_col = "total_mj_kg",
  out_col = "Z_gravity_total_lag1",
  lag_months = 1,
  missing_lag = "na",
  month_col = "month_id"
)

gravity_iv <- gravity_iv |>
  left_join(gravity_iv_total, by = c("cvegeo", "month_id"))

cat("\nOne-month-lag gravity IV computed:",
    nrow(gravity_iv), "municipality-month cells\n")
cat("Z_gravity_lag1 summary:\n")
print(summary(gravity_iv$Z_gravity_lag1))
cat("Months without prior-month MUCD data:",
    paste(gravity_iv |>
            filter(is.na(Z_gravity_lag1)) |>
            distinct(month_id) |>
            mutate(month_label = format(as.Date(sprintf("%04d-%02d-01",
                                                        (month_id - 1L) %/% 12L,
                                                        (month_id - 1L) %% 12L + 1L)),
                                        "%Y-%m")) |>
            pull(month_label),
          collapse = ", "), "\n")

# ============================================================================
# 4. MERGE WITH PURCHASE DATA AND RUN REGRESSIONS
# ============================================================================

# State-month aggregate seizure, lagged one month (simple alternative IV)
seizure_state_month <- seizure_muni_month |>
  mutate(state_code_s = str_sub(cvegeo_seizure, 1, 2)) |>
  group_by(state_code_s, month_id) |>
  summarise(state_seizure_kg_lag1 = sum(seizure_kg, na.rm = TRUE), .groups = "drop") |>
  mutate(month_id = month_id + 1L) |>
  select(state_code_s, month_id, state_seizure_kg_lag1)

# Merge
df <- purchase_data |>
  left_join(gravity_iv, by = c("cvegeo", "month_id")) |>
  mutate(
    state_code   = str_sub(cvegeo, 1, 2),
    ln_quantity  = log(quantity_g),
    ln_price     = log(price_mxn_g),
    quality_good = as.integer(quality == "good"),
    month_str    = format(purchase_date, "%Y-%m"),
    state_month  = paste0(str_sub(cvegeo, 1, 2), "_", format(purchase_date, "%Y-%m")),
    educ_f       = factor(education,
                          levels = c("Elementary school", "Middle school",
                                     "High school", "Bachelor", "Graduate")),
    gender_f     = factor(gender)
  ) |>
  left_join(seizure_state_month,
            by = c("state_code" = "state_code_s", "month_id" = "month_id"))

# Construct logged instruments; missing prior-month lags remain NA by design.
df <- df |>
  mutate(
    ln_Z_gravity_lag1    = log(1 + Z_gravity_lag1),
    ln_Z_grav_total_lag1 = log(1 + Z_gravity_total_lag1),
    ln_state_kg_lag1     = log(1 + state_seizure_kg_lag1)
  )

# Trimming
p01 <- quantile(df$price_mxn_g, 0.01)
p99 <- quantile(df$price_mxn_g, 0.99)
q01 <- quantile(df$quantity_g, 0.01)
q99 <- quantile(df$quantity_g, 0.99)

df <- df |>
  mutate(in_trim = price_mxn_g >= p01 & price_mxn_g <= p99 &
                   quantity_g  >= q01 & quantity_g  <= q99)

cat("\n========= MERGE DIAGNOSTICS =========\n")
cat("Merged sample:", nrow(df), "observations\n")
cat("Non-missing lagged gravity IV (seizures only):",
    sum(!is.na(df$Z_gravity_lag1)), "\n")
cat("Non-zero lagged gravity IV (seizures only):",
    sum(df$Z_gravity_lag1 > 0, na.rm = TRUE), "\n")
cat("Non-zero lagged gravity IV (total enforcement):",
    sum(df$Z_gravity_total_lag1 > 0, na.rm = TRUE), "\n")
cat("Non-zero lagged state seizure kg:",
    sum(df$state_seizure_kg_lag1 > 0, na.rm = TRUE), "\n")
cat("\nln_Z_gravity_lag1 summary:\n")
print(summary(df$ln_Z_gravity_lag1))

# ============================================================================
# 5. REGRESSIONS
# ============================================================================
cat("\n\n========= REGRESSIONS =========\n")

# ---- OLS baseline (state-by-month FE for comparability with m6) ----
ols_base <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
                  | state_month,
                  data = df, cluster = ~state_code)

# ---- IV 1: gravity (seizures only), state-by-month FE ----
iv_grav <- feols(
  ln_quantity ~ quality_good + age + gender_f + educ_f
  | state_month
  | ln_price ~ ln_Z_gravity_lag1,
  data = df, cluster = ~state_code
)

# ---- IV 2: gravity (total enforcement), state-by-month FE ----
iv_grav_total <- feols(
  ln_quantity ~ quality_good + age + gender_f + educ_f
  | state_month
  | ln_price ~ ln_Z_grav_total_lag1,
  data = df, cluster = ~state_code
)

# ---- IV 3: overidentified (gravity + state-level seizures) ----
iv_overid <- feols(
  ln_quantity ~ quality_good + age + gender_f + educ_f
  | state_month
  | ln_price ~ ln_Z_gravity_lag1 + ln_state_kg_lag1,
  data = df, cluster = ~state_code
)

# ---- IV 4: gravity, trimmed sample ----
iv_grav_trim <- feols(
  ln_quantity ~ quality_good + age + gender_f + educ_f
  | state_month
  | ln_price ~ ln_Z_gravity_lag1,
  data = filter(df, in_trim), cluster = ~state_code
)

# ---- IV 5-6: gravity by quality segment ----
iv_grav_good <- feols(
  ln_quantity ~ age + gender_f + educ_f
  | state_month
  | ln_price ~ ln_Z_gravity_lag1,
  data = filter(df, quality_good == 1), cluster = ~state_code
)

iv_grav_bad <- feols(
  ln_quantity ~ age + gender_f + educ_f
  | state_month
  | ln_price ~ ln_Z_gravity_lag1,
  data = filter(df, quality_good == 0), cluster = ~state_code
)

# ============================================================================
# 6. FIRST-STAGE DIAGNOSTICS
# ============================================================================
cat("\n\n========= FIRST-STAGE DIAGNOSTICS =========\n")

# Explicit first stage (with state-by-month FE)
fs_grav <- feols(ln_price ~ ln_Z_gravity_lag1 + quality_good + age + gender_f + educ_f
                 | state_month,
                 data = df, cluster = ~state_code)

cat("\n--- First stage: ln(price) ~ ln(1 + Z_gravity_lag1) | state x month FE ---\n")
cat("  Coef on ln_Z_gravity_lag1:",
    round(coef(fs_grav)["ln_Z_gravity_lag1"], 5), "\n")
cat("  SE:", round(se(fs_grav)["ln_Z_gravity_lag1"], 5), "\n")
cat("  t-stat:", round(coef(fs_grav)["ln_Z_gravity_lag1"] /
                       se(fs_grav)["ln_Z_gravity_lag1"], 2), "\n")

cat("\n  IV F-stat (seizures only, state-by-month FE):\n  ")
print(fitstat(iv_grav, "ivf"))

cat("\n  IV F-stat (total enforcement, state-by-month FE):\n  ")
print(fitstat(iv_grav_total, "ivf"))

cat("\n  IV F-stat (overidentified):\n  ")
print(fitstat(iv_overid, "ivf"))

# Sargan/Hansen overid test
cat("\n  Sargan test (overidentified):\n  ")
print(fitstat(iv_overid, "sargan"))

# ============================================================================
# 7. EXPORT TABLES
# ============================================================================
output_dir <- here("LaTeX", "tables")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- Table 5: Main IV comparison ----
tab5 <- list(
  "(1) OLS"              = ols_base,
  "(2) IV seizures"      = iv_grav,
  "(3) IV total enf."    = iv_grav_total,
  "(4) IV overid."       = iv_overid,
  "(5) IV trimmed"       = iv_grav_trim
)

cat("\n\n========== TABLE 5: IV WITH SEIZURE GRAVITY ==========\n")
msummary(tab5, output = "markdown",
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("ln_price" = "ln(Price)",
                      "fit_ln_price" = "ln(Price) [IV]",
                      "quality_good" = "Good quality",
                      "age" = "Age",
                      gender_coef_map,
                      educ_coef_map),
         gof_map = c("nobs", "r.squared", "adj.r.squared"))

msummary(tab5,
         output = file.path(output_dir, "table5_iv_seizure_gravity.tex"),
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("ln_price" = "ln(Price)",
                      "fit_ln_price" = "ln(Price) [IV]",
                      "quality_good" = "Good quality",
                      "age" = "Age",
                      gender_coef_map,
                      educ_coef_map),
         gof_map = c("nobs", "r.squared", "adj.r.squared"),
         title = "IV estimates: seizure-gravity instrument",
         notes = c(
           "Dependent variable: ln(quantity in grams).",
           "Instrument: ln(1 + Zm,t-1), where Zm,t-1 is the OSM driving-time-weighted sum of marijuana seizures in source municipalities in the prior month.",
           "Seizure data: MUCD (datosabiertosdrogas.mucd.org.mx), compiled from official transparency-request responses.",
           "The estimation sample uses one-month-lagged MUCD marijuana seizures from December 2018--August 2019.",
           "All columns include state-by-month FE and demographics; SEs clustered at state level.",
           "* p $<$ 0.1, ** p $<$ 0.05, *** p $<$ 0.01."
         ))
insert_label_after_caption(
  file.path(output_dir, "table5_iv_seizure_gravity.tex"),
  "tab:iv-gravity"
)
force_table_here(file.path(output_dir, "table5_iv_seizure_gravity.tex"))

# ---- Table 6: Quality-segmented IV ----
tab6 <- list(
  "(2) IV full"    = iv_grav,
  "(5) IV good"    = iv_grav_good,
  "(6) IV bad"     = iv_grav_bad,
  "(7) IV trimmed" = iv_grav_trim
)

cat("\n\n========== TABLE 6: IV BY QUALITY SEGMENT ==========\n")
msummary(tab6, output = "markdown",
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("fit_ln_price" = "ln(Price) [IV]",
                      "quality_good" = "Good quality",
                      "age" = "Age",
                      gender_coef_map,
                      educ_coef_map),
         gof_map = c("nobs", "r.squared", "adj.r.squared"))

msummary(tab6,
         output = file.path(output_dir, "table6_iv_seizure_quality.tex"),
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("fit_ln_price" = "ln(Price) [IV]",
                      "quality_good" = "Good quality",
                      "age" = "Age",
                      gender_coef_map,
                      educ_coef_map),
         gof_map = c("nobs", "r.squared", "adj.r.squared"),
         title = "IV estimates by quality segment (seizure-gravity instrument)",
         notes = c(
           "Dependent variable: ln(quantity in grams).",
           "Instrument: ln(1 + Zm,t-1), the one-month-lag OSM driving-time-weighted seizure-gravity index based on marijuana seizures only.",
           "The estimation sample uses one-month-lagged MUCD marijuana seizures from December 2018--August 2019.",
           "All columns include state-by-month FE and demographics. SEs clustered at state level.",
           "* p $<$ 0.1, ** p $<$ 0.05, *** p $<$ 0.01."
         ))
insert_label_after_caption(
  file.path(output_dir, "table6_iv_seizure_quality.tex"),
  "tab:iv-gravity-quality"
)
force_table_here(file.path(output_dir, "table6_iv_seizure_quality.tex"))

cat("\n\nAll tables saved to:", output_dir, "\n")
cat("Done.\n")
