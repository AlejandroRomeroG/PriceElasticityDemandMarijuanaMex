###############################################################################
#  RegressionAnalysis.R
#  Price Elasticity of Demand for Marijuana in Mexico
#  Author: Alejandro Romero González
#  -----------------------------------------------------------------------
#  This script estimates log–log demand regressions of the form:
#     ln(quantity_g) = alpha + beta * ln(price_mxn_g) + controls + FE + epsilon
#  where beta is the price elasticity of demand.
#
#  Specifications:
#    (1) Pooled OLS (baseline)
#    (2) + quality control
#    (3) + demographics (age, gender, education)
#    (4) + state fixed effects
#    (5) + month fixed effects
#    (6) + state × month FE
#    (7) Quality-segmented models (good vs bad separately)
#    (8) Price × quality interaction
#    (9) Robustness: trimmed sample (outlier removal)
#    (10) Robustness: heaping correction (round quantities only)
#    (11) IV estimation: OSM driving time to production hubs
#
#  All standard errors are heteroskedasticity-robust (HC1) unless
#  clustering is applied (clustered at state level where FE allow).
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
library(fixest)       # Fast fixed-effects and IV estimation
library(modelsummary)  # Publication-quality regression tables
library(lubridate)
library(stringr)
library(sf)

fixest::setFixest_estimation(fixef.rm = "infinite_coef")
options(
  modelsummary_factory_latex = "kableExtra",
  modelsummary_format_numeric_latex = "plain"
)
source(here("Code", "OSMTravelTime.R"))

# ---- 1. Load and prepare data --------------------------------------------
df <- read_parquet(here("Data", "Marijuana_Prices_in_Mexico_clean.parquet"))

# Construct log variables and time identifiers
df <- df |>
  mutate(
    ln_quantity  = log(quantity_g),
    ln_price     = log(price_mxn_g),
    quality_good = as.integer(quality == "good"),

    # State code (first 2 digits of 5-digit municipality code)
    state_code = str_pad(as.character(muni_code), width = 5, pad = "0") |>
                 str_sub(1, 2),

    # Time identifiers
    purchase_date = as.Date(purchase_date),
    month         = floor_date(purchase_date, "month"),
    month_str     = format(month, "%Y-%m"),

    # Education as unordered factor so controls enter as dummies, not polynomials
    educ_f = factor(education,
                    levels = c("Elementary school", "Middle school",
                               "High school", "Bachelor", "Graduate")),

    # Gender factor
    gender_f = factor(gender)
  )

cat("Analytic sample:", nrow(df), "observations\n")
cat("States:", n_distinct(df$state_code), " | Municipalities:", n_distinct(df$muni_code), "\n")
cat("Months:", n_distinct(df$month_str), "\n")
cat("Quality split: good =", sum(df$quality_good), "bad =", sum(!df$quality_good), "\n")

# ---- 2. Flag outliers for robustness checks ------------------------------
#  Winsorize at 1st/99th percentile of price_mxn_g and quantity_g
p01 <- quantile(df$price_mxn_g, 0.01)
p99 <- quantile(df$price_mxn_g, 0.99)
q01 <- quantile(df$quantity_g, 0.01)
q99 <- quantile(df$quantity_g, 0.99)

df <- df |>
  mutate(
    in_trim = price_mxn_g >= p01 & price_mxn_g <= p99 &
              quantity_g  >= q01 & quantity_g  <= q99,
    # Flag heaped (round) quantities: multiples of 5 or 28 (ounce)
    quantity_round = (quantity_g %% 5 == 0) | (abs(quantity_g - 28) < 0.5) |
                     (abs(quantity_g - 56) < 0.5)
  )

cat("Trimmed sample:", sum(df$in_trim), "observations\n")

# ---- 3. Compute OSM driving-time IV candidates ----------------------------
#  Key idea (à la Davis et al. 2015, Halcoussis et al. 2017):
#  Network driving time from major marijuana production / trafficking hubs to
#  the purchase municipality proxies for supply-side transportation costs.
#  Major production states: Sinaloa, Chihuahua, Guerrero, Durango, Michoacán
#  (the "Golden Triangle" + Guerrero + Michoacán).
#
#  We compute hub-to-municipality driving times with OSRM/OpenStreetMap and
#  cache the resulting matrix in Data/Derived. No straight-line fallback is used.

mun_centroids <- municipality_centroids_osm(
  here("Data/Shapefiles/conjunto_de_datos", "00mun.shp")
)

# Production-hub states:
# Sinaloa-25, Chihuahua-08, Guerrero-12, Durango-10, Michoacán-16.
# We construct one hub point per state by calculating the centroid of the
# projected full-state geometry from INEGI polygons.
hub_state_codes <- c("25", "08", "12", "10", "16")

hub_points <- hub_points_from_states(
  here("Data/Shapefiles/conjunto_de_datos", "00ent.shp"),
  hub_state_codes
)

df_cvegeo <- df |>
  mutate(cvegeo = str_pad(as.character(muni_code), width = 5, pad = "0")) |>
  distinct(cvegeo)

df_hub_time <- tryCatch(
  df_cvegeo |>
    left_join(mun_centroids, by = "cvegeo") |>
    nearest_hub_travel_time(
      hub_points = hub_points,
      cache_path = here("Data", "Derived", "osrm_statehub_travel_time_min.rds")
    ),
  error = function(e) {
    warning(
      "OSM driving-time IV unavailable; continuing without Table 4. ",
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

# Merge back
df <- df |>
  mutate(cvegeo = str_pad(as.character(muni_code), width = 5, pad = "0")) |>
  left_join(df_hub_time, by = "cvegeo")

has_osm_iv <- any(is.finite(df$ln_travel_time_hub))
if (has_osm_iv) {
  cat("OSM driving-time IV: median =",
      median(df$travel_time_nearest_hub_min, na.rm = TRUE),
      "minutes, range =",
      range(df$travel_time_nearest_hub_min, na.rm = TRUE), "\n")
}

# ---- 4. Baseline OLS specifications (fixest) ----------------------------
#  setFixest_fml: common controls shorthand
#  All SEs HC1 (robust) by default; cluster at state_code where state FE used.

# (1) Pooled OLS
m1 <- feols(ln_quantity ~ ln_price, data = df, vcov = "hetero")

# (2) + quality dummy
m2 <- feols(ln_quantity ~ ln_price + quality_good, data = df, vcov = "hetero")

# (3) + demographics
m3 <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f,
            data = df, vcov = "hetero")

# (4) + state FE (clustered at state)
m4 <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
            | state_code,
            data = df, cluster = ~state_code)

# (5) + month FE
m5 <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
            | state_code + month_str,
            data = df, cluster = ~state_code)

# (6) State × month FE (absorbs all state-time demand shocks)
#     NOTE: only feasible where cells are not singletons
df <- df |> mutate(state_month = paste0(state_code, "_", month_str))

m6 <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
            | state_month,
            data = df, cluster = ~state_code)

# ---- 5. Quality-segmented models ----------------------------------------
m7_good <- feols(ln_quantity ~ ln_price + age + gender_f + educ_f
                 | state_month,
                 data = filter(df, quality_good == 1),
                 cluster = ~state_code)

m7_bad  <- feols(ln_quantity ~ ln_price + age + gender_f + educ_f
                 | state_month,
                 data = filter(df, quality_good == 0),
                 cluster = ~state_code)

# (8) Interaction model (pooled)
m8 <- feols(ln_quantity ~ ln_price * quality_good + age + gender_f + educ_f
            | state_month,
            data = df, cluster = ~state_code)

# ---- 6. Robustness: trimmed sample & heaping ----------------------------
# (9) Trimmed (drop 1st/99th percentile outliers in price and quantity)
m9 <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
            | state_month,
            data = filter(df, in_trim),
            cluster = ~state_code)

# (10) Excluding heaped quantities (keep only non-round quantities)
m10 <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
             | state_month,
             data = filter(df, !quantity_round),
             cluster = ~state_code)

# ---- 7. Instrumental Variables -------------------------------------------
#  Instrument: ln(OSM driving time from nearest production hub)
#  Identifying assumptions:
#    Relevance: municipalities with longer driving times from production hubs
#      face higher transportation costs → higher retail prices.
#    Exclusion: driving time to production regions affects quantity demanded
#      only through its effect on price, conditional on state-by-month FE and controls.
#      This requires that routing time is not correlated with unobserved demand
#      factors (e.g., local preferences) after conditioning on state FE.
#
#  NOTE: Because driving time varies only at the municipality level and state-by-
#  month FE absorb all state-time shocks, the remaining variation is within-
#  state routing-time differences across municipalities.

if (has_osm_iv) {
  # (11a) IV: state-by-month FE + demographics (comparable to m6)
  m11a <- feols(ln_quantity ~ quality_good + age + gender_f + educ_f
                | state_month
                | ln_price ~ ln_travel_time_hub,
                data = df, cluster = ~state_code)

  # (11b) IV: trimmed sample
  m11b <- feols(ln_quantity ~ quality_good + age + gender_f + educ_f
                | state_month
                | ln_price ~ ln_travel_time_hub,
                data = filter(df, in_trim), cluster = ~state_code)

  # (11c) IV: quality-segmented (good only)
  m11c <- feols(ln_quantity ~ age + gender_f + educ_f
                | state_month
                | ln_price ~ ln_travel_time_hub,
                data = filter(df, quality_good == 1), cluster = ~state_code)

  # ---- 8. First-stage diagnostics ---------------------------------------
  fs11a <- feols(ln_price ~ ln_travel_time_hub + quality_good + age + gender_f + educ_f
                 | state_month,
                 data = df, cluster = ~state_code)

  cat("\n=== FIRST STAGE (OSM driving-time IV) ===\n")
  cat("Coefficient on ln_travel_time_hub:",
      coef(fs11a)["ln_travel_time_hub"], "\n")
  cat("t-stat:",
      coef(fs11a)["ln_travel_time_hub"] / se(fs11a)["ln_travel_time_hub"],
      "\n")
  cat("First-stage F (from IV model):\n")
  print(fitstat(m11a, "ivf"))
} else {
  cat("\n=== OSM driving-time IV skipped: no cached state-hub OSRM matrix available ===\n")
}

# ---- 9. Export regression tables -----------------------------------------
# Table 1: Main OLS specifications
tab1 <- list(
  "(1) Pooled"    = m1,
  "(2) +Quality"  = m2,
  "(3) +Demog."   = m3,
  "(4) +State FE" = m4,
  "(5) +Month FE" = m5,
  "(6) St×Mo FE"  = m6
)

# Table 2: Quality-segmented and interactions
tab2 <- list(
  "(7a) Good only" = m7_good,
  "(7b) Bad only"  = m7_bad,
  "(8) Interaction" = m8
)

# Table 3: Robustness
tab3 <- list(
  "(6) Benchmark" = m6,
  "(9) Trimmed"   = m9,
  "(10) No heap"  = m10
)

# Table 4: IV
if (has_osm_iv) {
  tab4 <- list(
    "(6) OLS benchmark" = m6,
    "(11a) IV full"   = m11a,
    "(11b) IV trim"   = m11b,
    "(11c) IV good"   = m11c
  )
}

# Print tables to console (LaTeX output)
cat("\n\n========== TABLE 1: MAIN OLS SPECIFICATIONS ==========\n")
msummary(tab1, output = "markdown",
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("ln_price" = "ln(Price)",
                      "quality_good" = "Good quality",
                      "age" = "Age"),
         gof_map = c("nobs", "r.squared", "adj.r.squared",
                      "FE: state_code", "FE: month_str", "FE: state_month"))

cat("\n\n========== TABLE 2: QUALITY-SEGMENTED MODELS ==========\n")
msummary(tab2, output = "markdown",
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("ln_price" = "ln(Price)",
                      "quality_good" = "Good quality",
                      "ln_price:quality_good" = "ln(Price) × Good",
                      "age" = "Age"),
         gof_map = c("nobs", "r.squared", "adj.r.squared"))

cat("\n\n========== TABLE 3: ROBUSTNESS ==========\n")
msummary(tab3, output = "markdown",
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("ln_price" = "ln(Price)",
                      "quality_good" = "Good quality",
                      "age" = "Age"),
         gof_map = c("nobs", "r.squared", "adj.r.squared"))

if (has_osm_iv) {
  cat("\n\n========== TABLE 4: INSTRUMENTAL VARIABLES ==========\n")
  msummary(tab4, output = "markdown",
           stars = c('*' = .1, '**' = .05, '***' = .01),
           coef_map = c("ln_price" = "ln(Price)",
                        "fit_ln_price" = "ln(Price) [IV]",
                        "quality_good" = "Good quality",
                        "age" = "Age"),
           gof_map = c("nobs", "r.squared", "adj.r.squared"))
}

# ---- 10. Save LaTeX tables for the paper ---------------------------------
output_dir <- here("LaTeX", "tables")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

tex_num <- function(x, digits = 3) {
  if (is.na(x) || !is.finite(x)) {
    return("NA")
  }
  paste0("\\num{", formatC(x, format = "f", digits = digits), "}")
}

tex_ci <- function(lower, upper, digits = 3) {
  if (any(!is.finite(c(lower, upper)))) {
    return("NA")
  }
  paste0("[", tex_num(lower, digits), ", ", tex_num(upper, digits), "]")
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

# Table 1
table1_path <- file.path(output_dir, "table1_main_ols.tex")
msummary(tab1,
         output = table1_path,
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("ln_price" = "ln(Price)",
                      "quality_good" = "Good quality",
                      "age" = "Age",
                      gender_coef_map,
                      educ_coef_map),
         gof_map = c("nobs", "r.squared", "adj.r.squared",
                      "FE: state_code", "FE: month_str", "FE: state_month"),
         title = "OLS estimates of the price elasticity of demand for marijuana",
         notes = c("Dependent variable: ln(quantity in grams).",
                   "Robust SEs in parentheses (cols 1--3); clustered at state level (cols 4--6).",
                   "* p < 0.1, ** p < 0.05, *** p < 0.01."))
insert_label_after_caption(table1_path, "tab:main-ols")
force_table_here(table1_path)
wrap_tabular_in_resizebox(table1_path)

# Table 2
table2_path <- file.path(output_dir, "table2_quality_segmented.tex")
msummary(tab2,
         output = table2_path,
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("ln_price" = "ln(Price)",
                      "quality_good" = "Good quality",
                      "ln_price:quality_good" = "ln(Price) x Good",
                      "age" = "Age",
                      gender_coef_map,
                      educ_coef_map),
         gof_map = c("nobs", "r.squared", "adj.r.squared"),
         title = "Quality-segmented demand estimates",
         notes = c("Dependent variable: ln(quantity in grams).",
                   "All columns include state-by-month FE and demographics. SEs clustered at state level.",
                   "* p < 0.1, ** p < 0.05, *** p < 0.01."))
insert_label_after_caption(table2_path, "tab:quality-segmented")
force_table_here(table2_path)

# Table 3
table3_path <- file.path(output_dir, "table3_robustness.tex")
msummary(tab3,
         output = table3_path,
         stars = c('*' = .1, '**' = .05, '***' = .01),
         coef_map = c("ln_price" = "ln(Price)",
                      "quality_good" = "Good quality",
                      "age" = "Age",
                      gender_coef_map,
                      educ_coef_map),
         gof_map = c("nobs", "r.squared", "adj.r.squared"),
	         title = "Robustness: outlier trimming and heaping",
	         notes = c("Dependent variable: ln(quantity in grams).",
	                   "All columns include state-by-month FE plus demographics. SEs clustered at state level.",
	                   "Trimmed = drop 1st/99th percentile of price and quantity.",
                   "No heap = exclude round-number quantities (multiples of 5 or near 28g).",
                   "* p < 0.1, ** p < 0.05, *** p < 0.01."))
insert_label_after_caption(table3_path, "tab:robustness")
force_table_here(table3_path)

table4_path <- file.path(output_dir, "table4_iv.tex")
if (has_osm_iv) {
  msummary(tab4,
           output = table4_path,
           stars = c('*' = .1, '**' = .05, '***' = .01),
           coef_map = c("ln_price" = "ln(Price)",
                        "fit_ln_price" = "ln(Price) [IV]",
                        "quality_good" = "Good quality",
                        "age" = "Age",
                        gender_coef_map,
                        educ_coef_map),
           gof_map = c("nobs", "r.squared", "adj.r.squared"),
           title = "Instrumental-variables estimates (OSM driving time to production hubs)",
           notes = c("Dependent variable: ln(quantity in grams).",
                     "Instrument: ln(OSM source-to-market driving time in minutes from nearest",
                     "constructed production-hub point to purchase municipality: Sinaloa, Chihuahua, Guerrero, Durango, Michoacan).",
                     "All columns include state-by-month FE and demographics. SEs clustered at state level.",
                     "For weak IV specifications, inference should rely on weak-instrument diagnostics rather than conventional significance stars.",
                     "* p < 0.1, ** p < 0.05, *** p < 0.01."))
  insert_label_after_caption(table4_path, "tab:iv-distance")
  force_table_here(table4_path)
} else {
  writeLines(
    c(
      "\\begin{table}[H]",
      "\\centering",
      "\\caption{Instrumental-variables estimates (OSM driving time to production hubs)}",
      "\\label{tab:iv-distance}",
      "\\centering",
      "\\begin{tabular}[t]{p{0.86\\linewidth}}",
      "\\toprule",
      "The OSM production-hub driving-time IV is omitted in this run because the updated full-state-centroid hub construction requires a regenerated state-hub OSRM cache.\\\\",
      "\\bottomrule",
      "\\end{tabular}",
      "\\end{table}"
    ),
    table4_path
  )
}

cat("\n\nAll tables saved to:", output_dir, "\n")
cat("Done.\n")
