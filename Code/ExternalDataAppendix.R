###############################################################################
# ExternalDataAppendix.R
# Descriptive diagnostics for MUCD enforcement data and OSM/OSRM driving times.
###############################################################################

project_lib <- file.path(getwd(), ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(here)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

output_dir <- here("LaTeX", "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fmt_int <- function(x) {
  ifelse(is.na(x), "--", formatC(round(x), format = "d", big.mark = ","))
}
fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits, big.mark = ","))
}
fmt_month_year <- function(x) format(x, "%B %Y")

write_tex <- function(path, lines) {
  writeLines(lines, path, useBytes = TRUE)
}

# ---------------------------------------------------------------------------
# MUCD enforcement data
# ---------------------------------------------------------------------------
mucd_dir <- here("Data", "MUCD")
mucd_files <- list.files(mucd_dir, pattern = "\\.csv$", full.names = TRUE)

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
    col_types = cols(.default = col_character()),
    locale = readr::locale(encoding = "latin1"),
    show_col_types = FALSE
  )
}

seizure_raw <- bind_rows(lapply(mucd_files, read_mucd))

seizures <- seizure_raw |>
  transmute(
    cvegeo_seizure = str_pad(as.character(idMunicipio), 5, pad = "0"),
    state_code_s   = str_sub(str_pad(as.character(idMunicipio), 5, pad = "0"), 1, 2),
    anio           = as.integer(Anio),
    mes            = as.integer(Mes),
    month_id       = as.integer(Anio) * 12L + as.integer(Mes),
    month_date     = as.Date(sprintf("%04d-%02d-01", as.integer(Anio), as.integer(Mes))),
    destruction    = readr::parse_number(
      as.character(Destrucciondecultivos_marihuana),
      na = c("", "NA")
    ),
    seizure_kg     = readr::parse_number(
      as.character(Aseguramientode_marihuana),
      na = c("", "NA")
    ),
    seed_quantity  = readr::parse_number(
      as.character(Aseguramientode_semillademarihuana),
      na = c("", "NA")
    )
  ) |>
  mutate(
    seizure_kg = replace_na(seizure_kg, 0)
  )

seizure_muni_month <- seizures |>
  filter(!is.na(cvegeo_seizure), nchar(cvegeo_seizure) == 5) |>
  group_by(cvegeo_seizure, month_id, month_date) |>
  summarise(
    state_code_s = first(state_code_s),
    seizure_kg = sum(seizure_kg, na.rm = TRUE),
    .groups = "drop"
  )

purchase_data <- arrow::read_parquet(
  here("Data", "Marijuana_Prices_in_Mexico_clean.parquet")
)
purchase_dates <- as.Date(purchase_data$purchase_date)
purchase_months <- sort(unique(
  as.integer(format(purchase_dates, "%Y")) * 12L +
    as.integer(format(purchase_dates, "%m"))
))
seizure_lag_months_used <- sort(unique(purchase_months - 1L))

seizure_muni_month_used <- seizure_muni_month |>
  filter(month_id %in% seizure_lag_months_used, seizure_kg > 0)
seizures_used <- seizure_muni_month_used$seizure_kg
seizure_sources_used <- seizure_muni_month_used |>
  distinct(cvegeo_seizure) |>
  pull(cvegeo_seizure)
mucd_months <- seizure_muni_month_used$month_date
mucd_period <- paste0(
  fmt_month_year(min(mucd_months, na.rm = TRUE)),
  "--",
  fmt_month_year(max(mucd_months, na.rm = TRUE))
)

mucd_rows <- c(
  paste0("Seizure lag period used & ", mucd_period, "\\\\"),
  paste0("Lag months used & ",
         fmt_int(n_distinct(seizure_muni_month_used$month_id)), "\\\\"),
  paste0("States with marijuana seizures & ",
         fmt_int(n_distinct(seizure_muni_month_used$state_code_s)), "\\\\"),
  paste0("Municipalities with marijuana seizures & ",
         fmt_int(n_distinct(seizure_muni_month_used$cvegeo_seizure)), "\\\\"),
  paste0("Municipality-month seizure records used & ",
         fmt_int(nrow(seizure_muni_month_used)), "\\\\"),
  paste0("Total marijuana seized in lag months (kg) & ",
         fmt_num(sum(seizure_muni_month_used$seizure_kg), 1), "\\\\"),
  paste0("Median seizure (kg) & ",
         fmt_num(median(seizures_used), 1), "\\\\"),
  paste0("90th percentile seizure (kg) & ",
         fmt_num(unname(quantile(seizures_used, 0.90)), 1), "\\\\"),
  paste0("Maximum seizure (kg) & ",
         fmt_num(max(seizures_used), 1), "\\\\")
)

write_tex(
  file.path(output_dir, "table_appendix_mucd.tex"),
  c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{MUCD marijuana-seizure records used in IV diagnostics}",
    "\\label{tab:appendix-mucd}",
    "\\begin{threeparttable}",
    "\\begin{tabular}{lr}",
    "\\toprule",
    "Statistic & Value\\\\",
    "\\midrule",
    mucd_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]\\small",
    "\\item \\textit{Notes:} MUCD denotes M\\'{e}xico Unido contra la Delincuencia. The table reports the municipality-month marijuana-seizure records in the one-month lag window used to construct the seizure-gravity instrument. Purchases in the estimation sample run from January through September 2019, so the lagged MUCD months used are December 2018 through August 2019. MUCD compiles official responses to public-information requests submitted to federal security and justice institutions.",
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{table}"
  )
)

# ---------------------------------------------------------------------------
# OSM/OSRM driving-time matrices
# ---------------------------------------------------------------------------
hub_cache_path <- here("Data", "Derived", "osrm_statehub_travel_time_min.rds")
hub_mat <- if (file.exists(hub_cache_path)) readRDS(hub_cache_path) else NULL
seizure_mat <- readRDS(here("Data", "Derived", "osrm_seizure_purchase_time_min.rds"))
seizure_mat_positive <- seizure_mat[
  intersect(rownames(seizure_mat), seizure_sources_used),
  ,
  drop = FALSE
]

matrix_stats <- function(mat, label) {
  if (is.null(mat)) {
    return(tibble(
      label = label,
      origins = NA_real_,
      destinations = NA_real_,
      od_pairs = NA_real_,
      mean = NA_real_,
      median = NA_real_,
      p10 = NA_real_,
      p90 = NA_real_,
      max = NA_real_
    ))
  }

  values <- as.numeric(mat)
  tibble(
    label = label,
    origins = nrow(mat),
    destinations = ncol(mat),
    od_pairs = length(values),
    mean = mean(values, na.rm = TRUE),
    median = median(values, na.rm = TRUE),
    p10 = unname(quantile(values, 0.10, na.rm = TRUE)),
    p90 = unname(quantile(values, 0.90, na.rm = TRUE)),
    max = max(values, na.rm = TRUE)
  )
}

travel_stats <- bind_rows(
  matrix_stats(hub_mat, "Production-hub states"),
  matrix_stats(seizure_mat_positive, "Seizure municipalities used")
)

travel_rows <- travel_stats |>
  mutate(
    row = paste0(
      label, " & ",
      fmt_int(origins), " & ",
      fmt_int(destinations), " & ",
      fmt_int(od_pairs), " & ",
      fmt_num(mean, 1), " & ",
      fmt_num(median, 1), " & ",
      fmt_num(p10, 1), " & ",
      fmt_num(p90, 1), " & ",
      fmt_num(max, 1), "\\\\"
    )
  ) |>
  pull(row)

osm_note <- if (is.null(hub_mat)) {
  "\\item \\textit{Notes:} Driving times are one-way source-to-market driving times in minutes computed with the Open Source Routing Machine (OSRM) on OpenStreetMap (OSM) road data. The production-hub matrix links five constructed hub-state points, derived from INEGI state geometries, to purchase municipalities; the seizure-gravity matrix reports MUCD marijuana-seizure municipalities in the one-month lag window linked to purchase municipalities. Dashes indicate that the updated state-centroid hub cache has not been generated in this run. OSRM/OSM is treated as a road-network snapshot queried and cached in April 2026, not as a time series."
} else {
  "\\item \\textit{Notes:} Driving times are one-way source-to-market driving times in minutes computed with the Open Source Routing Machine (OSRM) on OpenStreetMap (OSM) road data. The production-hub matrix links five constructed hub-state points, derived from INEGI state geometries, to purchase municipalities; the seizure-gravity matrix reports MUCD marijuana-seizure municipalities in the one-month lag window linked to purchase municipalities. OSRM/OSM is treated as a road-network snapshot queried and cached in April 2026, not as a time series."
}

write_tex(
  file.path(output_dir, "table_appendix_osm_osrm.tex"),
  c(
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    "\\caption{OSM/OSRM driving-time matrix diagnostics}",
    "\\label{tab:appendix-osm-osrm}",
    "\\begin{threeparttable}",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrrrrrrr}",
    "\\toprule",
    "Matrix & Origins & Dest. & Pairs & Mean & Median & P10 & P90 & Max\\\\",
    "\\midrule",
    travel_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    "\\begin{tablenotes}[flushleft]\\small",
    osm_note,
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{table}"
  )
)

cat("Appendix tables written to:", output_dir, "\n")
