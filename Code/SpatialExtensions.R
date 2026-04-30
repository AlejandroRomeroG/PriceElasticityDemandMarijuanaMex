## SpatialExtensions.R — Spatial robustness checks for marijuana demand elasticity
## Author: Claude  |  Date: 2026-03-18

project_lib <- file.path(getwd(), ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(here)
library(arrow)
library(dplyr)
library(stringr)
library(lubridate)
library(fixest)
library(sf)

fixest::setFixest_estimation(fixef.rm = "infinite_coef")

# ── Data prep ──────────────────────────────────────────────────────────────────
df <- read_parquet(here("Data", "Marijuana_Prices_in_Mexico_clean.parquet"))

df <- df |> mutate(
  ln_quantity  = log(quantity_g),
  ln_price     = log(price_mxn_g),
  quality_good = as.integer(quality == "good"),
  cvegeo       = str_pad(as.character(muni_code), 5, pad = "0"),
  state_code   = str_sub(cvegeo, 1, 2),
  purchase_date = as.Date(purchase_date),
  month_str    = format(floor_date(purchase_date, "month"), "%Y-%m"),
  educ_f       = factor(education,
                        levels = c("Elementary school", "Middle school",
                                   "High school", "Bachelor", "Graduate")),
  gender_f     = factor(gender)
)

cat("Observations:", nrow(df), "\n")
cat("States:", n_distinct(df$state_code), "\n\n")

# ── Baseline for reference ─────────────────────────────────────────────────────
baseline <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
                  | state_code, data = df, cluster = ~state_code)
cat("=== BASELINE (state FE, clustered at state) ===\n")
cat(sprintf("  ln_price coef = %.4f  (SE = %.4f)\n",
            coef(baseline)["ln_price"], se(baseline)["ln_price"]))
cat(sprintf("  N = %d,  R² = %.4f\n\n", baseline$nobs, r2(baseline, "r2")))

# ══════════════════════════════════════════════════════════════════════════════
# 1. LEAVE-ONE-STATE-OUT JACKKNIFE
# ══════════════════════════════════════════════════════════════════════════════
cat("=== 1. LEAVE-ONE-STATE-OUT JACKKNIFE ===\n")

states <- sort(unique(df$state_code))
jack_elast <- numeric(length(states))
names(jack_elast) <- states

for (i in seq_along(states)) {
  sub <- df |> filter(state_code != states[i])
  m <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
             | state_code, data = sub, cluster = ~state_code)
  jack_elast[i] <- coef(m)["ln_price"]
}

cat(sprintf("  Range of elasticities: [%.4f, %.4f]\n",
            min(jack_elast), max(jack_elast)))
cat(sprintf("  SD across jackknife replications: %.4f\n", sd(jack_elast)))
cat(sprintf("  Mean: %.4f\n", mean(jack_elast)))

# Which states move the estimate most?
dev <- jack_elast - mean(jack_elast)
top3 <- sort(abs(dev), decreasing = TRUE)[1:3]
cat("  Most influential states (by |deviation from mean|):\n")
for (nm in names(top3)) {
  cat(sprintf("    State %s: elasticity = %.4f  (dev = %+.4f)\n",
              nm, jack_elast[nm], dev[nm]))
}
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. SPATIAL LAG OF PRICE (state-month leave-one-out mean)
# ══════════════════════════════════════════════════════════════════════════════
cat("=== 2. SPATIAL LAG OF PRICE (state × month cell mean, excl. self) ===\n")

df <- df |>
  group_by(state_code, month_str) |>
  mutate(
    cell_sum    = sum(ln_price),
    cell_n      = n(),
    ln_price_lag = (cell_sum - ln_price) / (cell_n - 1)
  ) |>
  ungroup()

# Drop singletons (cell_n == 1 → no peer mean)
n_single <- sum(df$cell_n == 1)
cat(sprintf("  Singleton state-month cells dropped: %d obs\n", n_single))
df_lag <- df |> filter(cell_n > 1)

m_lag <- feols(ln_quantity ~ ln_price + ln_price_lag + quality_good + age +
                 gender_f + educ_f | state_code,
               data = df_lag, cluster = ~state_code)

cat(sprintf("  Own-price elasticity: %.4f  (SE = %.4f)\n",
            coef(m_lag)["ln_price"], se(m_lag)["ln_price"]))
cat(sprintf("  Spatial lag coef:     %.4f  (SE = %.4f)\n",
            coef(m_lag)["ln_price_lag"], se(m_lag)["ln_price_lag"]))
cat(sprintf("  N = %d,  R² = %.4f\n", m_lag$nobs, r2(m_lag, "r2")))

# Compare to baseline on same sample
m_lag0 <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
                | state_code, data = df_lag, cluster = ~state_code)
cat(sprintf("  Baseline on same sample: %.4f  (SE = %.4f)\n",
            coef(m_lag0)["ln_price"], se(m_lag0)["ln_price"]))
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# 3. DISTANCE BANDS TO PRODUCTION HUBS
# ══════════════════════════════════════════════════════════════════════════════
cat("=== 3. DISTANCE BANDS TO PRODUCTION HUBS ===\n")

# Read municipality shapefile
mun_sf <- st_read(here("Data", "Shapefiles", "conjunto_de_datos", "00mun.shp"),
                  quiet = TRUE)

# Compute centroids (suppress longlat warning)
suppressWarnings({
  centroids <- st_centroid(mun_sf)
})
centroids$cvegeo_shp <- paste0(centroids$CVE_ENT, centroids$CVE_MUN)

# Hub states
hub_codes <- c("25", "08", "12", "10", "16")  # Sinaloa, Chihuahua, Guerrero, Durango, Michoacán
hub_centroids <- centroids |> filter(CVE_ENT %in% hub_codes)

# For each municipality in the data, compute min distance to any hub municipality centroid
obs_munis <- unique(df$cvegeo)
obs_centroids <- centroids |> filter(cvegeo_shp %in% obs_munis)

cat(sprintf("  Municipalities in data: %d\n", length(obs_munis)))
cat(sprintf("  Matched to shapefile:   %d\n", nrow(obs_centroids)))

# Distance matrix (meters → km)
suppressWarnings({
  dist_mat <- st_distance(obs_centroids, hub_centroids)
})
min_dist_m <- apply(dist_mat, 1, min)
dist_df <- data.frame(
  cvegeo    = obs_centroids$cvegeo_shp,
  dist_hub_km = as.numeric(min_dist_m) / 1000,
  stringsAsFactors = FALSE
)

# Assign bands
dist_band_labels <- c(
  near = "0-200km",
  mid = "200-500km",
  far = "500+km"
)

dist_df <- dist_df |> mutate(
  dist_band = case_when(
    dist_hub_km <= 200  ~ "near",
    dist_hub_km <= 500  ~ "mid",
    TRUE                ~ "far"
  ),
  dist_band = factor(dist_band, levels = names(dist_band_labels))
)

df <- df |> left_join(dist_df, by = "cvegeo")
cat("  Distance band distribution:\n")
print(table(
  factor(df$dist_band, levels = names(dist_band_labels),
         labels = unname(dist_band_labels)),
  useNA = "ifany"
))
cat("\n")

# Drop obs with no match
df_dist <- df |> filter(!is.na(dist_band))

# Interaction model
m_dist <- feols(ln_quantity ~ ln_price:dist_band + quality_good + age +
                  gender_f + educ_f | state_code,
                data = df_dist, cluster = ~state_code)

cat("  Elasticities by distance band:\n")
dist_coefs <- coef(m_dist)
dist_ses   <- se(m_dist)
band_names <- grep("ln_price:", names(dist_coefs), value = TRUE)
for (bn in band_names) {
  band_code <- sub("^ln_price:dist_band", "", bn)
  band_label <- dist_band_labels[[band_code]]
  cat(sprintf("    %s : %.4f  (SE = %.4f)\n", band_label, dist_coefs[bn], dist_ses[bn]))
}
cat(sprintf("  N = %d\n", m_dist$nobs))

# Formal test: are band elasticities equal?
if (all(c("ln_price:dist_bandnear", "ln_price:dist_bandmid",
          "ln_price:dist_bandfar") %in% names(dist_coefs))) {
  wt <- tryCatch(
    wald(
      m_dist,
      c(
        "ln_price:dist_bandnear = ln_price:dist_bandmid",
        "ln_price:dist_bandnear = ln_price:dist_bandfar"
      )
    ),
    error = function(e) NULL
  )
  if (is.list(wt) && all(c("stat", "p") %in% names(wt))) {
    cat(sprintf("  Wald test (equal elasticities across bands): stat = %.3f, p = %.4f\n",
                wt$stat, wt$p))
  } else {
    cat("  Wald test unavailable under the current fixest return format.\n")
  }
}
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# 4. MUNICIPALITY-LEVEL CLUSTERING
# ══════════════════════════════════════════════════════════════════════════════
cat("=== 4. CLUSTERING AT MUNICIPALITY LEVEL (vs STATE) ===\n")

m_muni <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
                | state_code, data = df, cluster = ~cvegeo)

cat(sprintf("  Cluster = state:        ln_price SE = %.4f\n",
            se(baseline)["ln_price"]))
cat(sprintf("  Cluster = municipality: ln_price SE = %.4f\n",
            se(m_muni)["ln_price"]))
cat(sprintf("  Ratio (muni/state):     %.3f\n",
            se(m_muni)["ln_price"] / se(baseline)["ln_price"]))

cat("\n  Full comparison of key SEs:\n")
vars_compare <- c("ln_price", "quality_good", "age")
for (v in vars_compare) {
  cat(sprintf("    %-15s  state SE=%.4f   muni SE=%.4f\n",
              v, se(baseline)[v], se(m_muni)[v]))
}

cat(sprintf("\n  Number of state clusters:        %d\n",
            n_distinct(df$state_code)))
cat(sprintf("  Number of municipality clusters: %d\n",
            n_distinct(df$cvegeo)))

cat("\n=== DONE ===\n")
