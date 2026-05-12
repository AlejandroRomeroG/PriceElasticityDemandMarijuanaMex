## BootstrapInference.R
## Bootstrap, permutation, subsampling, and influence diagnostics
## for marijuana demand elasticity estimation

project_lib <- file.path(getwd(), ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(here); library(arrow); library(dplyr); library(stringr); library(lubridate); library(fixest)

fixest::setFixest_estimation(fixef.rm = "infinite_coef")

if (!requireNamespace("fwildclusterboot", quietly = TRUE)) {
  stop(
    "Package 'fwildclusterboot' is required for the wild cluster bootstrap. ",
    "Install it with install.packages('fwildclusterboot') before running this script.",
    call. = FALSE
  )
}

if (requireNamespace("dqrng", quietly = TRUE)) {
  dqrng::dqset.seed(42)
}

set.seed(42)

df <- read_parquet(here("Data","Marijuana_Prices_in_Mexico_clean.parquet"))
df <- df |> mutate(
  ln_quantity = log(quantity_g), ln_price = log(price_mxn_g),
  quality_good = as.integer(quality=="good"),
  cvegeo = str_pad(as.character(muni_code),5,pad="0"),
  state_code = str_sub(cvegeo,1,2),
  state_fe = as.integer(factor(state_code)),
  purchase_date = as.Date(purchase_date),
  month_str = format(floor_date(purchase_date,"month"),"%Y-%m"),
  state_month = paste0(state_code, "_", month_str),
  state_month_f = factor(state_month),
  month_f = factor(month_str),
  educ_f = factor(education, levels=c("Elementary school","Middle school","High school","Bachelor","Graduate")),
  gender_f = factor(gender)
)

cat("======================================================\n")
cat("  BOOTSTRAP INFERENCE FOR MARIJUANA DEMAND ELASTICITY\n")
cat("======================================================\n\n")
cat(sprintf("N = %d observations, %d unique states\n\n", nrow(df), length(unique(df$state_fe))))

# --- Benchmark estimation: state-by-month FE ---
baseline <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f + state_month_f,
                  data = df, cluster = ~state_fe)
obs_beta <- coef(baseline)["ln_price"]
obs_se   <- se(baseline)["ln_price"]
cat(sprintf("Baseline elasticity: %.4f  (SE = %.4f)\n\n", obs_beta, obs_se))

# ============================================================
# 1. WILD CLUSTER BOOTSTRAP (Webb weights, 9999 reps)
# ============================================================
cat("------------------------------------------------------\n")
cat("1. WILD CLUSTER BOOTSTRAP (Webb weights, 9999 replications)\n")
cat("------------------------------------------------------\n")

wild_boot <- fwildclusterboot::boottest(
  baseline,
  clustid = "state_fe",
  param = "ln_price",
  B = 9999,
  type = "webb",
  engine = "R"
)

cat("  Wild cluster bootstrap summary:\n")
print(summary(wild_boot))
cat("\n  Wild cluster bootstrap 95% CI:\n")
wild_ci <- as.numeric(stats::confint(wild_boot))
print(stats::confint(wild_boot))
cat("\n")

# ============================================================
# 2. PAIRS CLUSTER BOOTSTRAP (secondary sensitivity, 9999 reps)
# ============================================================
cat("------------------------------------------------------\n")
cat("2. PAIRS CLUSTER BOOTSTRAP (secondary sensitivity, 9999 replications)\n")
cat("------------------------------------------------------\n")

states <- unique(df$state_fe)
n_states <- length(states)
B <- 9999
boot_betas <- numeric(B)

for (b in seq_len(B)) {
  # resample states with replacement
  sampled_states <- sample(states, n_states, replace = TRUE)
  # build bootstrap dataset (handle duplicated states)
  boot_list <- lapply(seq_along(sampled_states), function(i) {
    sub <- df[df$state_fe == sampled_states[i], ]
    sub$boot_state <- paste0(sampled_states[i], "_", i)
    sub$boot_state_month <- paste0(sub$boot_state, "_", sub$month_str)
    sub
  })
  boot_df <- do.call(rbind, boot_list)

  fit <- tryCatch(
    feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f + factor(boot_state_month),
          data = boot_df, notes = FALSE),
    error = function(e) NULL
  )
  boot_betas[b] <- if (!is.null(fit)) coef(fit)["ln_price"] else NA_real_
}

boot_betas_clean <- boot_betas[!is.na(boot_betas)]
cat(sprintf("  Successful replications: %d / %d\n", length(boot_betas_clean), B))
cat(sprintf("  Bootstrap mean:    %.4f\n", mean(boot_betas_clean)))
cat(sprintf("  Bootstrap median:  %.4f\n", median(boot_betas_clean)))
cat(sprintf("  Bootstrap SD:      %.4f\n", sd(boot_betas_clean)))
q_boot <- quantile(boot_betas_clean, c(0.025, 0.975))
cat(sprintf("  95%% CI (percentile): [%.4f, %.4f]\n\n", q_boot[1], q_boot[2]))

output_dir <- here("LaTeX", "tables")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

fmt <- function(x, digits = 3) {
  if (is.na(x) || !is.finite(x)) return("--")
  formatC(x, format = "f", digits = digits)
}

fmt_ci <- function(lower, upper, digits = 3) {
  paste0("[", fmt(lower, digits), ", ", fmt(upper, digits), "]")
}

conv_ci <- obs_beta + c(-1.96, 1.96) * obs_se

writeLines(
  c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Bootstrap inference for the benchmark elasticity}",
    "\\label{tab:bootstrap-inference}",
    "\\centering",
    "\\begin{tabular}[t]{lccc}",
    "\\toprule",
    "Method & Estimate & SE & 95\\% CI\\\\",
    "\\midrule",
    paste0("State-clustered SE & ", fmt(obs_beta), " & ", fmt(obs_se), " & ",
           fmt_ci(conv_ci[1], conv_ci[2]), "\\\\"),
    paste0("Wild cluster bootstrap & ", fmt(obs_beta), " & -- & ",
           fmt_ci(wild_ci[1], wild_ci[2]), "\\\\"),
    paste0("Pairs cluster bootstrap & ", fmt(mean(boot_betas_clean)), " & ",
           fmt(sd(boot_betas_clean)), " & ", fmt_ci(q_boot[1], q_boot[2]), "\\\\"),
    "\\bottomrule",
    "\\multicolumn{4}{l}{\\rule{0pt}{1em}Benchmark model: state-by-month FE, quality, age, gender, and education controls.}\\\\",
    "\\multicolumn{4}{l}{\\rule{0pt}{1em}Wild bootstrap uses Webb weights, 9,999 replications, and state clusters.}\\\\",
    "\\multicolumn{4}{l}{\\rule{0pt}{1em}Pairs cluster bootstrap resamples states with replacement across 9,999 replications.}\\\\",
    "\\end{tabular}",
    "\\end{table}"
  ),
  file.path(output_dir, "table_bootstrap_inference.tex")
)

# ============================================================
# 3. PERMUTATION TEST (permute ln_price within states, 1000 reps)
# ============================================================
cat("------------------------------------------------------\n")
cat("3. PERMUTATION TEST (1000 permutations)\n")
cat("------------------------------------------------------\n")

P <- 1000
perm_betas <- numeric(P)

for (p in seq_len(P)) {
  df_perm <- df
  # permute ln_price within each state-month cell, the benchmark identifying margin
  df_perm <- df_perm |>
    group_by(state_month) |>
    mutate(ln_price = sample(ln_price, size = dplyr::n())) |>
    ungroup()

  fit_perm <- tryCatch(
    feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f + state_month_f,
          data = df_perm, notes = FALSE),
    error = function(e) NULL
  )
  perm_betas[p] <- if (!is.null(fit_perm)) coef(fit_perm)["ln_price"] else NA_real_
}

perm_betas_clean <- perm_betas[!is.na(perm_betas)]
perm_pval <- mean(abs(perm_betas_clean) >= abs(obs_beta))
cat(sprintf("  Observed |elasticity|: %.4f\n", abs(obs_beta)))
cat(sprintf("  Permutation mean:      %.4f\n", mean(perm_betas_clean)))
cat(sprintf("  Permutation SD:        %.4f\n", sd(perm_betas_clean)))
cat(sprintf("  Permutation p-value:   %.4f  (fraction |perm beta| >= |obs beta|)\n\n", perm_pval))

# ============================================================
# 4. SUBSAMPLING STABILITY (50%, 70%, 90% x 500 reps)
# ============================================================
cat("------------------------------------------------------\n")
cat("4. SUBSAMPLING STABILITY (500 reps per fraction)\n")
cat("------------------------------------------------------\n")

fracs <- c(0.50, 0.70, 0.90)
S <- 500

for (frac in fracs) {
  sub_betas <- numeric(S)
  n_sub <- round(nrow(df) * frac)

  for (s in seq_len(S)) {
    idx <- sample(seq_len(nrow(df)), n_sub, replace = FALSE)
    df_sub <- df[idx, ]

    fit_sub <- tryCatch(
      feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f + factor(state_month),
            data = df_sub, notes = FALSE),
      error = function(e) NULL
    )
    sub_betas[s] <- if (!is.null(fit_sub)) coef(fit_sub)["ln_price"] else NA_real_
  }

  sub_clean <- sub_betas[!is.na(sub_betas)]
  cat(sprintf("  Fraction = %.0f%%: mean = %.4f, SD = %.4f  (%d/%d successful)\n",
              frac * 100, mean(sub_clean), sd(sub_clean), length(sub_clean), S))
}
cat("\n")

# ============================================================
# 5. INFLUENCE DIAGNOSTICS (lm with state-month dummies)
# ============================================================
cat("------------------------------------------------------\n")
cat("5. INFLUENCE DIAGNOSTICS (Cook's distance, top 10)\n")
cat("------------------------------------------------------\n")

fit_lm <- lm(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f +
               factor(state_month),
             data = df)

h_vals <- hatvalues(fit_lm)
cooks  <- cooks.distance(fit_lm)

df$hat_value <- h_vals
df$cooks_d   <- cooks

top10 <- df |>
  arrange(desc(cooks_d)) |>
  slice_head(n = 10) |>
  select(price_mxn_g, quantity_g, state, quality, hat_value, cooks_d)

cat("\n  Top 10 most influential observations by Cook's distance:\n\n")
cat(sprintf("  %-6s  %-10s  %-10s  %-20s  %-8s  %-10s  %-12s\n",
            "Rank", "Price(MXN)", "Qty(g)", "State", "Quality", "Hat value", "Cook's D"))
cat("  ", paste(rep("-", 85), collapse=""), "\n", sep="")

for (i in 1:nrow(top10)) {
  cat(sprintf("  %-6d  %-10.2f  %-10.2f  %-20s  %-8s  %-10.4f  %-12.6f\n",
              i,
              top10$price_mxn_g[i],
              top10$quantity_g[i],
              top10$state[i],
              top10$quality[i],
              top10$hat_value[i],
              top10$cooks_d[i]))
}

cat("\n  Summary of Cook's distance distribution:\n")
cat(sprintf("    Mean:   %.6f\n", mean(cooks, na.rm = TRUE)))
cat(sprintf("    Median: %.6f\n", median(cooks, na.rm = TRUE)))
cat(sprintf("    Max:    %.6f\n", max(cooks, na.rm = TRUE)))
cat(sprintf("    Obs with Cook's D > 4/N (%.6f): %d\n",
            4/nrow(df), sum(cooks > 4/nrow(df), na.rm = TRUE)))

cat("\n======================================================\n")
cat("  DONE\n")
cat("======================================================\n")
