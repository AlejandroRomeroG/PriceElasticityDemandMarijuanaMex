## ============================================================
## Functional Forms & Bounding Exercises
## Price Elasticity of Demand for Marijuana in Mexico
## ============================================================

project_lib <- file.path(getwd(), ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(arrow)
library(dplyr)
library(stringr)
library(lubridate)
library(fixest)
library(AER)       # for tobit
library(MASS)      # for misc

fixest::setFixest_estimation(fixef.rm = "infinite_coef")

# ── Load & prep ──────────────────────────────────────────────
df <- read_parquet("Data/Marijuana_Prices_in_Mexico_clean.parquet") %>%
  mutate(
    ln_quantity  = log(quantity_g),
    ln_price     = log(price_mxn_g),
    quality_good = as.integer(quality == "good"),
    cvegeo       = str_pad(as.character(muni_code), 5, pad = "0"),
    state_code   = str_sub(cvegeo, 1, 2),
    month_str    = format(floor_date(as.Date(purchase_date), "month"), "%Y-%m"),
    state_month   = paste0(state_code, "_", month_str),
    educ_f       = factor(education,
                          levels = c("Elementary school","Middle school",
                                     "High school","Bachelor","Graduate")),
    gender_f     = factor(gender)
  )

cat("N =", nrow(df), "\n")
cat("Mean price (MXN/g):", round(mean(df$price_mxn_g), 3), "\n")
cat("Mean quantity (g):",   round(mean(df$quantity_g), 3), "\n")
cat("Mean ln_price:",       round(mean(df$ln_price), 4), "\n")
cat("Mean ln_quantity:",    round(mean(df$ln_quantity), 4), "\n\n")

mean_p <- mean(df$price_mxn_g)
mean_q <- mean(df$quantity_g)
mean_lnp <- mean(df$ln_price)
mean_lnq <- mean(df$ln_quantity)

# ── 0. Benchmark log-log comparison ──────────────────────────
cat("================================================================\n")
cat("0. BENCHMARK: Log-Log with State-by-Month FE\n")
cat("================================================================\n")
m0 <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
            | state_month, data = df, vcov = "hetero")
summary(m0)
cat("\n")

# ── 1. LINEAR DEMAND MODEL ──────────────────────────────────
cat("================================================================\n")
cat("1. LINEAR DEMAND MODEL: quantity_g ~ price_mxn_g + controls | state-by-month FE\n")
cat("================================================================\n")
m1 <- feols(quantity_g ~ price_mxn_g + quality_good + age + gender_f + educ_f
            | state_month, data = df, vcov = "hetero")
summary(m1)

b_linear <- coef(m1)["price_mxn_g"]
elas_linear <- b_linear * mean_p / mean_q
cat("\nImplied elasticity at the mean (b * p_bar / q_bar):\n")
cat(sprintf("  b = %.4f, mean_p = %.3f, mean_q = %.3f\n", b_linear, mean_p, mean_q))
cat(sprintf("  Elasticity = %.4f\n\n", elas_linear))

# ── 2. SEMI-LOG MODELS ──────────────────────────────────────
cat("================================================================\n")
cat("2a. LIN-LOG: quantity_g ~ ln_price + controls | state-by-month FE\n")
cat("================================================================\n")
m2a <- feols(quantity_g ~ ln_price + quality_good + age + gender_f + educ_f
             | state_month, data = df, vcov = "hetero")
summary(m2a)

b_linlog <- coef(m2a)["ln_price"]
elas_linlog <- b_linlog / mean_q
cat("\nImplied elasticity at the mean (b / q_bar):\n")
cat(sprintf("  b = %.4f, mean_q = %.3f\n", b_linlog, mean_q))
cat(sprintf("  Elasticity = %.4f\n\n", elas_linlog))

cat("================================================================\n")
cat("2b. LOG-LIN: ln_quantity ~ price_mxn_g + controls | state-by-month FE\n")
cat("================================================================\n")
m2b <- feols(ln_quantity ~ price_mxn_g + quality_good + age + gender_f + educ_f
             | state_month, data = df, vcov = "hetero")
summary(m2b)

b_loglin <- coef(m2b)["price_mxn_g"]
elas_loglin <- b_loglin * mean_p
cat("\nImplied elasticity at the mean (b * p_bar):\n")
cat(sprintf("  b = %.6f, mean_p = %.3f\n", b_loglin, mean_p))
cat(sprintf("  Elasticity = %.4f\n\n", elas_loglin))

# ── 3. BOX-COX STYLE GRID ───────────────────────────────────
cat("================================================================\n")
cat("3. BOX-COX GRID SEARCH (dependent variable transformation)\n")
cat("================================================================\n")

lambdas <- c(0, 0.25, 0.5, 0.75, 1)
bc_results <- data.frame(lambda = lambdas, coef_lnprice = NA, se = NA,
                         loglik = NA, AIC = NA, BIC = NA)

n <- nrow(df)

for (i in seq_along(lambdas)) {
  lam <- lambdas[i]
  if (lam == 0) {
    df$y_bc <- log(df$quantity_g)
  } else {
    df$y_bc <- (df$quantity_g^lam - 1) / lam
  }

  # Use feols for the regression
  m_bc <- feols(y_bc ~ ln_price + quality_good + age + gender_f + educ_f
                | state_month, data = df, vcov = "hetero")

  bc_results$coef_lnprice[i] <- coef(m_bc)["ln_price"]
  bc_results$se[i] <- sqrt(vcov(m_bc)["ln_price", "ln_price"])

  # Compute log-likelihood with Jacobian correction for Box-Cox
  resid_bc <- residuals(m_bc)
  sigma2 <- sum(resid_bc^2) / n
  # Log-lik of normal errors + Jacobian of transformation
  ll <- -n/2 * log(2 * pi) - n/2 * log(sigma2) - n/2
  # Jacobian: d(y_bc)/d(q) = q^(lambda-1), so sum log|Jacobian| = (lambda-1)*sum(log(q))
  ll <- ll + (lam - 1) * sum(log(df$quantity_g))

  k <- length(coef(m_bc)) + fixest::fixef(m_bc) %>% lengths() %>% sum() + 1  # +1 for sigma
  bc_results$loglik[i] <- ll
  bc_results$AIC[i] <- -2 * ll + 2 * k
  bc_results$BIC[i] <- -2 * ll + log(n) * k
}

cat("\nBox-Cox Grid Results:\n")
cat(sprintf("%-8s %12s %10s %12s %12s %12s\n",
            "Lambda", "Coef(ln_p)", "SE", "LogLik", "AIC", "BIC"))
cat(paste(rep("-", 70), collapse = ""), "\n")
for (i in 1:nrow(bc_results)) {
  cat(sprintf("%-8.2f %12.4f %10.4f %12.1f %12.1f %12.1f\n",
              bc_results$lambda[i], bc_results$coef_lnprice[i],
              bc_results$se[i], bc_results$loglik[i],
              bc_results$AIC[i], bc_results$BIC[i]))
}
best_lam <- bc_results$lambda[which.min(bc_results$BIC)]
cat(sprintf("\nBest lambda by BIC: %.2f\n", best_lam))
cat(sprintf("Best lambda by AIC: %.2f\n\n", bc_results$lambda[which.min(bc_results$AIC)]))

# ── 4. TOBIT / CENSORING CHECK ──────────────────────────────
cat("================================================================\n")
cat("4. TOBIT MODEL (left-censored at ln(1)=0)\n")
cat("================================================================\n")

# Tobit with AER::tobit (uses survreg underneath)
# Include state-by-month dummies explicitly
m_tobit <- tobit(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f +
                   factor(state_month),
                 left = 0, data = df)

cat("Tobit coefficient on ln_price:\n")
tobit_coef <- coef(m_tobit)["ln_price"]
tobit_se   <- sqrt(vcov(m_tobit)["ln_price", "ln_price"])
cat(sprintf("  Coefficient: %.4f (SE = %.4f)\n", tobit_coef, tobit_se))

# How many observations are at the censoring point?
n_censored <- sum(df$ln_quantity <= 0 + 1e-8)
cat(sprintf("  Observations at or below ln(1)=0: %d (%.1f%%)\n",
            n_censored, 100 * n_censored / nrow(df)))

# OLS comparison
m_ols_comp <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
                    | state_month, data = df, vcov = "hetero")
ols_coef <- coef(m_ols_comp)["ln_price"]
ols_se   <- sqrt(vcov(m_ols_comp)["ln_price", "ln_price"])

cat(sprintf("\n  OLS coefficient:   %.4f (SE = %.4f)\n", ols_coef, ols_se))
cat(sprintf("  Tobit coefficient: %.4f (SE = %.4f)\n", tobit_coef, tobit_se))
cat(sprintf("  Difference:        %.4f\n", tobit_coef - ols_coef))
if (abs(tobit_coef - ols_coef) < 0.02) {
  cat("  --> Tobit and OLS are very similar; censoring at 1g is not a concern.\n\n")
} else {
  cat("  --> Notable difference between Tobit and OLS; censoring may matter.\n\n")
}

# ── 5. OSTER (2019) BOUNDS FOR OMITTED VARIABLE BIAS ────────
cat("================================================================\n")
cat("5. OSTER (2019) BOUNDS FOR OMITTED VARIABLE BIAS\n")
cat("================================================================\n")

# Short regression: only ln_price (no controls, no FE)
m_short <- feols(ln_quantity ~ ln_price, data = df)
beta_short <- coef(m_short)["ln_price"]
R2_short   <- r2(m_short, type = "r2")

# Long regression: full controls + state-by-month FE
m_long <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
                | state_month, data = df)
beta_long <- coef(m_long)["ln_price"]
R2_long   <- r2(m_long, type = "r2")

cat(sprintf("Short regression (price only): beta = %.4f, R² = %.4f\n", beta_short, R2_short))
cat(sprintf("Long regression (full):        beta = %.4f, R² = %.4f\n", beta_long, R2_long))

# Oster's delta* formula (simplified version):
# delta* = (beta_long * (R_long - R_short)) / ((beta_short - beta_long) * (R_max - R_long))
# Under Oster's recommended R_max = min(1, 1.3 * R_long)

R_max_oster <- min(1, 1.3 * R2_long)
cat(sprintf("R_max (Oster recommendation 1.3*R²_long): %.4f\n", R_max_oster))

# delta* for beta = 0
delta_star <- (beta_long * (R2_long - R2_short)) /
  ((beta_short - beta_long) * (R_max_oster - R2_long))

cat(sprintf("\ndelta* (ratio of selection on unobservables to observables\n"))
cat(sprintf("        needed to drive elasticity to zero): %.4f\n", delta_star))

if (abs(delta_star) > 1) {
  cat("  --> |delta*| >> 1: The result is ROBUST to omitted variable bias.\n")
  cat("      Unobservables would need to be", round(abs(delta_star), 1),
      "times as important\n      as observables to eliminate the price effect.\n\n")
} else {
  cat("  --> |delta*| <= 1: Potential vulnerability to OVB.\n\n")
}

# Also compute the bias-adjusted beta under delta=1
# beta_adj = beta_long - delta * (beta_short - beta_long) * (R_max - R_long) / (R_long - R_short)
beta_adj_d1 <- beta_long - 1 * (beta_short - beta_long) * (R_max_oster - R2_long) / (R2_long - R2_short)
cat(sprintf("Bias-adjusted beta (delta=1, R_max=%.3f): %.4f\n", R_max_oster, beta_adj_d1))
cat("  (This is the elasticity if unobservables are equally important as observables.)\n\n")

output_dir <- file.path("LaTeX", "tables")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

fmt <- function(x, digits = 3) {
  if (is.na(x) || !is.finite(x)) return("--")
  formatC(x, format = "f", digits = digits)
}

writeLines(
  c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Oster coefficient-stability diagnostics}",
    "\\label{tab:oster-bounds}",
    "\\centering",
    "\\begin{tabular}[t]{p{0.62\\linewidth}r}",
    "\\toprule",
    "Statistic & Estimate\\\\",
    "\\midrule",
    paste0("Short-specification elasticity & ", fmt(beta_short), "\\\\"),
    paste0("Short-specification $R^2$ & ", fmt(R2_short), "\\\\"),
    paste0("Long-specification elasticity & ", fmt(beta_long), "\\\\"),
    paste0("Long-specification $R^2$ & ", fmt(R2_long), "\\\\"),
    paste0("$R_{\\max}=1.3R^2_{long}$ & ", fmt(R_max_oster), "\\\\"),
    paste0("Bias-adjusted elasticity ($\\delta=1$) & ", fmt(beta_adj_d1), "\\\\"),
    paste0("$|\\delta^*|$ to set elasticity to zero & ", fmt(abs(delta_star), 2), "\\\\"),
    "\\bottomrule",
    "\\multicolumn{2}{p{0.78\\linewidth}}{\\footnotesize \\rule{0pt}{1em}Short specification includes only ln(price); long specification adds controls plus state-by-month FE.}\\\\",
    "\\multicolumn{2}{p{0.78\\linewidth}}{\\footnotesize \\rule{0pt}{1em}The table reports the absolute value of $\\delta^*$; the raw calculation is negative because coefficients are negative.}\\\\",
    "\\end{tabular}",
    "\\end{table}"
  ),
  file.path(output_dir, "table_oster_bounds.tex")
)

# ── 6. MEASUREMENT ERROR BOUNDS (KLEPPER-LEAMER) ────────────
cat("================================================================\n")
cat("6. MEASUREMENT ERROR BOUNDS (Klepper-Leamer / Reverse Regression)\n")
cat("================================================================\n")

# Forward regression: ln_q = a + b_f * ln_p + controls | state-by-month FE
m_fwd <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
               | state_month, data = df)
b_forward <- coef(m_fwd)["ln_price"]
se_forward <- sqrt(vcov(m_fwd, vcov = "hetero")["ln_price", "ln_price"])

# Reverse regression: ln_p = a + b_r * ln_q + controls | state-by-month FE
m_rev <- feols(ln_price ~ ln_quantity + quality_good + age + gender_f + educ_f
               | state_month, data = df)
b_reverse <- coef(m_rev)["ln_quantity"]
se_reverse <- sqrt(vcov(m_rev, vcov = "hetero")["ln_quantity", "ln_quantity"])

# Reciprocal of reverse = upper bound (in absolute terms, lower bound since negative)
b_reverse_recip <- 1 / b_reverse

cat(sprintf("Forward OLS:       beta_f = %.4f (SE = %.4f)\n", b_forward, se_forward))
cat(sprintf("Reverse OLS:       beta_r = %.4f (SE = %.4f)\n", b_reverse, se_reverse))
cat(sprintf("1 / beta_reverse:         = %.4f\n\n", b_reverse_recip))

cat("Klepper-Leamer bounds on the true elasticity:\n")
lower_bound <- min(b_forward, b_reverse_recip)
upper_bound <- max(b_forward, b_reverse_recip)
cat(sprintf("  [ %.4f , %.4f ]\n", lower_bound, upper_bound))
cat("  (Forward OLS is attenuated if price has classical measurement error;\n")
cat("   the reciprocal of the reverse regression gives the other bound.)\n\n")

# ── SUMMARY TABLE ────────────────────────────────────────────
cat("================================================================\n")
cat("SUMMARY: Elasticity Estimates Across Specifications\n")
cat("================================================================\n")
cat(sprintf("%-35s %10s %10s\n", "Specification", "Elasticity", "SE"))
cat(paste(rep("-", 57), collapse = ""), "\n")
cat(sprintf("%-35s %10.4f %10.4f\n", "Log-log (benchmark St×Mo FE)", ols_coef, ols_se))
cat(sprintf("%-35s %10.4f %10s\n",   "Linear (at mean)", elas_linear, "—"))
cat(sprintf("%-35s %10.4f %10s\n",   "Lin-log (at mean)", elas_linlog, "—"))
cat(sprintf("%-35s %10.4f %10s\n",   "Log-lin (at mean)", elas_loglin, "—"))
cat(sprintf("%-35s %10.4f %10.4f\n", "Tobit (left-cens at 0)", tobit_coef, tobit_se))
cat(sprintf("%-35s %10.4f %10s\n",   "Oster bias-adj (delta=1)", beta_adj_d1, "—"))
cat(sprintf("%-35s %10s %10s\n",     "KL bounds",
            sprintf("[%.3f, %.3f]", lower_bound, upper_bound), "—"))
cat(sprintf("%-35s %10.2f %10s\n",   "Oster delta*", delta_star, "—"))
cat(sprintf("%-35s %10.2f %10s\n",   "Best Box-Cox lambda (BIC)", best_lam, "—"))
cat("\nDone.\n")
