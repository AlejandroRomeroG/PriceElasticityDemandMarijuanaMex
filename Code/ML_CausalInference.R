###############################################################################
#  ML_CausalInference.R
#  Machine Learning-Augmented Causal Inference for Price Elasticity
#  -----------------------------------------------------------------------
#  1. Double/Debiased Machine Learning (DML)
#  2. LASSO variable selection + Post-LASSO OLS
#  3. Random Forest variable importance
#  4. Binscatter / nonparametric demand curve
#  5. Sorted Effects (heterogeneity in elasticity by decile)
###############################################################################

# ---- 0. Packages -----------------------------------------------------------
needed <- c("arrow", "dplyr", "stringr", "lubridate", "glmnet", "ranger",
            "sandwich", "lmtest")
for (pkg in needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

library(arrow)
library(dplyr)
library(stringr)
library(lubridate)
library(glmnet)
library(ranger)
library(sandwich)
library(lmtest)

set.seed(42)

# ---- 1. Data prep ----------------------------------------------------------
df <- read_parquet("Data/Marijuana_Prices_in_Mexico_clean.parquet")

df <- df %>%
  mutate(
    ln_quantity   = log(quantity_g),
    ln_price      = log(price_mxn_g),
    quality_good  = as.integer(quality == "good"),
    cvegeo        = str_pad(as.character(muni_code), 5, pad = "0"),
    state_code    = str_sub(cvegeo, 1, 2),
    month_str     = format(floor_date(as.Date(purchase_date), "month"), "%Y-%m"),
    educ_f        = factor(education,
                           levels = c("Elementary school", "Middle school",
                                      "High school", "Bachelor", "Graduate")),
    gender_f      = factor(gender)
  )

cat("\n========== DATA SUMMARY ==========\n")
cat("N =", nrow(df), "\n")
cat("ln_price:    mean =", round(mean(df$ln_price), 3),
    " sd =", round(sd(df$ln_price), 3), "\n")
cat("ln_quantity: mean =", round(mean(df$ln_quantity), 3),
    " sd =", round(sd(df$ln_quantity), 3), "\n")

# Build model matrix of controls (no intercept, no ln_price)
# Controls: quality_good, age, gender_f, educ_f, state dummies
X_formula <- ~ quality_good + age + gender_f + educ_f + factor(state_code) - 1
X_mat <- model.matrix(X_formula, data = df)
cat("Control matrix dimensions:", nrow(X_mat), "x", ncol(X_mat), "\n")

###############################################################################
# 1. DOUBLE/DEBIASED MACHINE LEARNING (DML)
###############################################################################
cat("\n\n================================================================\n")
cat("  1. DOUBLE/DEBIASED MACHINE LEARNING (DML)\n")
cat("================================================================\n")

K <- 5  # folds
n <- nrow(df)
fold_id <- sample(rep(1:K, length.out = n))

# Storage for residuals
resid_Y <- numeric(n)  # ln_quantity residuals
resid_D <- numeric(n)  # ln_price residuals

cat("\nRunning 5-fold cross-fitting with Random Forest...\n")

for (k in 1:K) {
  train_idx <- which(fold_id != k)
  test_idx  <- which(fold_id == k)

  # First stage: predict ln_price from controls
  rf_D <- ranger(x = X_mat[train_idx, ], y = df$ln_price[train_idx],
                 num.trees = 500, min.node.size = 5, num.threads = 1)
  pred_D <- predict(rf_D, data = X_mat[test_idx, ])$predictions
  resid_D[test_idx] <- df$ln_price[test_idx] - pred_D

  # Second stage: predict ln_quantity from controls
  rf_Y <- ranger(x = X_mat[train_idx, ], y = df$ln_quantity[train_idx],
                 num.trees = 500, min.node.size = 5, num.threads = 1)
  pred_Y <- predict(rf_Y, data = X_mat[test_idx, ])$predictions
  resid_Y[test_idx] <- df$ln_quantity[test_idx] - pred_Y

  cat("  Fold", k, "done. Test size:", length(test_idx), "\n")
}

# DML estimate: regress Y-residuals on D-residuals (no intercept by Frisch-Waugh)
dml_fit <- lm(resid_Y ~ resid_D - 1)
dml_coef <- coef(dml_fit)
dml_se   <- sqrt(vcovHC(dml_fit, type = "HC1")[1, 1])

cat("\n--- DML Results (Random Forest, 5-fold) ---\n")
cat("  Elasticity estimate: ", round(dml_coef, 4), "\n")
cat("  Robust SE:           ", round(dml_se, 4), "\n")
cat("  t-statistic:         ", round(dml_coef / dml_se, 2), "\n")
cat("  95% CI:              [", round(dml_coef - 1.96 * dml_se, 4), ",",
    round(dml_coef + 1.96 * dml_se, 4), "]\n")

# Also run DML with LASSO
cat("\nRunning 5-fold cross-fitting with LASSO...\n")
resid_Y_lasso <- numeric(n)
resid_D_lasso <- numeric(n)

for (k in 1:K) {
  train_idx <- which(fold_id != k)
  test_idx  <- which(fold_id == k)

  # LASSO for ln_price
  cv_D <- cv.glmnet(X_mat[train_idx, ], df$ln_price[train_idx],
                     alpha = 1, nfolds = 5)
  pred_D <- predict(cv_D, newx = X_mat[test_idx, ], s = "lambda.min")
  resid_D_lasso[test_idx] <- df$ln_price[test_idx] - as.numeric(pred_D)

  # LASSO for ln_quantity
  cv_Y <- cv.glmnet(X_mat[train_idx, ], df$ln_quantity[train_idx],
                     alpha = 1, nfolds = 5)
  pred_Y <- predict(cv_Y, newx = X_mat[test_idx, ], s = "lambda.min")
  resid_Y_lasso[test_idx] <- df$ln_quantity[test_idx] - as.numeric(pred_Y)

  cat("  Fold", k, "done.\n")
}

dml_lasso_fit <- lm(resid_Y_lasso ~ resid_D_lasso - 1)
dml_lasso_coef <- coef(dml_lasso_fit)
dml_lasso_se   <- sqrt(vcovHC(dml_lasso_fit, type = "HC1")[1, 1])

cat("\n--- DML Results (LASSO, 5-fold) ---\n")
cat("  Elasticity estimate: ", round(dml_lasso_coef, 4), "\n")
cat("  Robust SE:           ", round(dml_lasso_se, 4), "\n")
cat("  t-statistic:         ", round(dml_lasso_coef / dml_lasso_se, 2), "\n")
cat("  95% CI:              [", round(dml_lasso_coef - 1.96 * dml_lasso_se, 4), ",",
    round(dml_lasso_coef + 1.96 * dml_lasso_se, 4), "]\n")


###############################################################################
# 2. LASSO VARIABLE SELECTION WITH INTERACTIONS
###############################################################################
cat("\n\n================================================================\n")
cat("  2. LASSO VARIABLE SELECTION (WITH INTERACTIONS)\n")
cat("================================================================\n")

# Build expanded model matrix with interactions
X_full_formula <- ~ (ln_price + quality_good + age + gender_f + educ_f +
                       factor(state_code))^2 - 1
X_full <- model.matrix(X_full_formula, data = df)
cat("Full interaction matrix dimensions:", nrow(X_full), "x", ncol(X_full), "\n")

# Fit LASSO with cross-validation
cv_lasso <- cv.glmnet(X_full, df$ln_quantity, alpha = 1, nfolds = 10)
cat("Lambda.min:", round(cv_lasso$lambda.min, 6), "\n")
cat("Lambda.1se:", round(cv_lasso$lambda.1se, 6), "\n")

# Coefficients at lambda.min
coefs_min <- coef(cv_lasso, s = "lambda.min")
selected_min <- rownames(coefs_min)[which(coefs_min[, 1] != 0)]
selected_min <- selected_min[selected_min != "(Intercept)"]

cat("\n--- Variables selected by LASSO (lambda.min) ---\n")
cat("  Number of variables selected:", length(selected_min), "\n")
# Show the selected variables with their coefficients
coef_vals <- coefs_min[c("(Intercept)", selected_min), 1]
for (i in seq_along(coef_vals)) {
  cat(sprintf("  %-50s %8.5f\n", names(coef_vals)[i], coef_vals[i]))
}

# Coefficients at lambda.1se (more parsimonious)
coefs_1se <- coef(cv_lasso, s = "lambda.1se")
selected_1se <- rownames(coefs_1se)[which(coefs_1se[, 1] != 0)]
selected_1se <- selected_1se[selected_1se != "(Intercept)"]

cat("\n--- Variables selected by LASSO (lambda.1se) ---\n")
cat("  Number of variables selected:", length(selected_1se), "\n")
coef_vals_1se <- coefs_1se[c("(Intercept)", selected_1se), 1]
for (i in seq_along(coef_vals_1se)) {
  cat(sprintf("  %-50s %8.5f\n", names(coef_vals_1se)[i], coef_vals_1se[i]))
}

# Post-LASSO OLS: refit OLS using only LASSO-selected variables
if (length(selected_min) > 0) {
  X_postlasso <- X_full[, selected_min, drop = FALSE]
  postlasso_fit <- lm(df$ln_quantity ~ X_postlasso)
  postlasso_robust <- coeftest(postlasso_fit, vcov = vcovHC(postlasso_fit, type = "HC1"))

  # Find the ln_price coefficient
  price_idx <- grep("ln_price$", rownames(postlasso_robust))
  if (length(price_idx) == 0) price_idx <- grep("ln_price", rownames(postlasso_robust))[1]

  cat("\n--- Post-LASSO OLS (lambda.min selection) ---\n")
  if (length(price_idx) > 0) {
    cat("  ln_price elasticity:  ", round(postlasso_robust[price_idx, 1], 4), "\n")
    cat("  Robust SE:            ", round(postlasso_robust[price_idx, 2], 4), "\n")
    cat("  t-stat:               ", round(postlasso_robust[price_idx, 3], 2), "\n")
  }
  cat("  R-squared:            ", round(summary(postlasso_fit)$r.squared, 4), "\n")
  cat("  Adj. R-squared:       ", round(summary(postlasso_fit)$adj.r.squared, 4), "\n")

  # Print ALL post-LASSO coefficients
  cat("\n  Full Post-LASSO OLS coefficients:\n")
  for (i in 1:nrow(postlasso_robust)) {
    cat(sprintf("    %-50s %8.4f  (SE=%6.4f, t=%6.2f)\n",
                rownames(postlasso_robust)[i],
                postlasso_robust[i, 1], postlasso_robust[i, 2], postlasso_robust[i, 3]))
  }
}


###############################################################################
# 3. RANDOM FOREST VARIABLE IMPORTANCE
###############################################################################
cat("\n\n================================================================\n")
cat("  3. RANDOM FOREST VARIABLE IMPORTANCE\n")
cat("================================================================\n")

# Build a clean data frame for RF
rf_data <- df %>%
  select(ln_quantity, ln_price, quality_good, age, gender_f, educ_f, state_code) %>%
  mutate(state_code = factor(state_code))

rf_full <- ranger(ln_quantity ~ ., data = rf_data,
                  num.trees = 1000, importance = "impurity",
                  min.node.size = 5, num.threads = 1)

cat("\nRandom Forest fit:\n")
cat("  OOB R-squared:", round(rf_full$r.squared, 4), "\n")
cat("  OOB MSE:      ", round(rf_full$prediction.error, 4), "\n")

# Variable importance
imp <- sort(rf_full$variable.importance, decreasing = TRUE)
cat("\n--- Variable Importance (Impurity-based) ---\n")
for (i in seq_along(imp)) {
  pct <- 100 * imp[i] / sum(imp)
  bar <- paste(rep("|", round(pct / 2)), collapse = "")
  cat(sprintf("  %-15s %8.2f  (%5.1f%%)  %s\n",
              names(imp)[i], imp[i], pct, bar))
}

# Also run with permutation importance
rf_perm <- ranger(ln_quantity ~ ., data = rf_data,
                  num.trees = 1000, importance = "permutation",
                  min.node.size = 5, num.threads = 1)

imp_perm <- sort(rf_perm$variable.importance, decreasing = TRUE)
cat("\n--- Variable Importance (Permutation-based) ---\n")
for (i in seq_along(imp_perm)) {
  pct <- 100 * imp_perm[i] / sum(imp_perm)
  bar <- paste(rep("|", round(pct / 2)), collapse = "")
  cat(sprintf("  %-15s %8.4f  (%5.1f%%)  %s\n",
              names(imp_perm)[i], imp_perm[i], pct, bar))
}


###############################################################################
# 4. BINSCATTER / NONPARAMETRIC DEMAND CURVE
###############################################################################
cat("\n\n================================================================\n")
cat("  4. BINSCATTER / NONPARAMETRIC DEMAND CURVE\n")
cat("================================================================\n")

# Residualize both ln_price and ln_quantity on controls + state FE
ctrl_formula <- ~ quality_good + age + gender_f + educ_f + factor(state_code)
X_ctrl <- model.matrix(ctrl_formula, data = df)

# Residualize ln_quantity
fit_y <- lm(df$ln_quantity ~ X_ctrl)
resid_qty <- residuals(fit_y)

# Residualize ln_price
fit_d <- lm(df$ln_price ~ X_ctrl)
resid_prc <- residuals(fit_d)

# Create 20 equal-sized bins of residualized ln_price
n_bins <- 20
bin_breaks <- quantile(resid_prc, probs = seq(0, 1, length.out = n_bins + 1))
bin_id <- cut(resid_prc, breaks = bin_breaks, include.lowest = TRUE, labels = FALSE)

binscatter <- data.frame(
  bin       = 1:n_bins,
  price_mid = tapply(resid_prc, bin_id, mean),
  qty_mean  = tapply(resid_qty, bin_id, mean),
  n_obs     = as.integer(table(bin_id))
)

cat("\n--- Binscatter: Residualized ln_price vs ln_quantity (20 bins) ---\n")
cat(sprintf("  %-4s  %10s  %10s  %5s\n", "Bin", "Price.mid", "Qty.mean", "N"))
cat(paste(rep("-", 40), collapse = ""), "\n")
for (i in 1:nrow(binscatter)) {
  cat(sprintf("  %-4d  %10.4f  %10.4f  %5d\n",
              binscatter$bin[i], binscatter$price_mid[i],
              binscatter$qty_mean[i], binscatter$n_obs[i]))
}

# Test for linearity: compare linear fit to quadratic
lin_fit  <- lm(qty_mean ~ price_mid, data = binscatter)
quad_fit <- lm(qty_mean ~ price_mid + I(price_mid^2), data = binscatter)

cat("\n--- Linearity test ---\n")
cat("  Linear slope:        ", round(coef(lin_fit)[2], 4), "\n")
cat("  Linear R-squared:    ", round(summary(lin_fit)$r.squared, 4), "\n")
cat("  Quadratic coefs:     ", round(coef(quad_fit)[2], 4), "(linear),",
    round(coef(quad_fit)[3], 4), "(squared)\n")
cat("  Quadratic R-squared: ", round(summary(quad_fit)$r.squared, 4), "\n")

# F-test for quadratic term
anova_test <- anova(lin_fit, quad_fit)
cat("  F-test for curvature: F =", round(anova_test$F[2], 3),
    ", p =", round(anova_test$`Pr(>F)`[2], 4), "\n")
if (anova_test$`Pr(>F)`[2] < 0.05) {
  cat("  => Significant curvature detected (log-log may not be fully linear)\n")
} else {
  cat("  => No significant curvature (log-log specification supported)\n")
}

# Correlation between binned price and quantity
cat("  Correlation (binned):", round(cor(binscatter$price_mid, binscatter$qty_mean), 4), "\n")


###############################################################################
# 5. SORTED EFFECTS (HETEROGENEITY IN ELASTICITY BY DECILE)
###############################################################################
cat("\n\n================================================================\n")
cat("  5. SORTED EFFECTS (ELASTICITY HETEROGENEITY BY DECILE)\n")
cat("================================================================\n")

# Auxiliary model: interact ln_price with all controls to get heterogeneous effects
df_sorted <- df %>%
  mutate(
    lnp_quality = ln_price * quality_good,
    lnp_age     = ln_price * age,
    lnp_male    = ln_price * as.integer(gender_f == "Male"),
    educ_num    = as.integer(educ_f),
    lnp_educ    = ln_price * educ_num
  )

# Fit the interacted model
sorted_fit <- lm(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f +
                    factor(state_code) +
                    lnp_quality + lnp_age + lnp_male + lnp_educ,
                  data = df_sorted)

cat("\nInteracted model summary (key coefficients):\n")
sorted_coefs <- coeftest(sorted_fit, vcov = vcovHC(sorted_fit, type = "HC1"))
# Print interaction terms
interaction_vars <- c("ln_price", "lnp_quality", "lnp_age", "lnp_male", "lnp_educ")
for (v in interaction_vars) {
  idx <- which(rownames(sorted_coefs) == v)
  if (length(idx) > 0) {
    cat(sprintf("  %-15s  %7.4f  (SE=%6.4f, t=%6.2f)\n",
                v, sorted_coefs[idx, 1], sorted_coefs[idx, 2], sorted_coefs[idx, 3]))
  }
}

# Compute individual-level predicted elasticity
# elasticity_i = beta_price + beta_lnp_quality * quality_good_i + ...
beta_p <- coef(sorted_fit)["ln_price"]
beta_pq <- coef(sorted_fit)["lnp_quality"]
beta_pa <- coef(sorted_fit)["lnp_age"]
beta_pm <- coef(sorted_fit)["lnp_male"]
beta_pe <- coef(sorted_fit)["lnp_educ"]

df_sorted$pred_elasticity <- beta_p +
  beta_pq * df_sorted$quality_good +
  beta_pa * df_sorted$age +
  beta_pm * as.integer(df_sorted$gender_f == "Male") +
  beta_pe * df_sorted$educ_num

cat("\n--- Distribution of predicted elasticities ---\n")
cat("  Mean:   ", round(mean(df_sorted$pred_elasticity), 4), "\n")
cat("  Median: ", round(median(df_sorted$pred_elasticity), 4), "\n")
cat("  SD:     ", round(sd(df_sorted$pred_elasticity), 4), "\n")
cat("  Min:    ", round(min(df_sorted$pred_elasticity), 4), "\n")
cat("  Max:    ", round(max(df_sorted$pred_elasticity), 4), "\n")

# Split into deciles and estimate local elasticity in each
df_sorted$decile <- cut(df_sorted$pred_elasticity,
                        breaks = quantile(df_sorted$pred_elasticity,
                                          probs = seq(0, 1, 0.1)),
                        include.lowest = TRUE, labels = 1:10)

cat("\n--- Sorted Effects: Elasticity by Decile ---\n")
cat(sprintf("  %-7s  %10s  %10s  %8s  %5s\n",
            "Decile", "Elasticity", "Rob.SE", "t-stat", "N"))
cat(paste(rep("-", 50), collapse = ""), "\n")

decile_results <- data.frame(decile = 1:10, elasticity = NA, se = NA, n = NA)

for (d in 1:10) {
  sub <- df_sorted[df_sorted$decile == d, ]
  if (nrow(sub) > 30) {
    # Simple OLS with state FE in each decile
    fit_d <- tryCatch({
      lm(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f,
         data = sub)
    }, error = function(e) NULL)

    if (!is.null(fit_d)) {
      robust_d <- tryCatch(
        coeftest(fit_d, vcov = vcovHC(fit_d, type = "HC1")),
        error = function(e) NULL
      )
      if (!is.null(robust_d)) {
        pidx <- which(rownames(robust_d) == "ln_price")
        decile_results$elasticity[d] <- robust_d[pidx, 1]
        decile_results$se[d]         <- robust_d[pidx, 2]
        decile_results$n[d]          <- nrow(sub)
        cat(sprintf("  %-7d  %10.4f  %10.4f  %8.2f  %5d\n",
                    d, robust_d[pidx, 1], robust_d[pidx, 2],
                    robust_d[pidx, 3], nrow(sub)))
      }
    }
  }
}

cat("\n--- Summary of Heterogeneity ---\n")
valid <- !is.na(decile_results$elasticity)
if (sum(valid) > 0) {
  cat("  Range of local elasticities: [",
      round(min(decile_results$elasticity[valid]), 4), ",",
      round(max(decile_results$elasticity[valid]), 4), "]\n")
  cat("  Spread (max - min):          ",
      round(max(decile_results$elasticity[valid]) - min(decile_results$elasticity[valid]), 4), "\n")
  cat("  Mean local elasticity:       ",
      round(mean(decile_results$elasticity[valid]), 4), "\n")
}


###############################################################################
# SUMMARY TABLE
###############################################################################
cat("\n\n================================================================\n")
cat("  SUMMARY: ALL ELASTICITY ESTIMATES\n")
cat("================================================================\n")
cat(sprintf("  %-40s  %8s  %8s\n", "Method", "Estimate", "SE"))
cat(paste(rep("-", 62), collapse = ""), "\n")
cat(sprintf("  %-40s  %8.4f  %8.4f\n", "OLS benchmark (St x Mo FE)", -0.652, 0.022))
cat(sprintf("  %-40s  %8.4f  %8.4f\n", "DML (Random Forest)", dml_coef, dml_se))
cat(sprintf("  %-40s  %8.4f  %8.4f\n", "DML (LASSO)", dml_lasso_coef, dml_lasso_se))
if (exists("postlasso_robust") && length(price_idx) > 0) {
  cat(sprintf("  %-40s  %8.4f  %8.4f\n", "Post-LASSO OLS",
              postlasso_robust[price_idx, 1], postlasso_robust[price_idx, 2]))
}
cat(sprintf("  %-40s  %8.4f  %8s\n", "Binscatter slope (nonparametric)",
            coef(lin_fit)[2], "—"))
if (sum(valid) > 0) {
  cat(sprintf("  %-40s  [%6.3f, %6.3f]\n", "Sorted effects range",
              min(decile_results$elasticity[valid]),
              max(decile_results$elasticity[valid])))
}

cat("\nDone.\n")
