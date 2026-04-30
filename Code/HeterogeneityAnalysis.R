################################################################################
# Deep Heterogeneity Analysis — Price Elasticity of Demand for Marijuana (Mexico)
################################################################################

project_lib <- file.path(getwd(), ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(arrow)
library(dplyr)
library(stringr)
library(lubridate)
library(fixest)
# marginaleffects removed — compute AME manually

fixest::setFixest_estimation(fixef.rm = "infinite_coef")

# ── Load and prep ─────────────────────────────────────────────────────────────
df <- read_parquet("Data/Marijuana_Prices_in_Mexico_clean.parquet") |>
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
    gender_f      = factor(gender),
    age_group     = cut(age, breaks = c(17, 24, 34, 44, Inf),
                        labels = c("18-24", "25-34", "35-44", "45+")),
    region        = case_when(
      state_code %in% c("02","03","26","08","05","19","28") ~ "Border",
      state_code %in% c("09","15","21","22","13","29","17") ~ "Central",
      TRUE ~ "Other"
    ),
    region_f      = factor(region, levels = c("Other", "Border", "Central"))
  )

cat("=== Sample overview ===\n")
cat("N =", nrow(df), "\n")
cat("Quality: good =", sum(df$quality_good), " bad =", sum(1 - df$quality_good), "\n")
cat("Gender:", table(df$gender_f), "\n")
cat("Age groups:", table(df$age_group), "\n\n")

################################################################################
# 1. TRIPLE INTERACTION: ln_price × quality_good × gender_f  (state FE)
################################################################################
cat("================================================================\n")
cat("1. TRIPLE INTERACTION MODEL: ln_price * quality_good * gender_f\n")
cat("================================================================\n\n")

m1 <- feols(ln_quantity ~ ln_price * quality_good * gender_f + age + educ_f
            | state_code, data = df, vcov = "HC1")
summary(m1)

# Implied elasticities for each quality×gender cell
b <- coef(m1)
cat("\n--- Implied price elasticities by quality × gender ---\n")

# Reference: bad quality, male (gender_f reference level)
ref_gender <- levels(df$gender_f)[1]  # alphabetical first
cat("Reference gender level:", ref_gender, "\n\n")

# Get coefficient names
cn <- names(b)
cat("Coefficient names:\n"); print(cn); cat("\n")

# Build elasticities manually
e_bad_ref   <- b["ln_price"]
e_good_ref  <- b["ln_price"] + b["ln_price:quality_good"]
# Interaction with gender
gname <- grep("^ln_price:gender_f", cn, value = TRUE)
gname_qual  <- grep("^ln_price:quality_good:gender_f", cn, value = TRUE)

if (length(gname) > 0) {
  e_bad_other  <- b["ln_price"] + b[gname[1]]
  e_good_other <- b["ln_price"] + b["ln_price:quality_good"] +
                  b[gname[1]] + b[gname_qual[1]]
} else {
  # gender might come first
  gname <- grep("gender_f.*:ln_price", cn, value = TRUE)
  gname_qual <- grep("quality_good:gender_f.*:ln_price|ln_price:quality_good:gender_f", cn, value = TRUE)
  e_bad_other  <- b["ln_price"] + b[gname[1]]
  e_good_other <- b["ln_price"] + b["ln_price:quality_good"] +
                  b[gname[1]] + b[gname_qual[1]]
}

other_gender <- setdiff(levels(df$gender_f), ref_gender)
cat(sprintf("  Bad quality,  %s:  %.4f\n", ref_gender, e_bad_ref))
cat(sprintf("  Good quality, %s:  %.4f\n", ref_gender, e_good_ref))
cat(sprintf("  Bad quality,  %s: %.4f\n", other_gender, e_bad_other))
cat(sprintf("  Good quality, %s: %.4f\n", other_gender, e_good_other))
cat("\n")

################################################################################
# 2. AGE-GROUP SEGMENTED ELASTICITIES
################################################################################
cat("================================================================\n")
cat("2. AGE-GROUP SEGMENTED ELASTICITIES\n")
cat("================================================================\n\n")

age_results <- list()
for (ag in levels(df$age_group)) {
  sub <- df |> filter(age_group == ag)
  # Need enough state variation
  n_states <- n_distinct(sub$state_code)
  if (nrow(sub) > 30 & n_states > 2) {
    m <- feols(ln_quantity ~ ln_price + quality_good | state_code,
               data = sub, vcov = "HC1")
    age_results[[ag]] <- list(
      n     = nrow(sub),
      beta  = coef(m)["ln_price"],
      se    = sqrt(vcov(m)["ln_price", "ln_price"]),
      r2    = r2(m, type = "ar2")
    )
    cat(sprintf("Age %s: elasticity = %.4f (SE=%.4f), N=%d, adj.R2=%.3f\n",
                ag, age_results[[ag]]$beta, age_results[[ag]]$se,
                age_results[[ag]]$n, age_results[[ag]]$r2))
  } else {
    cat(sprintf("Age %s: too few obs (N=%d, states=%d)\n", ag, nrow(sub), n_states))
  }
}

# Formal test: pool with interaction
cat("\nFormal test (interaction model):\n")
m2_test <- feols(ln_quantity ~ ln_price * age_group + quality_good | state_code,
                 data = df, vcov = "HC1")
summary(m2_test)
cat("\n")

################################################################################
# 3. EDUCATION GRADIENT
################################################################################
cat("================================================================\n")
cat("3. EDUCATION GRADIENT\n")
cat("================================================================\n\n")

educ_results <- list()
for (el in levels(df$educ_f)) {
  sub <- df |> filter(educ_f == el)
  n_states <- n_distinct(sub$state_code)
  if (nrow(sub) > 30 & n_states > 2) {
    m <- feols(ln_quantity ~ ln_price + quality_good | state_code,
               data = sub, vcov = "HC1")
    educ_results[[el]] <- list(
      n    = nrow(sub),
      beta = coef(m)["ln_price"],
      se   = sqrt(vcov(m)["ln_price", "ln_price"])
    )
    cat(sprintf("  %-20s: elasticity = %.4f (SE=%.4f), N=%d\n",
                el, educ_results[[el]]$beta, educ_results[[el]]$se,
                educ_results[[el]]$n))
  } else {
    cat(sprintf("  %-20s: too few obs (N=%d)\n", el, nrow(sub)))
  }
}

# Interaction model
cat("\nInteraction model (ln_price × educ_f):\n")
m3_int <- feols(ln_quantity ~ ln_price * educ_f + quality_good + age | state_code,
                data = df, vcov = "HC1")
summary(m3_int)
cat("\n")

################################################################################
# 4. EXTENSIVE vs INTENSIVE MARGIN DECOMPOSITION
################################################################################
cat("================================================================\n")
cat("4. EXTENSIVE vs INTENSIVE MARGIN DECOMPOSITION\n")
cat("================================================================\n\n")

med_q <- median(df$quantity_g)
cat("Median quantity (grams):", med_q, "\n\n")
df <- df |> mutate(large_purchase = as.integer(quantity_g > med_q))

cat("Large purchases:", sum(df$large_purchase),
    " Small:", sum(1 - df$large_purchase), "\n\n")

# 4a. Logit: extensive margin (probability of large purchase)
cat("--- 4a. Logit: Pr(large_purchase) on ln_price ---\n")
m4a <- feglm(large_purchase ~ ln_price + quality_good + age + gender_f + educ_f
             | state_code, data = df, family = binomial(link = "logit"),
             vcov = "HC1")
summary(m4a)

# Average marginal effect of ln_price
# Average marginal effect computed manually
phat <- predict(m4a, type = "response")
ame_lnprice <- mean(coef(m4a)["ln_price"] * phat * (1 - phat))
cat(sprintf("  Avg marginal effect of ln_price on Pr(large): %.4f\n\n", ame_lnprice))

# 4b. OLS on large purchases only (intensive margin)
cat("--- 4b. OLS: intensive margin (large purchases only) ---\n")
m4b <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
             | state_code, data = df |> filter(large_purchase == 1), vcov = "HC1")
summary(m4b)

# 4c. OLS on small purchases only
cat("--- 4c. OLS: small purchases only ---\n")
m4c <- feols(ln_quantity ~ ln_price + quality_good + age + gender_f + educ_f
             | state_code, data = df |> filter(large_purchase == 0), vcov = "HC1")
summary(m4c)

cat("\n--- Margin decomposition summary ---\n")
cat(sprintf("  Extensive margin (logit coef on ln_price): %.4f (SE=%.4f)\n",
            coef(m4a)["ln_price"], sqrt(vcov(m4a)["ln_price","ln_price"])))
cat(sprintf("  Intensive margin (large only, elasticity): %.4f (SE=%.4f), N=%d\n",
            coef(m4b)["ln_price"], sqrt(vcov(m4b)["ln_price","ln_price"]),
            nobs(m4b)))
cat(sprintf("  Small purchases only (elasticity):         %.4f (SE=%.4f), N=%d\n",
            coef(m4c)["ln_price"], sqrt(vcov(m4c)["ln_price","ln_price"]),
            nobs(m4c)))
cat("\n")

################################################################################
# 5. PRICE SEGMENTATION (TERCILES)
################################################################################
cat("================================================================\n")
cat("5. PRICE SEGMENTATION (TERCILES)\n")
cat("================================================================\n\n")

df <- df |> mutate(
  price_tercile = cut(price_mxn_g, breaks = quantile(price_mxn_g, c(0, 1/3, 2/3, 1)),
                      include.lowest = TRUE,
                      labels = c("Cheap", "Medium", "Expensive"))
)

cat("Price tercile cutoffs (MXN/g):\n")
print(quantile(df$price_mxn_g, c(0, 1/3, 2/3, 1)))
cat("\n")

for (seg in c("Cheap", "Medium", "Expensive")) {
  sub <- df |> filter(price_tercile == seg)
  n_states <- n_distinct(sub$state_code)
  if (nrow(sub) > 30 & n_states > 2) {
    m <- feols(ln_quantity ~ ln_price + quality_good + age | state_code,
               data = sub, vcov = "HC1")
    cat(sprintf("  %-10s: elasticity = %.4f (SE=%.4f), N=%d, price range=[%.2f, %.2f]\n",
                seg, coef(m)["ln_price"], sqrt(vcov(m)["ln_price","ln_price"]),
                nrow(sub), min(sub$price_mxn_g), max(sub$price_mxn_g)))
  } else {
    cat(sprintf("  %-10s: too few obs or states (N=%d, states=%d)\n",
                seg, nrow(sub), n_states))
  }
}
cat("\n")

################################################################################
# 6. GEOGRAPHIC HETEROGENEITY (Border / Central / Other)
################################################################################
cat("================================================================\n")
cat("6. GEOGRAPHIC HETEROGENEITY\n")
cat("================================================================\n\n")

cat("Region distribution:\n")
print(table(df$region_f))
cat("\n")

for (reg in c("Border", "Central", "Other")) {
  sub <- df |> filter(region_f == reg)
  n_states <- n_distinct(sub$state_code)
  if (nrow(sub) > 30 & n_states > 1) {
    m <- feols(ln_quantity ~ ln_price + quality_good + age
               | state_code, data = sub, vcov = "HC1")
    cat(sprintf("  %-8s: elasticity = %.4f (SE=%.4f), N=%d, states=%d\n",
                reg, coef(m)["ln_price"], sqrt(vcov(m)["ln_price","ln_price"]),
                nrow(sub), n_states))
  } else {
    cat(sprintf("  %-8s: too few obs (N=%d, states=%d)\n", reg, nrow(sub), n_states))
  }
}

# Interaction model
cat("\nInteraction model (ln_price × region):\n")
m6_int <- feols(ln_quantity ~ ln_price * region_f + quality_good + age
                | state_code, data = df, vcov = "HC1")
summary(m6_int)

cat("\nImplied regional elasticities from interaction model:\n")
b6 <- coef(m6_int)
e_other   <- b6["ln_price"]
e_border  <- b6["ln_price"] + b6["ln_price:region_fBorder"]
e_central <- b6["ln_price"] + b6["ln_price:region_fCentral"]
cat(sprintf("  Other:   %.4f\n", e_other))
cat(sprintf("  Border:  %.4f\n", e_border))
cat(sprintf("  Central: %.4f\n", e_central))

cat("\n\n=== HETEROGENEITY ANALYSIS COMPLETE ===\n")
