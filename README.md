# Price Elasticity of Demand for Marijuana in Mexico

This repository contains a fully reproducible workflow to clean, explore, and analyze a transaction-level dataset on marijuana purchases in Mexico, with the goal of estimating the **price elasticity of demand**.

All stages of the pipeline (**cleaning**, **exploratory analysis**, and **econometric estimation**) are implemented in **Quarto notebooks using R**.

---

## Project overview

**Research question:** How sensitive is marijuana demand to changes in price in Mexico?

**Core outputs (planned):**
- Cleaned analytic dataset(s)
- Exploratory plots and descriptive tables
- Econometric estimates of price elasticity (baseline + robustness checks)

---

## Data

### Source
Bejarano Romero, Raul; Alonso Arechar, Antonio; Enciso Higuera, Froylán (2020), **“Dataset on marijuana prices in Mexico”**, Mendeley Data, V1. DOI: 10.17632/hk4vydrg4v.1

### Raw file in this repository
- `Data/Marijuana_Prices_in_Mexico_raw-3_utf8_bom.csv`

**Important note:** The only modification applied to the raw dataset prior to analysis was converting the file **encoding to UTF-8** (to ensure consistent parsing of accents and special characters). No rows, values, or variables were altered outside the scripted cleaning steps documented in the notebooks.
