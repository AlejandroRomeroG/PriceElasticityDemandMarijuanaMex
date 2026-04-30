# The Price Elasticity of Demand for Marijuana in Mexico

This repository contains the replication package and first draft of an academic paper estimating the transaction-level price elasticity of marijuana demand in Mexico.

The main manuscript is compiled from `LaTeX/PriceElasticityDemandMarijuanaMexico.tex`. The full reproducibility entry point is:

```bash
Rscript run_all.R
```

The master script cleans the raw survey data, regenerates the descriptive and econometric tables, rebuilds appendix diagnostics, and compiles the manuscript PDF.

## Research Design

The paper uses an anonymous crowdsourced survey of marijuana purchases in Mexico fielded in 2019. The cleaned analytic sample contains 3,293 completed transactions across all 32 states and 351 municipalities.

The preferred specification is a log-log demand regression with product-quality controls, demographics, and state-by-month fixed effects. The benchmark coefficient is interpreted as a transaction-level elasticity: it describes how the quantity purchased in an observed transaction varies with the recorded unit price. It does not directly identify participation, purchase frequency, total individual consumption, or legal-market demand after regulation.

## Repository Structure

```text
Code/      Cleaning, table generation, econometric scripts, and OSRM helpers
Data/      Raw survey data, cleaned parquet, MUCD files, OSRM caches, INEGI layers
LaTeX/     Manuscript source, bibliography, generated tables, and compiled PDF
Literature/Reference PDFs used while drafting
run_all.R  Master replication script
```

## Data Inputs

The main survey source is:

Bejarano Romero, Raul; Alonso Arechar, Antonio; Enciso Higuera, Froylán (2020), "Dataset on marijuana prices in Mexico", Mendeley Data, V1. DOI: `10.17632/hk4vydrg4v.1`.

The repository includes the UTF-8 encoded raw survey file used by the cleaning script:

```text
Data/Marijuana_Prices_in_Mexico_raw-3_utf8_bom.csv
```

External inputs used for IV diagnostics are:

- `Data/MUCD/`: municipality-month marijuana seizure records from Mexico Unido contra la Delincuencia.
- `Data/Derived/`: cached OSRM/OpenStreetMap driving-time matrices used by the IV scripts.
- `Data/Shapefiles/conjunto_de_datos/00mun.*` and `00ent.*`: INEGI municipality and state geometries used to compute centroids.

The OSRM caches are included so the main replication can run without querying a public routing server. If a cache is missing and you want to rebuild it, set `OSRM_ALLOW_NETWORK=1` before running the relevant scripts. This may be slow and depends on OSRM server availability.

## Software Requirements

The pipeline uses R and LaTeX. The cleaning and descriptive files are Quarto notebooks, but the master script executes their code chunks directly through `knitr::purl()` so the replication does not depend on the Quarto CLI.

- R with packages: `arrow`, `dplyr`, `fixest`, `modelsummary`, `sf`, `readr`, `tidyr`, `lubridate`, `stringr`, `AER`, `fwildclusterboot`, `kableExtra`, `knitr`.
- `latexmk` with `biber` for compiling the manuscript bibliography.

The scripts do not install packages automatically. If a local `.Rlib/` directory exists, `run_all.R` prepends it to `.libPaths()`.

## Reproducing the Draft

Run the full pipeline from the repository root:

```bash
Rscript run_all.R
```

Expected key outputs:

- `Data/Marijuana_Prices_in_Mexico_clean.parquet`
- `LaTeX/tables/*.tex`
- `LaTeX/PriceElasticityDemandMarijuanaMexico.pdf`
- step-by-step logs under `Logs/replication_<timestamp>/`

Useful flags:

```bash
RUN_BOOTSTRAP=0 Rscript run_all.R
COMPILE_PDF=0 Rscript run_all.R
RUN_EXPLORATORY=1 Rscript run_all.R
OSRM_ALLOW_NETWORK=1 Rscript run_all.R
```

`RUN_EXPLORATORY=1` also runs exploratory scripts that are not part of the current paper tables. The default pipeline is the reproducible path for the first draft manuscript.
