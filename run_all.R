#!/usr/bin/env Rscript

###############################################################################
# Master replication script
# Price Elasticity of Demand for Marijuana in Mexico
#
# Run from the project root:
#   Rscript run_all.R
#
# Useful flags:
#   RUN_BOOTSTRAP=0     skip bootstrap inference during quick checks
#   COMPILE_PDF=0       regenerate tables without compiling the paper
#   OSRM_ALLOW_NETWORK=1 allow OSRM queries if a cached matrix is missing
###############################################################################

find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "PriceElasticityDemandMarijuanaMexico.Rproj"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find project root. Run this script inside the project.", call. = FALSE)
    }
    path <- parent
  }
}

root <- find_project_root()
setwd(root)

project_lib <- file.path(root, ".Rlib")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

flag <- function(name, default = TRUE) {
  value <- Sys.getenv(name, unset = if (isTRUE(default)) "1" else "0")
  value %in% c("1", "true", "TRUE", "yes", "YES")
}

run_bootstrap <- flag("RUN_BOOTSTRAP", TRUE)
compile_pdf <- flag("COMPILE_PDF", TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_dir <- file.path(root, "Logs", paste0("replication_", timestamp))
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

step_counter <- 0L
started_at <- Sys.time()

message("Project root: ", root)
message("Logs: ", log_dir)

check_packages <- function(pkgs) {
  missing <- pkgs[
    !vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]
  if (length(missing) > 0) {
    stop(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      "\nInstall them before running the replication pipeline.",
      call. = FALSE
    )
  }
}

required_packages <- c(
  "AER", "MASS", "arrow", "dplyr", "fixest", "fwildclusterboot",
  "here", "kableExtra", "knitr", "lubridate",
  "modelsummary", "readr", "sf", "stringr",
  "tidyr"
)

check_packages(unique(required_packages))

check_command <- function(command) {
  if (!nzchar(Sys.which(command))) {
    stop("Required command not found on PATH: ", command, call. = FALSE)
  }
}

if (compile_pdf) check_command("latexmk")

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

run_command <- function(label, command, args = character(), workdir = root) {
  step_counter <<- step_counter + 1L
  log_path <- file.path(log_dir, sprintf("%02d_%s.log", step_counter, safe_name(label)))
  message(sprintf("\n[%02d] %s", step_counter, label))
  message("     command: ", paste(c(command, args), collapse = " "))
  message("     log: ", log_path)

  old <- setwd(workdir)
  on.exit(setwd(old), add = TRUE)

  t0 <- Sys.time()
  status <- system2(command, args = args, stdout = log_path, stderr = log_path)
  elapsed <- difftime(Sys.time(), t0, units = "secs")

  if (!identical(status, 0L)) {
    stop(
      "Step failed: ", label,
      "\nExit status: ", status,
      "\nSee log: ", log_path,
      call. = FALSE
    )
  }

  message(sprintf("     done in %.1f seconds", as.numeric(elapsed)))
  invisible(log_path)
}

rscript <- file.path(R.home("bin"), "Rscript")

run_qmd <- function(label, qmd_path) {
  runner <- file.path(log_dir, paste0(safe_name(label), "_runner.R"))
  writeLines(
    c(
      "tmp <- tempfile(fileext = '.R')",
      paste0(
        "knitr::purl(", deparse(qmd_path), ", output = tmp, ",
        "documentation = 0, quiet = TRUE)"
      ),
      "source(tmp, echo = FALSE)"
    ),
    runner
  )
  run_command(label, rscript, runner)
}

run_qmd("Clean transaction data", "Code/Cleaning.qmd")
run_qmd("Descriptive statistics", "Code/DescriptiveAnalysis.qmd")

paper_scripts <- c(
  "Code/RegressionAnalysis.R",
  "Code/IV_SeizureGravity.R",
  "Code/RegressionExtensions.R",
  "Code/MainIVDiagnostics.R",
  "Code/FunctionalForms_Bounds.R",
  "Code/ExternalDataAppendix.R"
)

if (run_bootstrap) {
  paper_scripts <- append(
    paper_scripts,
    "Code/BootstrapInference.R",
    after = match("Code/RegressionExtensions.R", paper_scripts)
  )
}

for (script in paper_scripts) {
  run_command(basename(script), rscript, script)
}

if (compile_pdf) {
  run_command(
    "Compile manuscript PDF",
    "latexmk",
    c("-pdf", "-interaction=nonstopmode", "-halt-on-error",
      "PriceElasticityDemandMarijuanaMexico.tex"),
    workdir = file.path(root, "LaTeX")
  )
}

elapsed_total <- difftime(Sys.time(), started_at, units = "mins")
message("\nReplication completed successfully.")
message(sprintf("Total runtime: %.1f minutes", as.numeric(elapsed_total)))
message("Final manuscript: ", file.path(root, "LaTeX", "PriceElasticityDemandMarijuanaMexico.pdf"))
