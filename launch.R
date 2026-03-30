#!/usr/bin/env Rscript
# ============================================================================
# Rangifer Assignment Tool — Cross-Platform Launcher
# ============================================================================
# Checks for required R packages, installs any that are missing, and launches
# the Shiny app. Works on Mac, Windows, and Linux.
#
# Usage:
#   Rscript launch.R
#   Or double-click launch.command (Mac) / launch.bat (Windows)
# ============================================================================

cat("\n")
cat("===================================================\n")
cat("  Rangifer Assignment Tool\n")
cat("  Probabilistic caribou subspecies/ecotype/herd\n")
cat("  assignment from SNP chip genotype data\n")
cat("===================================================\n\n")

# ---- Set working directory to script location ----
# Works whether called from Rscript, RStudio, or sourced
get_script_dir <- function() {
  # Method 1: Rscript --file= argument
  args <- commandArgs(FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("--file=", "", file_arg[1]))))
  }
  # Method 2: RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    return(normalizePath(dirname(rstudioapi::getSourceEditorContext()$path)))
  }
  # Method 3: sys.frame (sourced)
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)
    if (!is.null(f$ofile)) {
      return(normalizePath(dirname(f$ofile)))
    }
  }
  # Fallback: current directory
  return(normalizePath("."))
}

app_dir <- get_script_dir()
setwd(app_dir)
cat("App directory:", app_dir, "\n\n")

# ---- Check R version ----
rv <- getRversion()
if (rv < "4.1.0") {
  stop("R >= 4.1.0 is required (you have ", as.character(rv), ")")
}
cat("R version:", as.character(rv), "\n")

# ---- Check and install dependencies ----
required_packages <- c(
  "shiny",       # Web framework
  "data.table",  # Fast data manipulation
  "adegenet",    # DAPC models (predict.dapc)
  "e1071",       # SVM models (predict.svm)
  "jsonlite",    # JSON parsing
  "DT",          # Interactive tables
  "plotly",      # Interactive plots
  "bslib",       # Bootstrap theming
  "shinyjs"      # Enable/disable UI elements
)

cat("\nChecking packages...\n")
missing <- character()
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("  ", pkg, "OK\n")
  } else {
    cat("  ", pkg, "MISSING\n")
    missing <- c(missing, pkg)
  }
}

if (length(missing) > 0) {
  cat("\nInstalling", length(missing), "missing package(s):",
      paste(missing, collapse = ", "), "\n")
  cat("This may take a few minutes on first run...\n\n")

  install.packages(missing, repos = "https://cloud.r-project.org",
                   dependencies = TRUE, quiet = FALSE)

  # Verify installation
  still_missing <- character()
  for (pkg in missing) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      still_missing <- c(still_missing, pkg)
    }
  }
  if (length(still_missing) > 0) {
    stop("Failed to install: ", paste(still_missing, collapse = ", "),
         "\nPlease install these manually with: ",
         'install.packages(c("', paste(still_missing, collapse = '", "'), '"))')
  }
  cat("All packages installed successfully.\n")
}

# ---- Check data files exist ----
cat("\nChecking data files...\n")
required_files <- c(
  "data/hierarchy_definition.json",
  "data/chip_snp_map.rds",
  "data/snp_info.rds",
  "data/x_chromosome_snps.rds"
)

model_dirs <- c("data/models/dapc", "data/node_data")
for (d in model_dirs) {
  if (!dir.exists(d) || length(list.files(d)) == 0) {
    stop("Missing data directory: ", d,
         "\nRun prep/prepare_app_data.R first to bundle pipeline results.")
  }
}
for (f in required_files) {
  if (!file.exists(f)) {
    stop("Missing data file: ", f,
         "\nRun prep/prepare_app_data.R first to bundle pipeline results.")
  }
}
cat("  Data files OK\n")

# ---- Launch ----
cat("\n")
cat("===================================================\n")
cat("  Starting Rangifer Assignment Tool...\n")
cat("  The app will open in your default web browser.\n")
cat("  Press Ctrl+C (or Esc) in this window to stop.\n")
cat("===================================================\n\n")

shiny::runApp(app_dir, launch.browser = TRUE)
