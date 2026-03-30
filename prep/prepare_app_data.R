#!/usr/bin/env Rscript
# ============================================================================
# Prepare App Data
# ============================================================================
# One-time script to bundle pipeline results into app-ready RDS/JSON files.
# This is used during development to rebuild the data/ directory from the
# Nextflow pipeline output. End users do not need to run this script.
#
# Usage: Rscript prep/prepare_app_data.R /path/to/pipeline/base
# ============================================================================

library(data.table)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript prep/prepare_app_data.R /path/to/pipeline/base")
}
base_dir <- normalizePath(args[1])

results_dir <- file.path(base_dir, "results/caribou_results/rangifer_assignment_output")
app_data_dir <- file.path(base_dir, "data")

cat("Base dir:", base_dir, "\n")
cat("Results dir:", results_dir, "\n")
cat("App data dir:", app_data_dir, "\n\n")

# Ensure output dirs exist
dir.create(file.path(app_data_dir, "models/dapc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(app_data_dir, "models/svm"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(app_data_dir, "models/popfinder"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(app_data_dir, "node_data"), recursive = TRUE, showWarnings = FALSE)

# ---- 1. Hierarchy definition ----
cat("1. Copying hierarchy_definition.json...\n")
file.copy(
  file.path(results_dir, "structure_discovery/hierarchy_definition.json"),
  file.path(app_data_dir, "hierarchy_definition.json"),
  overwrite = TRUE
)

# ---- 2. Chip SNP map ----
cat("2. Building chip_snp_map.rds...\n")
chip_map <- fread(
  file.path(results_dir, "data_integration/chip_snp_map.tsv"),
  header = FALSE, col.names = c("snp_name", "chrom", "pos", "snp_id")
)
chip_map <- unique(chip_map[, .(snp_name, snp_id)], by = "snp_name")
saveRDS(chip_map, file.path(app_data_dir, "chip_snp_map.rds"))
cat("  ", nrow(chip_map), "unique SNP name mappings\n")

# ---- 3. SNP info (REF/ALT lookup) ----
cat("3. Building snp_info.rds...\n")
snp_info <- fread(file.path(results_dir, "data_integration/snp_info.tsv"), header = TRUE)
saveRDS(snp_info, file.path(app_data_dir, "snp_info.rds"))
cat("  ", nrow(snp_info), "SNPs with REF/ALT\n")

# ---- 4. Load hierarchy for node enumeration ----
hierarchy <- fromJSON(file.path(app_data_dir, "hierarchy_definition.json"))

# Build list of all node names
all_nodes <- c("root")
for (l1_group in names(hierarchy$level2)) {
  node_info <- hierarchy$level2[[l1_group]]
  all_nodes <- c(all_nodes, node_info$node_name)
}
for (l2_key in names(hierarchy$level3)) {
  node_info <- hierarchy$level3[[l2_key]]
  all_nodes <- c(all_nodes, node_info$node_name)
}
cat("\n4. Processing", length(all_nodes), "nodes...\n")

# Calibrate SVM probability thresholds from leave-one-out cross-validation.
# Upper = 5th percentile of correct-prediction probabilities (floor 0.70).
# Lower = 20th percentile (floor 0.50).
calibrate_threshold <- function(loocv_file) {
  if (!file.exists(loocv_file)) return(list(upper = 0.95, lower = 0.80))
  loocv <- fread(loocv_file, header = TRUE)
  correct_probs <- loocv[correct == TRUE]$posterior_max
  if (length(correct_probs) < 3) return(list(upper = 0.95, lower = 0.80))
  upper <- max(quantile(correct_probs, 0.05, na.rm = TRUE), 0.70)
  lower <- max(quantile(correct_probs, 0.20, na.rm = TRUE), 0.50)
  list(upper = as.numeric(upper), lower = as.numeric(lower))
}

node_data_dir <- file.path(results_dir, "classification/node_data/nodes")
dapc_dir <- file.path(results_dir, "classification/dapc")
assignpop_dir <- file.path(results_dir, "classification/assignpop")
popfinder_dir <- file.path(results_dir, "classification/popfinder")

for (node_name in all_nodes) {
  cat("  Processing node:", node_name, "\n")

  # --- Node metadata ---
  node_path <- file.path(node_data_dir, node_name)
  if (!dir.exists(node_path)) {
    cat("    WARNING: node data dir not found, skipping\n")
    next
  }

  ref_geno <- fread(file.path(node_path, "reference_genotypes.tsv"), header = TRUE)
  ref_labels <- fread(file.path(node_path, "reference_labels.tsv"), header = TRUE)
  meta <- fromJSON(file.path(node_path, "node_meta.json"))

  geno_mat <- as.matrix(ref_geno[, -1])
  snp_names <- colnames(ref_geno)[-1]
  ref_means <- colMeans(geno_mat, na.rm = TRUE)
  groups <- if (is.list(meta$groups)) unlist(meta$groups) else meta$groups
  n_groups <- length(groups)

  # DAPC LOOCV thresholds
  dapc_results_dir <- file.path(dapc_dir, paste0(node_name, "_dapc_results"))
  thresholds <- calibrate_threshold(file.path(dapc_results_dir, "loocv_results.tsv"))

  # DAPC accuracy
  dapc_summary_file <- file.path(dapc_results_dir, "summary.tsv")
  loocv_accuracy <- NA_real_
  if (file.exists(dapc_summary_file)) {
    s <- fread(dapc_summary_file, header = TRUE)
    loocv_accuracy <- s$loocv_accuracy[1]
  }

  # SVM thresholds (from assignPOP SVM LOOCV)
  ap_dir_name <- paste0(node_name, "_fst_filtered_assignpop_results")
  ap_path <- file.path(assignpop_dir, ap_dir_name)
  svm_thresholds <- list(upper = 0.95, lower = 0.80)
  if (dir.exists(ap_path)) {
    svm_loocv_file <- file.path(ap_path, "loocv_svm.tsv")
    if (file.exists(svm_loocv_file)) {
      svm_loocv <- fread(svm_loocv_file, header = TRUE)
      correct_probs <- svm_loocv[correct == TRUE]$probability
      if (length(correct_probs) >= 3) {
        svm_thresholds$upper <- as.numeric(max(quantile(correct_probs, 0.05, na.rm = TRUE), 0.70))
        svm_thresholds$lower <- as.numeric(max(quantile(correct_probs, 0.20, na.rm = TRUE), 0.50))
      }
    }
  }

  # Popfinder thresholds
  pf_dir_name <- paste0(node_name, "_popfinder_results")
  pf_path <- file.path(popfinder_dir, pf_dir_name)
  pf_thresholds <- list(upper = 0.95, lower = 0.80)
  if (dir.exists(pf_path)) {
    pf_thresholds <- calibrate_threshold(file.path(pf_path, "loocv_results.tsv"))
  }

  node_obj <- list(
    node_name = node_name,
    snp_names = snp_names,
    ref_means = ref_means,
    groups = groups,
    n_groups = n_groups,
    level = meta$level,
    loocv_accuracy = loocv_accuracy,
    thresholds = thresholds,
    svm_thresholds = svm_thresholds,
    pf_thresholds = pf_thresholds
  )
  saveRDS(node_obj, file.path(app_data_dir, "node_data", paste0(node_name, ".rds")))

  # --- DAPC model (multi-group nodes only) ---
  if (n_groups >= 2) {
    dapc_model_file <- file.path(dapc_results_dir, "dapc_model.rds")
    if (file.exists(dapc_model_file)) {
      file.copy(dapc_model_file,
                file.path(app_data_dir, "models/dapc", paste0(node_name, ".rds")),
                overwrite = TRUE)
      cat("    Copied DAPC model\n")
    }
  }

  # --- SVM model (multi-group nodes only) ---
  if (n_groups >= 2 && dir.exists(ap_path)) {
    svm_model_file <- file.path(ap_path, "svm_model.rds")
    svm_snp_file <- file.path(ap_path, "svm_snp_names.rds")
    svm_means_file <- file.path(ap_path, "svm_ref_means.rds")
    if (file.exists(svm_model_file)) {
      file.copy(svm_model_file,
                file.path(app_data_dir, "models/svm", paste0(node_name, ".rds")),
                overwrite = TRUE)
      cat("    Copied SVM model\n")
    } else {
      cat("    WARNING: svm_model.rds not found (pipeline re-run needed)\n")
    }
    if (file.exists(svm_snp_file)) {
      node_obj$svm_snp_names <- readRDS(svm_snp_file)
    }
    if (file.exists(svm_means_file)) {
      node_obj$svm_ref_means <- readRDS(svm_means_file)
    }
    # Re-save node_obj with SVM metadata
    saveRDS(node_obj, file.path(app_data_dir, "node_data", paste0(node_name, ".rds")))
  }

  # --- Popfinder NN weights (multi-group nodes only) ---
  if (n_groups >= 2 && dir.exists(pf_path)) {
    nn_file <- file.path(pf_path, "nn_weights.json")
    if (file.exists(nn_file)) {
      file.copy(nn_file,
                file.path(app_data_dir, "models/popfinder", paste0(node_name, ".json")),
                overwrite = TRUE)
      cat("    Copied popfinder NN weights\n")
    } else {
      cat("    WARNING: nn_weights.json not found (pipeline re-run needed)\n")
    }
  }
}

cat("\nData preparation complete.\n")
cat("App data directory:", app_data_dir, "\n")
