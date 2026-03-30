# ============================================================================
# Rangifer Assignment Shiny App — Global Setup
# ============================================================================

library(shiny)
library(data.table)
library(adegenet)
library(e1071)
library(jsonlite)
library(DT)
library(plotly)
library(bslib)
library(shinyjs)

# Remove upload size limit (default is 5 MB)
options(shiny.maxRequestSize = Inf)

# Source utilities and modules
source("R/utils_genotype.R")
source("R/utils_hierarchy.R")
source("R/mod_classify.R")

# ---- Load hierarchy ----
hierarchy <- fromJSON("data/hierarchy_definition.json")
routing <- build_routing_table(hierarchy)

# Build display-name mapping for "Not applicable" ecotypes
# These are monotypic subspecies where the ecotype level is meaningless;
# replace with the actual herd name from L3 for cleaner display.
ecotype_display_names <- list()
for (l1_group in names(hierarchy$level2)) {
  l2_info <- hierarchy$level2[[l1_group]]
  groups <- if (is.list(l2_info$groups)) unlist(l2_info$groups) else l2_info$groups
  if (length(groups) == 1 && groups == "Not applicable") {
    # Find the corresponding L3 herd name
    l3_key <- paste0(gsub(" ", "_", l1_group), "_Not_applicable")
    l3_info <- hierarchy$level3[[l3_key]]
    if (!is.null(l3_info)) {
      herd <- if (is.list(l3_info$groups)) unlist(l3_info$groups) else l3_info$groups
      ecotype_display_names[[l1_group]] <- herd[1]
    }
  }
}
# Result: e.g., "Rangifer tarandus fennicus" -> "Finnish Forest",
#                "Rangifer tarandus pearyi" -> "Peary"

# ---- Load SNP mapping data ----
chip_snp_map <- readRDS("data/chip_snp_map.rds")
snp_info <- readRDS("data/snp_info.rds")
ref_lookup <- setNames(snp_info$ref, snp_info$snp_id)
alt_lookup <- setNames(snp_info$alt, snp_info$snp_id)

# ---- Load X-chromosome SNP list ----
x_snp_ids <- if (file.exists("data/x_chromosome_snps.rds")) {
  readRDS("data/x_chromosome_snps.rds")
} else {
  character(0)
}
cat("X-chromosome SNPs:", length(x_snp_ids), "\n")

# ---- Load DAPC models ----
dapc_files <- list.files("data/models/dapc", pattern = "\\.rds$", full.names = TRUE)
dapc_models <- setNames(
  lapply(dapc_files, readRDS),
  tools::file_path_sans_ext(basename(dapc_files))
)
cat("Loaded", length(dapc_models), "DAPC models\n")

# ---- Load SVM models ----
svm_files <- list.files("data/models/svm", pattern = "\\.rds$", full.names = TRUE)
svm_models <- if (length(svm_files) > 0) {
  setNames(
    lapply(svm_files, readRDS),
    tools::file_path_sans_ext(basename(svm_files))
  )
} else {
  list()
}
cat("Loaded", length(svm_models), "SVM models\n")

# ---- Load popfinder NN weights ----
nn_files <- list.files("data/models/popfinder", pattern = "\\.rds$", full.names = TRUE)
nn_weights <- if (length(nn_files) > 0) {
  setNames(
    lapply(nn_files, readRDS),
    tools::file_path_sans_ext(basename(nn_files))
  )
} else {
  list()
}
cat("Loaded", length(nn_weights), "popfinder NN weight sets\n")

# ---- Load node metadata ----
node_files <- list.files("data/node_data", pattern = "\\.rds$", full.names = TRUE)
node_data <- setNames(
  lapply(node_files, readRDS),
  tools::file_path_sans_ext(basename(node_files))
)
cat("Loaded", length(node_data), "node metadata objects\n")

# ---- Summary ----
cat("\n=== Rangifer Assignment App Ready ===\n")
cat("Hierarchy: 3 levels,", length(hierarchy$level1$groups), "L1 groups\n")
cat("SNP map:", nrow(chip_snp_map), "entries\n")
cat("Models: DAPC=", length(dapc_models), " SVM=", length(svm_models),
    " NN=", length(nn_weights), "\n")
