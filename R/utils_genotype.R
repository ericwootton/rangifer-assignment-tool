# ============================================================================
# Genotype Conversion Utilities
# ============================================================================

#' Safely extract columns from a data.table as a matrix.
#' Bypasses data.table's [.data.table dispatch which can fail with
#' special characters in column names (e.g. ":" in CHROM:POS IDs)
#' or with single-row tables from dcast.
#' @param dt A data.table
#' @param cols Character vector of column names to extract
#' @return A numeric matrix with named columns
dt_cols_to_matrix <- function(dt, cols) {
  m <- matrix(NA_real_, nrow = nrow(dt), ncol = length(cols))
  colnames(m) <- cols
  for (j in seq_along(cols)) {
    m[, j] <- as.numeric(dt[[cols[j]]])
  }
  m
}

# SRY marker names on the RangiferSNP63k_v2 chip (Y-linked)
SRY_MARKERS <- c("SRY-1", "SRY-2", "SRY-3", "SRY-5", "SRY-6", "SRY-7", "SRY-8")

#' Determine genetic sex from SRY marker calls in FinalReport data
#' Males (XY): SRY probes hybridize to Y chromosome -> valid genotype calls
#' Females (XX): no Y chromosome -> missing calls on all SRY probes
#' @param fr data.table of raw FinalReport with SNPName, clean_id, Allele1Top, Allele2Top
#' @return data.table with sample_id, genetic_sex, sry_calls, sry_total, sex_confidence
determine_sex_from_sry <- function(fr) {
  sry_data <- fr[SNPName %in% SRY_MARKERS]

  if (nrow(sry_data) == 0) {
    return(data.table(
      sample_id = character(),
      genetic_sex = character(),
      sry_calls = integer(),
      sry_total = integer(),
      sex_confidence = character()
    ))
  }

  # A valid call = both alleles are not "-"
  sry_data[, has_call := Allele1Top != "-" & Allele2Top != "-"]

  sex_dt <- sry_data[, .(
    sry_calls = sum(has_call),
    sry_total = .N
  ), by = clean_id]
  setnames(sex_dt, "clean_id", "sample_id")

  # Classification:
  #   >= 5 of 7 SRY called -> Male
  #   <= 1 of 7 SRY called -> Female
  #   2-4 -> Ambiguous
  sex_dt[, genetic_sex := fifelse(
    sry_calls >= 5, "Male",
    fifelse(sry_calls <= 1, "Female", "Ambiguous")
  )]

  sex_dt[, sex_confidence := fifelse(
    (genetic_sex == "Male" & sry_calls >= 6) |
    (genetic_sex == "Female" & sry_calls == 0), "high",
    fifelse(genetic_sex == "Ambiguous", "low", "moderate")
  )]

  sex_dt
}

#' Determine genetic sex from X-chromosome heterozygosity ratio
#' Works for both FinalReport and 012 matrix inputs.
#' Males (XY) are hemizygous on X -> very low X het rate (ratio ~0.06-0.31)
#' Females (XX) are diploid on X -> normal X het rate (ratio ~0.85-1.06)
#' @param geno_012 data.table with sample_id + SNP columns (012-coded)
#' @param x_snp_ids Character vector of X-linked SNP IDs (CHROM:POS format)
#' @return data.table with sample_id, x_het_ratio, x_het_sex, x_het_confidence
determine_sex_from_xhet <- function(geno_012, x_snp_ids) {
  all_snps <- colnames(geno_012)[-1]
  x_in <- intersect(x_snp_ids, all_snps)
  auto_in <- setdiff(all_snps, x_snp_ids)

  if (length(x_in) < 50) {
    # Too few X-linked SNPs for reliable ratio
    return(data.table(
      sample_id = geno_012$sample_id,
      x_het_ratio = NA_real_,
      x_het_sex = NA_character_,
      x_het_confidence = NA_character_,
      n_x_snps = length(x_in)
    ))
  }

  x_mat <- dt_cols_to_matrix(geno_012, x_in)
  a_mat <- dt_cols_to_matrix(geno_012, auto_in)

  x_called <- rowSums(!is.na(x_mat))
  a_called <- rowSums(!is.na(a_mat))
  x_het <- ifelse(x_called > 0, rowSums(x_mat == 1, na.rm = TRUE) / x_called, NA_real_)
  a_het <- ifelse(a_called > 0, rowSums(a_mat == 1, na.rm = TRUE) / a_called, NA_real_)
  ratio <- ifelse(a_het > 0, x_het / a_het, NA_real_)

  # Classification thresholds (validated: males 0.06-0.31, females 0.85-1.06)
  sex <- fifelse(
    is.na(ratio), NA_character_,
    fifelse(ratio < 0.35, "Male",
      fifelse(ratio >= 0.70, "Female", "Ambiguous")
    )
  )

  confidence <- fifelse(
    is.na(ratio), NA_character_,
    fifelse(
      (sex == "Male" & ratio < 0.30) | (sex == "Female" & ratio > 0.75), "high",
      fifelse(sex == "Ambiguous", "low", "moderate")
    )
  )

  data.table(
    sample_id = geno_012$sample_id,
    x_het_ratio = round(ratio, 4),
    x_het_sex = sex,
    x_het_confidence = confidence,
    n_x_snps = length(x_in)
  )
}

complement_base <- function(base) {
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  comp[base]
}

#' Convert allele pairs to 012 coding for a single SNP
#' @param allele_pairs Character vector of 2-char allele pairs ("AG", "CC", etc.)
#' @param ref Reference allele
#' @param alt Alternate allele
#' @return Integer vector: 0 (hom ref), 1 (het), 2 (hom alt), NA (missing/failed)
alleles_to_012 <- function(allele_pairs, ref, alt) {
  a1 <- substr(allele_pairs, 1, 1)
  a2 <- substr(allele_pairs, 2, 2)

  # Skip strand-ambiguous A/T or C/G SNPs
  is_ambiguous <- (ref == "A" & alt == "T") | (ref == "T" & alt == "A") |
                  (ref == "C" & alt == "G") | (ref == "G" & alt == "C")
  if (is_ambiguous) return(rep(NA_integer_, length(allele_pairs)))

  # Gather observed alleles for strand detection
  non_na <- which(!is.na(allele_pairs) & allele_pairs != "" & allele_pairs != "NA")
  if (length(non_na) == 0) return(rep(NA_integer_, length(allele_pairs)))

  sample_check <- non_na[seq_len(min(20, length(non_na)))]
  all_alleles <- unique(c(a1[sample_check], a2[sample_check]))
  all_alleles <- all_alleles[!is.na(all_alleles) & all_alleles != ""]

  use_complement <- FALSE
  if (all(all_alleles %in% c(ref, alt))) {
    # Direct match
  } else if (all(complement_base(all_alleles) %in% c(ref, alt))) {
    use_complement <- TRUE
    a1 <- complement_base(a1)
    a2 <- complement_base(a2)
  } else {
    return(rep(NA_integer_, length(allele_pairs)))
  }

  # Count ALT alleles
  result <- ifelse(
    is.na(a1) | is.na(a2) | a1 == "" | a2 == "",
    NA_integer_,
    as.integer((a1 == alt) + (a2 == alt))
  )

  # Validate: alleles should be ref or alt
  invalid <- !is.na(a1) & !is.na(a2) & a1 != "" & a2 != "" &
             !(a1 %in% c(ref, alt) & a2 %in% c(ref, alt))
  result[invalid] <- NA_integer_

  result
}

#' Parse Illumina FinalReport CSV and convert to 012 genotype matrix
#' @param file_path Path to FinalReport CSV
#' @param chip_snp_map data.table with snp_name, snp_id columns
#' @param ref_lookup Named vector: snp_id -> REF allele
#' @param alt_lookup Named vector: snp_id -> ALT allele
#' @return list(geno_012 = data.table, stats = list)
parse_finalreport <- function(file_path, chip_snp_map, ref_lookup, alt_lookup) {
  # Find [Data] section
  header_lines <- readLines(file_path, n = 30)
  skip_n <- grep("^\\[Data\\]", header_lines)
  if (length(skip_n) == 0) stop("Cannot find [Data] section in FinalReport")

  fr <- data.table::fread(file_path, skip = skip_n, header = TRUE)

  # Normalize column names
  if ("Sample ID" %in% colnames(fr)) setnames(fr, "Sample ID", "SampleID")
  if ("SNP Name" %in% colnames(fr)) setnames(fr, "SNP Name", "SNPName")
  if ("Allele1 - Top" %in% colnames(fr)) setnames(fr, "Allele1 - Top", "Allele1Top")
  if ("Allele2 - Top" %in% colnames(fr)) setnames(fr, "Allele2 - Top", "Allele2Top")

  # Extract clean sample ID
  fr[, clean_id := sapply(strsplit(as.character(SampleID), "#"), function(x) {
    if (length(x) >= 2) x[2] else x[1]
  })]

  # Sex determination from SRY markers (before filtering to mapped SNPs)
  sex_result <- determine_sex_from_sry(fr)

  # Map SNP names to CHROM:POS
  fr <- merge(fr, chip_snp_map, by.x = "SNPName", by.y = "snp_name", all.x = TRUE)
  n_unmapped <- sum(is.na(fr$snp_id))
  fr <- fr[!is.na(snp_id)]

  # Encode allele pair
  fr[, allele_pair := ifelse(
    Allele1Top == "-" | Allele2Top == "-",
    NA_character_,
    paste0(Allele1Top, Allele2Top)
  )]

  # Deduplicate (multiple probe designs for same position)
  fr <- unique(fr, by = c("clean_id", "snp_id"))

  # Pivot to wide
  allele_wide <- dcast(fr, clean_id ~ snp_id, value.var = "allele_pair")
  setnames(allele_wide, "clean_id", "sample_id")

  # Convert allele pairs to 012
  snp_cols <- colnames(allele_wide)[-1]
  shared_snps <- intersect(snp_cols, names(ref_lookup))

  n_direct <- 0L; n_complement <- 0L; n_failed <- 0L; n_ambiguous <- 0L

  # Build as a plain list first, then convert to data.table once at the end.
  # Adding 46k+ columns one-by-one to a data.table exceeds its internal
  # column allocation (truelength) and corrupts the object.
  geno_list <- vector("list", length(shared_snps) + 1L)
  names(geno_list) <- c("sample_id", shared_snps)
  geno_list[["sample_id"]] <- allele_wide$sample_id

  for (snp in shared_snps) {
    ref <- ref_lookup[snp]
    alt <- alt_lookup[snp]
    ap <- allele_wide[[snp]]

    coded <- alleles_to_012(ap, ref, alt)
    geno_list[[snp]] <- coded

    # Track stats
    non_na <- which(!is.na(ap) & ap != "" & ap != "NA")
    if (length(non_na) == 0) next

    is_amb <- (ref == "A" & alt == "T") | (ref == "T" & alt == "A") |
              (ref == "C" & alt == "G") | (ref == "G" & alt == "C")
    if (is_amb) { n_ambiguous <- n_ambiguous + 1L; next }

    sample_alleles <- unique(c(
      substr(ap[non_na[seq_len(min(20, length(non_na)))]], 1, 1),
      substr(ap[non_na[seq_len(min(20, length(non_na)))]], 2, 2)
    ))
    sample_alleles <- sample_alleles[!is.na(sample_alleles) & sample_alleles != ""]

    if (all(sample_alleles %in% c(ref, alt))) {
      n_direct <- n_direct + 1L
    } else if (all(complement_base(sample_alleles) %in% c(ref, alt))) {
      n_complement <- n_complement + 1L
    } else {
      n_failed <- n_failed + 1L
    }
  }

  geno_012 <- as.data.table(geno_list)

  stats <- list(
    n_samples = nrow(geno_012),
    n_snps_mapped = length(shared_snps),
    n_unmapped = n_unmapped,
    n_direct = n_direct,
    n_complement = n_complement,
    n_ambiguous = n_ambiguous,
    n_failed = n_failed
  )

  list(geno_012 = geno_012, stats = stats, sex = sex_result)
}

#' Compute per-sample QC metrics
#' @param geno_012 data.table with sample_id + SNP columns (012-coded)
#' @return data.table with sample_id, call_rate, het_rate, quality_flags
compute_sample_qc <- function(geno_012) {
  snp_cols <- colnames(geno_012)[-1]
  n_snps <- length(snp_cols)

  geno_mat <- dt_cols_to_matrix(geno_012, snp_cols)

  call_rate <- rowSums(!is.na(geno_mat)) / n_snps
  het_rate <- rowSums(geno_mat == 1, na.rm = TRUE) / rowSums(!is.na(geno_mat))

  flags <- character(nrow(geno_012))
  flags[call_rate < 0.70] <- "low_call_rate"
  flags[het_rate > 0.50] <- ifelse(
    flags[het_rate > 0.50] == "", "possible_mixed",
    paste0(flags[het_rate > 0.50], ";possible_mixed")
  )
  flags[het_rate < 0.15] <- ifelse(
    flags[het_rate < 0.15] == "", "possible_non_rangifer",
    paste0(flags[het_rate < 0.15], ";possible_non_rangifer")
  )

  data.table(
    sample_id = geno_012$sample_id,
    call_rate = round(call_rate, 4),
    het_rate = round(het_rate, 4),
    quality_flags = ifelse(flags == "", "pass", flags)
  )
}
