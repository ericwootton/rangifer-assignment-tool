# ============================================================================
# Hierarchical 3-Method Classification Engine
# ============================================================================

#' Predict at a single DAPC node
#' @param node_name Node identifier
#' @param geno_012 data.table with sample_id + SNP columns
#' @param nd Node data (snp_names, ref_means, groups, etc.)
#' @param model DAPC model object
#' @return list(assignment, max_posterior, posteriors)
predict_dapc <- function(node_name, geno_012, nd, model) {
  model_snps <- as.character(nd$snp_names)
  available <- intersect(model_snps, colnames(geno_012))

  sample_mat <- matrix(NA_real_, nrow = nrow(geno_012), ncol = length(model_snps))
  colnames(sample_mat) <- model_snps
  rownames(sample_mat) <- geno_012$sample_id

  # Fill available SNPs
  if (length(available) > 0) {
    sample_mat[, available] <- dt_cols_to_matrix(geno_012, available)
  }

  # Mean-impute using reference means
  ref_means <- as.numeric(nd$ref_means)
  for (j in seq_along(model_snps)) {
    na_mask <- is.na(sample_mat[, j])
    if (any(na_mask)) {
      sample_mat[na_mask, j] <- round(ref_means[j])
    }
  }

  pred <- predict(model, newdata = sample_mat)

  list(
    assignment = as.character(pred$assign),
    max_posterior = apply(pred$posterior, 1, max),
    posteriors = pred$posterior
  )
}

#' Predict at a single SVM node
#' @param geno_012 data.table with sample_id + SNP columns
#' @param nd Node data (must have svm_snp_names, svm_ref_means)
#' @param model SVM model object
#' @return list(assignment, max_probability, probabilities) or NULL if model unavailable
predict_svm <- function(geno_012, nd, model) {
  if (is.null(model) || is.null(nd$svm_snp_names)) return(NULL)

  svm_snps <- as.character(nd$svm_snp_names)
  available <- intersect(svm_snps, colnames(geno_012))

  sample_mat <- matrix(NA_real_, nrow = nrow(geno_012), ncol = length(svm_snps))
  colnames(sample_mat) <- svm_snps
  rownames(sample_mat) <- geno_012$sample_id

  if (length(available) > 0) {
    sample_mat[, available] <- dt_cols_to_matrix(geno_012, available)
  }

  # Mean-impute
  ref_means <- as.numeric(nd$svm_ref_means)
  if (length(ref_means) == 0 || all(is.na(ref_means))) ref_means <- as.numeric(nd$ref_means)
  for (j in seq_along(svm_snps)) {
    na_mask <- is.na(sample_mat[, j])
    if (any(na_mask)) {
      m <- if (j <= length(ref_means)) ref_means[j] else 1
      sample_mat[na_mask, j] <- round(m)
    }
  }

  svm_df <- data.frame(sample_mat, check.names = FALSE)

  p <- tryCatch(
    predict(model, svm_df, probability = TRUE),
    error = function(e) NULL
  )
  if (is.null(p)) return(NULL)

  probs <- attr(p, "probabilities")
  list(
    assignment = as.character(p),
    max_probability = apply(probs, 1, max),
    probabilities = probs
  )
}

#' Predict using popfinder NN weights (pure R forward pass)
#' @param geno_012 data.table with sample_id + SNP columns
#' @param nn_weights List with W1,b1,W2,b2,W3,b3,scaler_mean,scaler_scale,classes,snp_names
#' @return list(assignment, max_probability, probabilities) or NULL
predict_nn <- function(geno_012, nn_weights) {
  if (is.null(nn_weights)) return(NULL)

  nn_snps <- as.character(nn_weights$snp_names)
  available <- intersect(nn_snps, colnames(geno_012))

  sample_mat <- matrix(0, nrow = nrow(geno_012), ncol = length(nn_snps))
  colnames(sample_mat) <- nn_snps

  if (length(available) > 0) {
    sample_mat[, available] <- dt_cols_to_matrix(geno_012, available)
  }

  # Replace NA with column mean (from scaler_mean which is the training mean)
  scaler_mean <- as.numeric(unlist(nn_weights$scaler_mean))
  scaler_scale <- as.numeric(unlist(nn_weights$scaler_scale))
  for (j in seq_along(nn_snps)) {
    na_mask <- is.na(sample_mat[, j])
    if (any(na_mask)) {
      sample_mat[na_mask, j] <- scaler_mean[j]
    }
  }

  # StandardScaler transform: (x - mean) / scale
  x <- sweep(sample_mat, 2, scaler_mean, "-")
  x <- sweep(x, 2, scaler_scale, "/")
  # Handle zero-scale columns
  x[is.nan(x)] <- 0
  x[is.infinite(x)] <- 0

  # Forward pass through 3-layer network
  # JSON loaded with simplifyVector=FALSE stores matrices as list-of-lists;
  # unlist + matrix is the safest conversion to a numeric matrix.
  W1 <- matrix(as.numeric(unlist(nn_weights$W1)), nrow = length(nn_weights$W1), byrow = TRUE)
  b1 <- as.numeric(unlist(nn_weights$b1))
  W2 <- matrix(as.numeric(unlist(nn_weights$W2)), nrow = length(nn_weights$W2), byrow = TRUE)
  b2 <- as.numeric(unlist(nn_weights$b2))
  W3 <- matrix(as.numeric(unlist(nn_weights$W3)), nrow = length(nn_weights$W3), byrow = TRUE)
  b3 <- as.numeric(unlist(nn_weights$b3))

  # Layer 1: Linear + ReLU
  h1 <- sweep(x %*% t(W1), 2, b1, "+")
  h1 <- pmax(h1, 0)

  # Layer 2: Linear + ReLU
  h2 <- sweep(h1 %*% t(W2), 2, b2, "+")
  h2 <- pmax(h2, 0)

  # Layer 3: Linear (logits)
  logits <- sweep(h2 %*% t(W3), 2, b3, "+")

  # Softmax with numerical stability
  logits_stable <- logits - apply(logits, 1, max)
  exp_logits <- exp(logits_stable)
  probs <- exp_logits / rowSums(exp_logits)

  classes <- as.character(unlist(nn_weights$classes))
  colnames(probs) <- classes

  pred_idx <- apply(probs, 1, which.max)
  list(
    assignment = classes[pred_idx],
    max_probability = apply(probs, 1, max),
    probabilities = probs
  )
}

#' Apply SVM-primary ensemble decision
#' @param svm_result SVM prediction result (or NULL)
#' @param dapc_result DAPC prediction result
#' @param pf_result Popfinder prediction result (or NULL)
#' @param svm_thresh SVM calibrated thresholds
#' @return data.table with assignment, tier, method details
ensemble_decision <- function(svm_result, dapc_result, pf_result, svm_thresh) {
  n <- length(dapc_result$assignment)

  result <- data.table(
    svm_group = if (!is.null(svm_result)) svm_result$assignment else rep(NA_character_, n),
    svm_probability = if (!is.null(svm_result)) svm_result$max_probability else rep(NA_real_, n),
    dapc_group = dapc_result$assignment,
    dapc_posterior = dapc_result$max_posterior,
    nn_group = if (!is.null(pf_result)) pf_result$assignment else rep(NA_character_, n),
    nn_probability = if (!is.null(pf_result)) pf_result$max_probability else rep(NA_real_, n),
    assignment = NA_character_,
    tier = NA_character_,
    ensemble_prob = NA_real_
  )

  for (i in seq_len(n)) {
    svm_g <- result$svm_group[i]
    dapc_g <- result$dapc_group[i]
    pf_g <- result$nn_group[i]
    svm_p <- result$svm_probability[i]

    # Primary assignment: SVM if available, otherwise DAPC
    if (!is.na(svm_g)) {
      result$assignment[i] <- svm_g
      dapc_agrees <- !is.na(dapc_g) && dapc_g == svm_g
      pf_agrees <- !is.na(pf_g) && pf_g == svm_g
      pf_ran <- !is.na(pf_g)

      if (pf_ran) {
        if (dapc_agrees && pf_agrees) {
          result$tier[i] <- "Strongly Supported"
        } else if (dapc_agrees || pf_agrees) {
          result$tier[i] <- "Supported"
        } else {
          result$tier[i] <- "Provisional"
        }
      } else {
        if (dapc_agrees) {
          result$tier[i] <- "Strongly Supported"
        } else if (!is.na(svm_p) && svm_p >= svm_thresh$upper) {
          result$tier[i] <- "Supported"
        } else {
          result$tier[i] <- "Provisional"
        }
      }
    } else {
      # No SVM: fall back to DAPC only
      result$assignment[i] <- dapc_g
      result$tier[i] <- "Provisional"
    }
  }

  # Ensemble probability = mean of each method's probability for the ASSIGNED group
  # (not max probability per method, which inflates confidence when methods disagree)
  for (i in seq_len(n)) {
    asgn <- result$assignment[i]
    if (is.na(asgn)) next

    vals <- numeric(0)
    if (!is.null(svm_result) && asgn %in% colnames(svm_result$probabilities)) {
      vals <- c(vals, svm_result$probabilities[i, asgn])
    }
    if (asgn %in% colnames(dapc_result$posteriors)) {
      vals <- c(vals, dapc_result$posteriors[i, asgn])
    }
    if (!is.null(pf_result) && asgn %in% colnames(pf_result$probabilities)) {
      vals <- c(vals, pf_result$probabilities[i, asgn])
    }
    result$ensemble_prob[i] <- if (length(vals) > 0) mean(vals) else NA_real_
  }

  result
}

#' Run full hierarchical classification
#' @param geno_012 data.table with sample_id + SNP columns (012-coded)
#' @param node_data Named list of node metadata
#' @param dapc_models Named list of DAPC models
#' @param svm_models Named list of SVM models
#' @param nn_weights Named list of popfinder NN weight objects
#' @param hierarchy Parsed hierarchy definition
#' @param routing Routing table from build_routing_table()
#' @param progress Optional shiny progress callback
#' @return list(results = data.table, posterior_details = list)
classify_samples <- function(geno_012, node_data, dapc_models, svm_models,
                             nn_weights, hierarchy, routing,
                             ecotype_display_names = list(), progress = NULL) {
  samples <- geno_012$sample_id
  n_samples <- length(samples)

  results <- data.table(
    sample_id = samples,
    level1_assignment = NA_character_, level1_tier = NA_character_,
    level1_svm_group = NA_character_, level1_svm_probability = NA_real_,
    level1_dapc_group = NA_character_, level1_dapc_posterior = NA_real_,
    level1_nn_group = NA_character_, level1_nn_probability = NA_real_,
    level1_ensemble_prob = NA_real_,
    level2_assignment = NA_character_, level2_tier = NA_character_,
    level2_svm_group = NA_character_, level2_svm_probability = NA_real_,
    level2_dapc_group = NA_character_, level2_dapc_posterior = NA_real_,
    level2_nn_group = NA_character_, level2_nn_probability = NA_real_,
    level2_ensemble_prob = NA_real_,
    level3_assignment = NA_character_, level3_tier = NA_character_,
    level3_svm_group = NA_character_, level3_svm_probability = NA_real_,
    level3_dapc_group = NA_character_, level3_dapc_posterior = NA_real_,
    level3_nn_group = NA_character_, level3_nn_probability = NA_real_,
    level3_ensemble_prob = NA_real_,
    deepest_confident_assignment = NA_character_
  )

  posterior_details <- list()

  # ---- Level 1: Root ----
  if (!is.null(progress)) progress$set(message = "Classifying Level 1 (Subspecies)...", value = 0.1)

  root_nd <- node_data[["root"]]
  root_dapc <- predict_dapc("root", geno_012, root_nd, dapc_models[["root"]])
  root_svm <- predict_svm(geno_012, root_nd, svm_models[["root"]])
  root_pf <- predict_nn(geno_012, nn_weights[["root"]])

  root_ensemble <- ensemble_decision(root_svm, root_dapc, root_pf, root_nd$svm_thresholds)

  results$level1_assignment <- root_ensemble$assignment
  results$level1_tier <- root_ensemble$tier
  results$level1_svm_group <- root_ensemble$svm_group
  results$level1_svm_probability <- root_ensemble$svm_probability
  results$level1_dapc_group <- root_ensemble$dapc_group
  results$level1_dapc_posterior <- root_ensemble$dapc_posterior
  results$level1_nn_group <- root_ensemble$nn_group
  results$level1_nn_probability <- root_ensemble$nn_probability
  results$level1_ensemble_prob <- root_ensemble$ensemble_prob

  posterior_details[["root"]] <- list(
    dapc = root_dapc$posteriors,
    svm = if (!is.null(root_svm)) root_svm$probabilities else NULL,
    pf = if (!is.null(root_pf)) root_pf$probabilities else NULL
  )

  # ---- Level 2 ----
  if (!is.null(progress)) progress$set(message = "Classifying Level 2 (Ecotype)...", value = 0.4)

  confident_l1 <- c("Strongly Supported", "Supported")
  unique_l1 <- unique(results$level1_assignment[!is.na(results$level1_assignment)])

  for (l1_group in unique_l1) {
    l2_node_name <- get_l2_node(l1_group, routing)
    if (is.na(l2_node_name) || is.null(l2_node_name)) next

    l2_nd <- node_data[[l2_node_name]]
    if (is.null(l2_nd)) next

    mask <- which(results$level1_assignment == l1_group &
                  results$level1_tier %in% confident_l1)
    if (length(mask) == 0) next

    if (l2_nd$n_groups < 2) {
      # Single-group terminal node
      results$level2_assignment[mask] <- l2_nd$groups[1]
      results$level2_tier[mask] <- "Strongly Supported"
      results$level2_dapc_group[mask] <- l2_nd$groups[1]
      results$level2_dapc_posterior[mask] <- 1.0
      results$level2_svm_group[mask] <- l2_nd$groups[1]
      results$level2_svm_probability[mask] <- 1.0
      results$level2_nn_group[mask] <- l2_nd$groups[1]
      results$level2_nn_probability[mask] <- 1.0
      results$level2_ensemble_prob[mask] <- 1.0
    } else {
      subset_geno <- geno_012[mask, ]
      l2_dapc <- predict_dapc(l2_node_name, subset_geno, l2_nd, dapc_models[[l2_node_name]])
      l2_svm <- predict_svm(subset_geno, l2_nd, svm_models[[l2_node_name]])
      l2_pf <- predict_nn(subset_geno, nn_weights[[l2_node_name]])
      l2_ensemble <- ensemble_decision(l2_svm, l2_dapc, l2_pf, l2_nd$svm_thresholds)

      results$level2_assignment[mask] <- l2_ensemble$assignment
      results$level2_tier[mask] <- l2_ensemble$tier
      results$level2_svm_group[mask] <- l2_ensemble$svm_group
      results$level2_svm_probability[mask] <- l2_ensemble$svm_probability
      results$level2_dapc_group[mask] <- l2_ensemble$dapc_group
      results$level2_dapc_posterior[mask] <- l2_ensemble$dapc_posterior
      results$level2_nn_group[mask] <- l2_ensemble$nn_group
      results$level2_nn_probability[mask] <- l2_ensemble$nn_probability
      results$level2_ensemble_prob[mask] <- l2_ensemble$ensemble_prob

      posterior_details[[l2_node_name]] <- list(
        dapc = l2_dapc$posteriors,
        svm = if (!is.null(l2_svm)) l2_svm$probabilities else NULL,
        pf = if (!is.null(l2_pf)) l2_pf$probabilities else NULL,
        sample_ids = subset_geno$sample_id
      )
    }
  }

  # ---- Level 3 ----
  if (!is.null(progress)) progress$set(message = "Classifying Level 3 (Herd)...", value = 0.7)

  confident_l2 <- c("Strongly Supported", "Supported")
  l2_assigned <- which(!is.na(results$level2_assignment) &
                       results$level2_tier %in% confident_l2)

  for (idx in l2_assigned) {
    l1_asgn <- results$level1_assignment[idx]
    l2_asgn <- results$level2_assignment[idx]
    l3_node_name <- get_l3_node(l1_asgn, l2_asgn, routing)
    if (is.na(l3_node_name) || is.null(l3_node_name)) next

    l3_nd <- node_data[[l3_node_name]]
    if (is.null(l3_nd)) next

    if (l3_nd$n_groups < 2) {
      results$level3_assignment[idx] <- l3_nd$groups[1]
      results$level3_tier[idx] <- "Strongly Supported"
      results$level3_dapc_group[idx] <- l3_nd$groups[1]
      results$level3_dapc_posterior[idx] <- 1.0
      results$level3_svm_group[idx] <- l3_nd$groups[1]
      results$level3_svm_probability[idx] <- 1.0
      results$level3_nn_group[idx] <- l3_nd$groups[1]
      results$level3_nn_probability[idx] <- 1.0
      results$level3_ensemble_prob[idx] <- 1.0
    }
  }

  # Batch L3 multi-group predictions by node
  l3_multigroup_nodes <- unique(na.omit(sapply(l2_assigned, function(idx) {
    l3_name <- get_l3_node(results$level1_assignment[idx],
                           results$level2_assignment[idx], routing)
    if (!is.null(l3_name) && !is.null(node_data[[l3_name]]) &&
        node_data[[l3_name]]$n_groups >= 2) l3_name else NA_character_
  })))

  for (l3_node_name in l3_multigroup_nodes) {
    l3_nd <- node_data[[l3_node_name]]

    # Find which samples route to this L3 node
    mask <- which(sapply(seq_len(nrow(results)), function(idx) {
      if (is.na(results$level2_assignment[idx]) ||
          !(results$level2_tier[idx] %in% confident_l2)) return(FALSE)
      l3_name <- get_l3_node(results$level1_assignment[idx],
                             results$level2_assignment[idx], routing)
      !is.na(l3_name) && l3_name == l3_node_name
    }))

    if (length(mask) == 0) next

    subset_geno <- geno_012[mask, ]
    l3_dapc <- predict_dapc(l3_node_name, subset_geno, l3_nd, dapc_models[[l3_node_name]])
    l3_svm <- predict_svm(subset_geno, l3_nd, svm_models[[l3_node_name]])
    l3_pf <- predict_nn(subset_geno, nn_weights[[l3_node_name]])
    l3_ensemble <- ensemble_decision(l3_svm, l3_dapc, l3_pf, l3_nd$svm_thresholds)

    results$level3_assignment[mask] <- l3_ensemble$assignment
    results$level3_tier[mask] <- l3_ensemble$tier
    results$level3_svm_group[mask] <- l3_ensemble$svm_group
    results$level3_svm_probability[mask] <- l3_ensemble$svm_probability
    results$level3_dapc_group[mask] <- l3_ensemble$dapc_group
    results$level3_dapc_posterior[mask] <- l3_ensemble$dapc_posterior
    results$level3_nn_group[mask] <- l3_ensemble$nn_group
    results$level3_nn_probability[mask] <- l3_ensemble$nn_probability
    results$level3_ensemble_prob[mask] <- l3_ensemble$ensemble_prob

    posterior_details[[l3_node_name]] <- list(
      dapc = l3_dapc$posteriors,
      svm = if (!is.null(l3_svm)) l3_svm$probabilities else NULL,
      pf = if (!is.null(l3_pf)) l3_pf$probabilities else NULL,
      sample_ids = subset_geno$sample_id
    )
  }

  # ---- Replace "Not applicable" ecotypes with herd names ----
  # Monotypic subspecies (fennicus, pearyi) have "Not applicable" at L2;
  # replace with the actual herd name for cleaner display.
  if (length(ecotype_display_names) > 0) {
    for (i in seq_len(nrow(results))) {
      l1 <- results$level1_assignment[i]
      if (!is.na(l1) && l1 %in% names(ecotype_display_names)) {
        display_name <- ecotype_display_names[[l1]]
        if (!is.na(results$level2_assignment[i]) &&
            results$level2_assignment[i] == "Not applicable") {
          set(results, i, "level2_assignment", display_name)
        }
        if (!is.na(results$level2_dapc_group[i]) &&
            results$level2_dapc_group[i] == "Not applicable") {
          set(results, i, "level2_dapc_group", display_name)
        }
        if (!is.na(results$level2_svm_group[i]) &&
            results$level2_svm_group[i] == "Not applicable") {
          set(results, i, "level2_svm_group", display_name)
        }
        if (!is.na(results$level2_nn_group[i]) &&
            results$level2_nn_group[i] == "Not applicable") {
          set(results, i, "level2_nn_group", display_name)
        }
      }
    }
  }

  # ---- Deepest confident assignment ----
  if (!is.null(progress)) progress$set(message = "Compiling results...", value = 0.9)

  confident_tiers <- c("Strongly Supported", "Supported")
  results[!is.na(level1_tier) & level1_tier %in% confident_tiers,
          deepest_confident_assignment := level1_assignment]
  results[!is.na(level2_tier) & level2_tier %in% confident_tiers,
          deepest_confident_assignment := level2_assignment]
  results[!is.na(level3_tier) & level3_tier %in% confident_tiers,
          deepest_confident_assignment := level3_assignment]

  list(results = results, posterior_details = posterior_details)
}
