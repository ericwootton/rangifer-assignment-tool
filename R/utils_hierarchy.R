# ============================================================================
# Hierarchy Traversal Utilities
# ============================================================================

#' Build routing table from hierarchy definition
#' Maps assignments at each level to the next level's node name
#' @param hierarchy Parsed hierarchy_definition.json
#' @return list with l1_to_l2 and l2_to_l3 named character vectors
build_routing_table <- function(hierarchy) {
  l1_to_l2 <- character()
  l2_to_l3 <- character()

  # Level 1 -> Level 2: map each subspecies to its L2 node

  for (l1_group in names(hierarchy$level2)) {
    node_info <- hierarchy$level2[[l1_group]]
    l1_to_l2[l1_group] <- node_info$node_name
  }

  # Level 2 -> Level 3: map each L2 ecotype to its L3 node
  # Key format: "{L2_node_name}::{ecotype_group}"
  for (l2_key in names(hierarchy$level3)) {
    node_info <- hierarchy$level3[[l2_key]]
    parent_l2_node <- paste0("L2_", gsub(" ", "_", node_info$grandparent_group))
    routing_key <- paste0(parent_l2_node, "::", node_info$parent_group)
    l2_to_l3[routing_key] <- node_info$node_name
  }

  list(l1_to_l2 = l1_to_l2, l2_to_l3 = l2_to_l3)
}

#' Get the L2 node name for a given L1 assignment
get_l2_node <- function(l1_assignment, routing) {
  routing$l1_to_l2[l1_assignment]
}

#' Get the L3 node name for a given L1 + L2 assignment
get_l3_node <- function(l1_assignment, l2_assignment, routing) {
  l2_node <- paste0("L2_", gsub(" ", "_", l1_assignment))
  key <- paste0(l2_node, "::", l2_assignment)
  routing$l2_to_l3[key]
}

#' Build a data.frame for collapsibleTree from hierarchy definition
#' @param hierarchy Parsed hierarchy_definition.json
#' @param results Optional assignment results to add counts
#' @param ecotype_display_names Named list mapping subspecies -> display name for "Not applicable" ecotypes
#' @return data.frame with pathString column for tree building
build_tree_data <- function(hierarchy, results = NULL, ecotype_display_names = list()) {
  rows <- list()

  l1_groups <- if (is.list(hierarchy$level1$groups)) {
    unlist(hierarchy$level1$groups)
  } else {
    hierarchy$level1$groups
  }
  l1_counts <- if (is.list(hierarchy$level1$sample_counts)) {
    unlist(hierarchy$level1$sample_counts)
  } else {
    hierarchy$level1$sample_counts
  }

  for (i in seq_along(l1_groups)) {
    l1 <- l1_groups[i]
    l1_ref_n <- l1_counts[i]

    # Get L2 info
    l2_info <- hierarchy$level2[[l1]]
    l2_groups <- if (is.list(l2_info$groups)) unlist(l2_info$groups) else l2_info$groups
    l2_counts <- if (is.list(l2_info$sample_counts)) unlist(l2_info$sample_counts) else l2_info$sample_counts

    for (j in seq_along(l2_groups)) {
      l2 <- l2_groups[j]
      l2_ref_n <- l2_counts[j]

      # Replace "Not applicable" with display name
      if (l2 == "Not applicable" && l1 %in% names(ecotype_display_names)) {
        l2 <- ecotype_display_names[[l1]]
      }

      # Get L3 info
      l3_key <- paste0(gsub(" ", "_", l1), "_", gsub(" ", "_", l2_groups[j]))
      l3_info <- hierarchy$level3[[l3_key]]

      if (!is.null(l3_info)) {
        l3_groups <- if (is.list(l3_info$groups)) unlist(l3_info$groups) else l3_info$groups
        l3_counts <- if (is.list(l3_info$sample_counts)) unlist(l3_info$sample_counts) else l3_info$sample_counts

        for (k in seq_along(l3_groups)) {
          rows[[length(rows) + 1]] <- data.frame(
            level1 = l1, level2 = l2, level3 = l3_groups[k],
            ref_n = l3_counts[k],
            stringsAsFactors = FALSE
          )
        }
      } else {
        rows[[length(rows) + 1]] <- data.frame(
          level1 = l1, level2 = l2, level3 = l2,
          ref_n = l2_ref_n,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  tree_df <- do.call(rbind, rows)

  # Add assignment counts if results provided
  if (!is.null(results) && nrow(results) > 0) {
    tree_df$assigned_n <- 0
    for (r in seq_len(nrow(tree_df))) {
      mask <- results$level1_assignment == tree_df$level1[r] &
              results$level2_assignment == tree_df$level2[r] &
              results$level3_assignment == tree_df$level3[r]
      mask[is.na(mask)] <- FALSE
      tree_df$assigned_n[r] <- sum(mask)
    }
  }

  tree_df
}
