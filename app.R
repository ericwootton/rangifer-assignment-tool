# ============================================================================
# Rangifer Assignment Shiny App — Main Entry Point
# ============================================================================

source("global.R")

# ============================================================================
# UI
# ============================================================================
ui <- page_navbar(
  title = "Rangifer Assignment Tool",
  theme = bs_theme(bootswatch = "flatly"),
  # ---- Tab 1: Upload & QC ----
  header = tagList(
    useShinyjs(),
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"))
  ),
  nav_panel("Upload & QC",
    layout_sidebar(
      sidebar = sidebar(
        width = 350,
        h4("Input Data"),
        fileInput("genotype_file", "Upload File:",
          accept = c(".csv", ".tsv", ".txt")
        ),
        actionButton("run_assignment", "Run Assignment",
          class = "btn-primary btn-lg", width = "100%"
        ),
        hr(),
        h5("Status"),
        verbatimTextOutput("status_text")
      ),

      card(
        card_header("Sample QC Summary"),
        DTOutput("qc_table")
      ),
      card(
        card_header("Data Processing Summary"),
        verbatimTextOutput("processing_summary")
      )
    )
  ),

  # ---- Tab 2: All Results ----
  nav_panel("All Results",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        downloadButton("download_results", "Download CSV",
                   style = "background-color:#2c3e50; border-color:#2c3e50; color:#fff;"),
        hr(),
        h5("Confidence Tiers"),
        tags$div(
          tags$div(style = "color: #27ae60; font-weight: bold;", "Strongly Supported"),
          tags$p(class = "text-muted", style = "font-size: 0.85em; margin-left: 10px;",
            "All available methods agree on the assignment."),
          tags$div(style = "color: #f39c12; font-weight: bold;", "Supported"),
          tags$p(class = "text-muted", style = "font-size: 0.85em; margin-left: 10px;",
            "SVM + at least one validator agree, or SVM probability",
            " exceeds the calibrated upper threshold."),
          tags$div(style = "color: #e74c3c; font-weight: bold;", "Provisional"),
          tags$p(class = "text-muted", style = "font-size: 0.85em; margin-left: 10px;",
            "No validator supports the SVM assignment. Treat with caution.")
        ),
        hr(),
        h5("Gating Rule"),
        tags$p(class = "text-muted", style = "font-size: 0.85em;",
          "Samples advance to deeper levels (L2, L3) only if the parent",
          " level is Strongly Supported or Supported."),
        hr(),
        h5("Summary"),
        verbatimTextOutput("tier_summary")
      ),
      card(
        card_header("Assignment Results"),
        DTOutput("results_table")
      )
    )
  ),

  # ---- Tab 3: Individual Results ----
  nav_panel("Individual Results",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 280,
        selectInput("detail_sample", "Select Sample:", choices = NULL),
        hr(),
        htmlOutput("sample_summary")
      ),
      card(
        fill = FALSE,
        card_header("Probability Report (%)"),
        card_body(
          tags$p(class = "text-muted", style = "font-size: 0.85em;",
            "Probabilities (0\u2013100%) from each classification method.",
            " Ensemble = mean across available methods."),
          DTOutput("probability_report")
        )
      ),
      card(
        fill = FALSE,
        card_header("Level 1 \u2014 Subspecies"),
        uiOutput("l1_plot_ui")
      ),
      card(
        fill = FALSE,
        card_header("Level 2 \u2014 Ecotype"),
        uiOutput("l2_plot_ui")
      ),
      card(
        fill = FALSE,
        card_header("Level 3 \u2014 Herd"),
        uiOutput("l3_plot_ui")
      )
    )
  ),

  # ---- Tab 4: Hierarchy ----
  nav_panel("Hierarchy",
    card(
      full_screen = TRUE,
      card_header("Classification Hierarchy"),
      card_body(
        class = "p-0",
        plotlyOutput("hierarchy_sunburst", height = "600px", width = "100%")
      )
    )
  ),

  # ---- Tab 5: About ----
  nav_panel("About",
    card(
      card_header("About This Tool"),
      card_body(
        h4("Rangifer Assignment Tool"),
        p("Probabilistic assignment of caribou/reindeer (", tags$em("Rangifer tarandus"), ")",
          " samples to subspecies, ecotype, and herd using SNP genotype data."),
        p("Source code:",
          tags$a(href = "https://github.com/ericwootton/rangifer-assignment-tool",
                 "github.com/ericwootton/rangifer-assignment-tool", target = "_blank"),
          tags$br(),
          "Contact:",
          tags$a(href = "mailto:eric.wootton@mail.mcgill.ca", "eric.wootton@mail.mcgill.ca")),
        br(),
        h5(tags$b("Preparing Your Input Data")),
        p("This tool accepts an Illumina FinalReport file generated from the",
          " RangiferSNP63K chip. To create this file from your sequencing data:"),
        tags$ol(
          tags$li("Open your project (.bsc file) in Illumina GenomeStudio, using the",
                  tags$b(" Genotyping Module")),
          tags$li("In the top menu, go to ", tags$b("Analysis > Reports > Report Wizard")),
          tags$li("Select ", tags$b("Final Report"), " from the list of report types"),
          tags$li("Click through using the default options to generate the report")),
        p("The resulting file should contain a [Data] section with SNP Name, Sample ID,",
          " and Allele1/2 - Top columns. Upload it directly on the ",
          tags$b("Upload & QC"), " tab."),
        br(),
        h5(tags$b("Methods")),
        p("The assignment pipeline is based on the RangiferSNP63K genotyping platform",
          " (Carrier et al., 2022), a caribou-specific chip targeting ~63,000 SNPs",
          " distributed across the caribou genome (~1 per 50 kb), including variants on",
          " the X chromosome and SRY gene for sex determination."),
        p("The reference panel of", nrow(node_data[["root"]]$ref_means),
          "quality-controlled samples was assembled by combining genotypes from",
          " multiple versions of the SNP chip with genotypes extracted from",
          " whole-genome sequencing (WGS) datasets, including publicly available data",
          " from the NCBI Sequence Read Archive and samples from the NWT Cumulative",
          " Impact Monitoring Program (NWT CIMP). After integration across platforms,",
          " batch correction removed SNPs with significant platform-driven effects,",
          " retaining 11,677 SNPs for downstream analysis."),
        p("Samples are classified through a 3-level hierarchy",
          " (subspecies, ecotype, herd/population). At each node, a node-specific",
          " quality filter requires a minimum call rate of 80% and removal of",
          " monomorphic SNPs. Assignment uses an ensemble of three methods:"),
        tags$ul(
          tags$li(tags$b("SVM"), " (Support Vector Machine) \u2014 Primary classifier.",
                  " Operates on the top 5,000 most informative SNPs per node, ranked by",
                  " pairwise F", tags$sub("ST"), ". Consistently achieved the highest",
                  " cross-validation accuracy across all nodes."),
          tags$li(tags$b("DAPC"), " (Discriminant Analysis of Principal Components) \u2014",
                  " Validator. Uses all available SNPs per node with the number of",
                  " retained PCs selected via cross-validation."),
          tags$li(tags$b("Neural Network"), " (popfinder architecture) \u2014",
                  " Validator. Applied only at nodes where each group has at least",
                  " 5 samples.")
        ),
        p("At nodes where the neural network cannot be trained (due to small group",
          " sizes), the ensemble falls back to ", tags$b("SVM + DAPC Only"),
          " mode. In this mode, confidence is determined by agreement between those",
          " two methods and, when they disagree, by calibrated SVM probability",
          " thresholds (see below)."),
        br(),
        h5(tags$b("Confidence Tiers & Thresholds")),
        p("Confidence tiers reflect agreement across the available methods:"),
        tags$table(class = "table table-sm table-bordered",
          tags$thead(tags$tr(
            tags$th("Tier"),
            tags$th("SVM + Both Validators"),
            tags$th("SVM + DAPC Only")
          )),
          tags$tbody(
            tags$tr(
              tags$td(tags$span(style = "color: #27ae60; font-weight: bold;",
                                "Strongly Supported")),
              tags$td("All 3 methods agree"),
              tags$td("SVM and DAPC agree")
            ),
            tags$tr(
              tags$td(tags$span(style = "color: #f39c12; font-weight: bold;",
                                "Supported")),
              tags$td("SVM + 1 validator agree"),
              tags$td("SVM prob \u2265 upper threshold")
            ),
            tags$tr(
              tags$td(tags$span(style = "color: #e74c3c; font-weight: bold;",
                                "Provisional")),
              tags$td("No validator agrees"),
              tags$td("SVM & DAPC disagree, SVM prob < threshold")
            )
          )
        ),
        p("Each node has calibrated SVM probability thresholds derived from",
          " leave-one-out cross-validation. When only SVM and DAPC are available",
          " and they disagree, the SVM upper threshold (95th percentile of",
          " correct-prediction probabilities) determines whether the result is",
          " elevated from Provisional to Supported."),
        p("Assignment at deeper levels requires ", tags$b("Strongly Supported"),
          " or ", tags$b("Supported"),
          " confidence at the parent level."),
        br(),
        h5(tags$b("Sex Determination")),
        p("Genetic sex is determined using two complementary methods:"),
        tags$ol(
          tags$li(tags$b("SRY markers"), " \u2014 7 Y-linked probes on the chip.",
                  " Males show valid calls (6\u20137/7), females show missing calls (0\u20131/7)."),
          tags$li(tags$b("X-chromosome heterozygosity ratio"),
                  " \u2014 Ratio of heterozygosity on 844 X-linked SNPs",
                  " vs. autosomal SNPs. Males (XY) are hemizygous on X (ratio 0.06\u20130.31);",
                  " females (XX) are diploid (ratio 0.85\u20131.06).")
        ),
        p("SRY provides the primary call with X-het as support."),
        br(),
        h5(tags$b("References")),
        p(class = "text-muted", style = "font-size: 0.85em;",
          "Carrier, A., et al. (2022). Design and validation of a 63K genome-wide",
          " SNP-genotyping platform for caribou/reindeer (Rangifer tarandus).",
          tags$em(" BMC Genomics"), ", 23(1), 687.")
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================
server <- function(input, output, session) {

  # Reactive values
  rv <- reactiveValues(
    geno_012 = NULL,
    qc = NULL,
    sex = NULL,
    results = NULL,
    posterior_details = NULL,
    stats = NULL,
    status = "Ready. Upload a file and click 'Run Assignment'."
  )

  output$status_text <- renderText(rv$status)

  # ---- Process uploaded file ----
  observeEvent(input$run_assignment, {
    req(input$genotype_file)

    shinyjs::disable("run_assignment")
    on.exit(shinyjs::enable("run_assignment"))

    withProgress(message = "Processing...", value = 0, {
      tryCatch({
        file_path <- input$genotype_file$datapath

        incProgress(0.1, detail = "Parsing FinalReport...")
        rv$status <- "Parsing FinalReport..."

        result <- parse_finalreport(file_path, chip_snp_map, ref_lookup, alt_lookup)
        rv$geno_012 <- result$geno_012
        rv$stats <- result$stats
        rv$sex <- result$sex

        rv$status <- paste0("Parsed ", result$stats$n_samples, " samples, ",
                            result$stats$n_snps_mapped, " SNPs mapped")

        # QC
        incProgress(0.1, detail = "Computing QC metrics...")
        rv$qc <- compute_sample_qc(rv$geno_012)

        # Sex determination: X-het ratio (works for both input formats)
        xhet_result <- determine_sex_from_xhet(rv$geno_012, x_snp_ids)

        # Merge SRY + X-het into QC
        # SRY takes priority (binary, unambiguous), X-het as backup
        sex_combined <- merge(rv$sex, xhet_result, by = "sample_id", all = TRUE)
        sex_combined[, genetic_sex := fifelse(
          !is.na(genetic_sex), genetic_sex,  # SRY result
          x_het_sex                          # X-het fallback
        )]
        sex_combined[, sex_confidence := fifelse(
          !is.na(sry_calls) & sry_calls >= 5 | (!is.na(sry_calls) & sry_calls <= 1),
          sex_confidence,                    # SRY confidence
          x_het_confidence                   # X-het confidence
        )]
        rv$qc <- merge(rv$qc, sex_combined, by = "sample_id", all.x = TRUE)

        # Check SNP overlap with root model
        root_snps <- node_data[["root"]]$snp_names
        input_snps <- colnames(rv$geno_012)[-1]
        overlap <- length(intersect(root_snps, input_snps))
        overlap_pct <- round(overlap / length(root_snps) * 100, 1)
        rv$status <- paste0(rv$status, " | SNP overlap with model: ",
                            overlap, "/", length(root_snps), " (", overlap_pct, "%)")

        if (overlap_pct < 50) {
          showNotification(
            paste0("Warning: Only ", overlap_pct, "% SNP overlap with the reference panel. ",
                   "Results may be unreliable."),
            type = "warning", duration = 10
          )
        }

        # Classification
        incProgress(0.2, detail = "Running hierarchical classification...")
        rv$status <- "Running classification..."

        progress <- Progress$new(session)
        on.exit(progress$close())

        cls <- classify_samples(
          rv$geno_012, node_data, dapc_models, svm_models,
          nn_weights, hierarchy, routing, ecotype_display_names, progress
        )

        rv$results <- merge(rv$qc, cls$results, by = "sample_id", all.y = TRUE)
        rv$posterior_details <- cls$posterior_details

        # Update sample selector
        updateSelectInput(session, "detail_sample", choices = rv$results$sample_id)

        rv$status <- paste0("Assignment complete: ", nrow(rv$results), " samples classified")
        incProgress(1, detail = "Done!")

      }, error = function(e) {
        rv$status <- paste0("Error: ", e$message)
        showNotification(paste("Error:", e$message), type = "error", duration = 15)
      })
    })
  })

  # ---- QC Table ----
  output$qc_table <- renderDT({
    req(rv$qc)
    datatable(rv$qc, options = list(paging = FALSE, info = FALSE, scrollX = TRUE),
              rownames = FALSE) %>%
      formatStyle("call_rate",
        backgroundColor = styleInterval(c(0.70, 0.80), c("#e74c3c", "#f39c12", "#27ae60"))) %>%
      formatStyle("het_rate",
        backgroundColor = styleInterval(c(0.15, 0.50), c("#f39c12", "#27ae60", "#e74c3c")))
  })

  # ---- Processing Summary ----
  output$processing_summary <- renderText({
    req(rv$stats)
    s <- rv$stats
    sex_line <- ""
    if (!is.null(rv$qc) && "genetic_sex" %in% colnames(rv$qc)) {
      n_m <- sum(rv$qc$genetic_sex == "Male", na.rm = TRUE)
      n_f <- sum(rv$qc$genetic_sex == "Female", na.rm = TRUE)
      n_a <- sum(rv$qc$genetic_sex == "Ambiguous", na.rm = TRUE)
      sex_line <- paste0("\n\nSex determination (SRY + X-het ratio):\n",
                         "  Males: ", n_m, "  Females: ", n_f,
                         if (n_a > 0) paste0("  Ambiguous: ", n_a) else "")
    }
    paste0(
      "Samples: ", s$n_samples, "\n",
      "SNPs mapped: ", s$n_snps_mapped, "\n",
      "Direct allele match: ", s$n_direct, "\n",
      "Complement match: ", s$n_complement, "\n",
      "Strand-ambiguous (skipped): ", s$n_ambiguous, "\n",
      "Failed to match: ", s$n_failed, "\n",
      "Unmapped SNPs: ", s$n_unmapped,
      sex_line
    )
  })

  # ---- Results Table ----
  output$results_table <- renderDT({
    req(rv$results)

    display_cols <- c("sample_id", "genetic_sex", "x_het_ratio", "call_rate",
                      "level1_assignment", "level1_tier",
                      "level1_ensemble_prob",
                      "level2_assignment", "level2_tier",
                      "level2_ensemble_prob",
                      "level3_assignment", "level3_tier",
                      "level3_ensemble_prob",
                      "deepest_confident_assignment", "quality_flags")
    display_cols <- intersect(display_cols, colnames(rv$results))

    dt <- as.data.frame(rv$results)[display_cols]

    # Format ensemble probabilities as percentage strings in R (avoids DT JS issues)
    for (col in c("level1_ensemble_prob", "level2_ensemble_prob", "level3_ensemble_prob")) {
      if (col %in% names(dt)) {
        dt[[col]] <- ifelse(is.na(dt[[col]]), "",
                            paste0(round(dt[[col]] * 100, 1), "%"))
      }
    }

    tier_vals <- c("Strongly Supported", "Supported", "Provisional")
    tier_cols <- c("#d5f5e3", "#fdebd0", "#fadbd8")

    tbl <- datatable(dt,
              options = list(paging = FALSE, info = FALSE, scrollX = TRUE),
              rownames = FALSE, selection = "single")

    for (prefix in c("level1", "level2", "level3")) {
      tier_col <- paste0(prefix, "_tier")
      if (tier_col %in% display_cols)
        tbl <- tbl %>% formatStyle(tier_col,
          backgroundColor = styleEqual(tier_vals, tier_cols))
    }
    tbl
  })

  output$tier_summary <- renderText({
    req(rv$results)
    r <- rv$results
    paste0(
      "Level 1:\n",
      "  Strongly Supported: ", sum(r$level1_tier == "Strongly Supported", na.rm = TRUE), "\n",
      "  Supported: ", sum(r$level1_tier == "Supported", na.rm = TRUE), "\n",
      "  Provisional: ", sum(r$level1_tier == "Provisional", na.rm = TRUE), "\n\n",
      "Level 2:\n",
      "  Strongly Supported: ", sum(r$level2_tier == "Strongly Supported", na.rm = TRUE), "\n",
      "  Supported: ", sum(r$level2_tier == "Supported", na.rm = TRUE), "\n",
      "  Provisional: ", sum(r$level2_tier == "Provisional", na.rm = TRUE), "\n\n",
      "Level 3:\n",
      "  Strongly Supported: ", sum(r$level3_tier == "Strongly Supported", na.rm = TRUE), "\n",
      "  Supported: ", sum(r$level3_tier == "Supported", na.rm = TRUE), "\n",
      "  Provisional: ", sum(r$level3_tier == "Provisional", na.rm = TRUE)
    )
  })

  output$download_results <- downloadHandler(
    filename = function() {
      paste0("rangifer_assignments_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      fwrite(rv$results, file)
    }
  )

  # ---- Hierarchy Sunburst ----
  output$hierarchy_sunburst <- renderPlotly({
    tree_df <- build_tree_data(hierarchy, rv$results, ecotype_display_names)

    # Build sunburst data
    ids <- character()
    labels <- character()
    parents <- character()
    values <- numeric()
    colors <- character()

    # Root
    ids <- c(ids, "Rangifer")
    labels <- c(labels, "Rangifer tarandus")
    parents <- c(parents, "")
    total_ref <- sum(tree_df$ref_n)
    values <- c(values, total_ref)
    colors <- c(colors, "#95a5a6")

    # L1
    l1_groups <- unique(tree_df$level1)
    for (l1 in l1_groups) {
      subset <- tree_df[tree_df$level1 == l1, ]
      n_ref <- sum(subset$ref_n)
      n_assigned <- if ("assigned_n" %in% colnames(subset)) sum(subset$assigned_n) else 0
      lbl <- if (n_assigned > 0) paste0(l1, " (", n_assigned, ")") else l1

      ids <- c(ids, l1)
      labels <- c(labels, lbl)
      parents <- c(parents, "Rangifer")
      values <- c(values, n_ref)
      colors <- c(colors, "#3498db")
    }

    # L2
    l2_combos <- unique(tree_df[, c("level1", "level2")])
    for (r in seq_len(nrow(l2_combos))) {
      l1 <- l2_combos$level1[r]
      l2 <- l2_combos$level2[r]
      subset <- tree_df[tree_df$level1 == l1 & tree_df$level2 == l2, ]
      n_ref <- sum(subset$ref_n)
      n_assigned <- if ("assigned_n" %in% colnames(subset)) sum(subset$assigned_n) else 0
      lbl <- if (n_assigned > 0) paste0(l2, " (", n_assigned, ")") else l2
      id <- paste0(l1, "/", l2)

      ids <- c(ids, id)
      labels <- c(labels, lbl)
      parents <- c(parents, l1)
      values <- c(values, n_ref)
      colors <- c(colors, "#2ecc71")
    }

    # L3
    for (r in seq_len(nrow(tree_df))) {
      l1 <- tree_df$level1[r]
      l2 <- tree_df$level2[r]
      l3 <- tree_df$level3[r]
      n_ref <- tree_df$ref_n[r]
      n_assigned <- if ("assigned_n" %in% colnames(tree_df)) tree_df$assigned_n[r] else 0
      lbl <- if (n_assigned > 0) paste0(l3, " (", n_assigned, ")") else l3
      parent_id <- paste0(l1, "/", l2)
      id <- paste0(l1, "/", l2, "/", l3)

      ids <- c(ids, id)
      labels <- c(labels, lbl)
      parents <- c(parents, parent_id)
      values <- c(values, n_ref)
      colors <- c(colors, "#e67e22")
    }

    plot_ly(
      ids = ids,
      labels = labels,
      parents = parents,
      values = values,
      type = "sunburst",
      branchvalues = "total",
      marker = list(colors = colors),
      textinfo = "label",
      hovertemplate = "%{label}<br>Ref samples: %{value}<extra></extra>"
    ) %>%
      layout(
        margin = list(t = 10, l = 10, r = 10, b = 10)
      ) %>%
      config(responsive = TRUE)
  })

  # ---- Sample Detail ----
  output$sample_summary <- renderUI({
    req(rv$results, input$detail_sample)
    r <- rv$results[sample_id == input$detail_sample]
    if (nrow(r) == 0) return(tags$p("Sample not found"))

    tier_color <- function(tier) {
      switch(as.character(tier),
        "Strongly Supported" = "#27ae60",
        "Supported" = "#f39c12",
        "Provisional" = "#e74c3c",
        "#999999"
      )
    }

    make_level_block <- function(label, assignment, tier, ensemble_prob) {
      if (is.na(assignment)) {
        return(tags$div(style = "margin-bottom: 12px;",
          tags$div(style = "font-weight: 600; color: #999;", label),
          tags$div(style = "color: #999;", "\u2014")
        ))
      }
      pct <- if (!is.na(ensemble_prob)) paste0(round(ensemble_prob * 100, 1), "%") else "\u2014"
      tags$div(style = "margin-bottom: 12px;",
        tags$div(style = "font-weight: 600;", label),
        tags$div(style = "font-size: 1.05em;", assignment),
        tags$div(
          tags$span(style = paste0("color:", tier_color(tier), "; font-weight:bold;"), tier),
          tags$span(style = "color: #666; margin-left: 6px;", pct)
        )
      )
    }

    tagList(
      make_level_block("Subspecies",
        r$level1_assignment, r$level1_tier, r$level1_ensemble_prob),
      make_level_block("Ecotype",
        r$level2_assignment, r$level2_tier, r$level2_ensemble_prob),
      make_level_block("Herd",
        r$level3_assignment, r$level3_tier, r$level3_ensemble_prob)
    )
  })

  # Probability report table showing per-group percentages from all methods
  output$probability_report <- renderDT({
    req(rv$posterior_details, rv$results, input$detail_sample)
    sid <- input$detail_sample
    r <- rv$results[sample_id == sid]
    if (nrow(r) == 0) return(NULL)

    build_level_rows <- function(level_label, node_name) {
      pd <- rv$posterior_details[[node_name]]
      if (is.null(pd)) return(NULL)
      sample_ids <- if (!is.null(pd$sample_ids)) pd$sample_ids else rv$geno_012$sample_id
      row_idx <- which(sample_ids == sid)
      if (length(row_idx) == 0) return(NULL)

      # Use DAPC groups as canonical set (always available)
      if (is.null(pd$dapc)) return(NULL)
      groups <- colnames(pd$dapc)
      dapc_vals <- round(as.numeric(pd$dapc[row_idx, ]) * 100, 1)

      svm_vals <- rep(NA_real_, length(groups))
      if (!is.null(pd$svm)) {
        for (i in seq_along(groups)) {
          if (groups[i] %in% colnames(pd$svm))
            svm_vals[i] <- round(as.numeric(pd$svm[row_idx, groups[i]]) * 100, 1)
        }
      }

      nn_vals <- rep(NA_real_, length(groups))
      if (!is.null(pd$pf)) {
        for (i in seq_along(groups)) {
          if (groups[i] %in% colnames(pd$pf))
            nn_vals[i] <- round(as.numeric(pd$pf[row_idx, groups[i]]) * 100, 1)
        }
      }

      ensemble_vals <- numeric(length(groups))
      for (i in seq_along(groups)) {
        vals <- c(svm_vals[i], dapc_vals[i], nn_vals[i])
        vals <- vals[!is.na(vals)]
        ensemble_vals[i] <- if (length(vals) > 0) round(mean(vals), 1) else NA_real_
      }

      data.frame(
        Level = level_label,
        Group = groups,
        `SVM (%)` = svm_vals,
        `DAPC (%)` = dapc_vals,
        `NN (%)` = nn_vals,
        `Ensemble (%)` = ensemble_vals,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }

    rows <- list()
    l1_rows <- build_level_rows("L1: Subspecies", "root")
    if (!is.null(l1_rows)) rows <- c(rows, list(l1_rows))

    if (!is.na(r$level1_assignment)) {
      l2_node <- get_l2_node(r$level1_assignment, routing)
      if (!is.null(l2_node) && !is.na(l2_node)) {
        l2_rows <- build_level_rows("L2: Ecotype", l2_node)
        if (!is.null(l2_rows)) rows <- c(rows, list(l2_rows))
      }
    }

    if (!is.na(r$level1_assignment) && !is.na(r$level2_assignment)) {
      l3_node <- get_l3_node(r$level1_assignment, r$level2_assignment, routing)
      if (!is.null(l3_node) && !is.na(l3_node)) {
        l3_rows <- build_level_rows("L3: Herd", l3_node)
        if (!is.null(l3_rows)) rows <- c(rows, list(l3_rows))
      }
    }

    if (length(rows) == 0) return(NULL)
    report <- do.call(rbind, rows)

    datatable(report,
      options = list(paging = FALSE, searching = FALSE, info = FALSE,
                     ordering = FALSE, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatStyle("Ensemble (%)", fontWeight = "bold")
  })

  # Helper to make a grouped bar chart showing all 3 methods for a given level
  make_posterior_plot <- function(sid, level_prefix) {
    req(rv$posterior_details)

    r <- rv$results[sample_id == sid]

    node_name <- if (level_prefix == "root") {
      "root"
    } else if (grepl("L2", level_prefix)) {
      l1_asgn <- r$level1_assignment
      if (is.na(l1_asgn)) return(plotly_empty())
      get_l2_node(l1_asgn, routing)
    } else {
      l1_asgn <- r$level1_assignment
      l2_asgn <- r$level2_assignment
      if (is.na(l1_asgn) || is.na(l2_asgn)) return(plotly_empty())
      get_l3_node(l1_asgn, l2_asgn, routing)
    }

    if (is.null(node_name) || is.na(node_name)) {
      # No node found — show 100% bar for the assignment if one exists
      level_col <- switch(level_prefix,
        "root" = "level1_assignment", "L2" = "level2_assignment", "L3" = "level3_assignment")
      asgn <- r[[level_col]]
      if (is.na(asgn)) return(plotly_empty())
      groups <- asgn
      svm_probs <- dapc_probs <- nn_probs <- 1
    } else {
      pd <- rv$posterior_details[[node_name]]
      if (is.null(pd)) {
        # Single-group node — no posteriors stored, show 100% for the assignment
        level_col <- switch(level_prefix,
          "root" = "level1_assignment", "L2" = "level2_assignment", "L3" = "level3_assignment")
        asgn <- r[[level_col]]
        if (is.na(asgn)) return(plotly_empty())
        groups <- asgn
        svm_probs <- dapc_probs <- nn_probs <- 1
      } else {
        sample_ids <- if (!is.null(pd$sample_ids)) pd$sample_ids else rv$geno_012$sample_id
        row_idx <- which(sample_ids == sid)
        if (length(row_idx) == 0) return(plotly_empty())

        # Collect probabilities from all methods
        groups <- NULL
        dapc_probs <- svm_probs <- nn_probs <- NULL

        if (!is.null(pd$dapc)) {
          groups <- colnames(pd$dapc)
          dapc_probs <- as.numeric(pd$dapc[row_idx, ])
        }
        if (!is.null(pd$svm)) {
          if (is.null(groups)) groups <- colnames(pd$svm)
          svm_probs <- rep(0, length(groups))
          for (g in intersect(groups, colnames(pd$svm))) {
            svm_probs[which(groups == g)] <- as.numeric(pd$svm[row_idx, g])
          }
        }
        if (!is.null(pd$pf)) {
          if (is.null(groups)) groups <- colnames(pd$pf)
          nn_probs <- rep(0, length(groups))
          for (g in intersect(groups, colnames(pd$pf))) {
            nn_probs[which(groups == g)] <- as.numeric(pd$pf[row_idx, g])
          }
        }

        if (is.null(groups)) return(plotly_empty())
      }
    }

    p <- plot_ly()

    if (!is.null(nn_probs)) {
      p <- p %>% add_trace(
        y = groups, x = round(nn_probs * 100, 1), name = "Neural Net",
        type = "bar", orientation = "h",
        marker = list(color = "#16a085"),
        text = paste0(round(nn_probs * 100, 1), "%"), textposition = "outside",
        hovertemplate = "%{y}: %{x:.1f}%<extra>Neural Net</extra>"
      )
    }

    if (!is.null(dapc_probs)) {
      p <- p %>% add_trace(
        y = groups, x = round(dapc_probs * 100, 1), name = "DAPC",
        type = "bar", orientation = "h",
        marker = list(color = "#2980b9"),
        text = paste0(round(dapc_probs * 100, 1), "%"), textposition = "outside",
        hovertemplate = "%{y}: %{x:.1f}%<extra>DAPC</extra>"
      )
    }

    if (!is.null(svm_probs)) {
      p <- p %>% add_trace(
        y = groups, x = round(svm_probs * 100, 1), name = "SVM",
        type = "bar", orientation = "h",
        marker = list(color = "#8e44ad"),
        text = paste0(round(svm_probs * 100, 1), "%"), textposition = "outside",
        hovertemplate = "%{y}: %{x:.1f}%<extra>SVM</extra>"
      )
    }

    p %>% layout(
      barmode = "group",
      xaxis = list(title = "Probability (%)", range = c(0, 115)),
      yaxis = list(title = "", categoryorder = "total ascending",
                   ticksuffix = "  "),
      margin = list(l = 220),
      legend = list(orientation = "h", y = -0.2, traceorder = "reversed")
    )
  }

  # Compute number of groups for a given level to set dynamic plot height
  get_n_groups <- function(sid, level_prefix) {
    req(rv$results, rv$posterior_details)
    r <- rv$results[sample_id == sid]
    if (nrow(r) == 0) return(0)
    node_name <- if (level_prefix == "root") {
      "root"
    } else if (grepl("L2", level_prefix)) {
      l1 <- r$level1_assignment
      if (is.na(l1)) return(0)
      get_l2_node(l1, routing)
    } else {
      l1 <- r$level1_assignment; l2 <- r$level2_assignment
      if (is.na(l1) || is.na(l2)) return(0)
      get_l3_node(l1, l2, routing)
    }
    if (is.null(node_name) || is.na(node_name)) {
      # Single-group node with no posteriors — check if assignment exists
      level_col <- switch(level_prefix,
        "root" = "level1_assignment", "L2" = "level2_assignment", "L3" = "level3_assignment")
      return(if (!is.na(r[[level_col]])) 1 else 0)
    }
    pd <- rv$posterior_details[[node_name]]
    if (is.null(pd) || is.null(pd$dapc)) {
      level_col <- switch(level_prefix,
        "root" = "level1_assignment", "L2" = "level2_assignment", "L3" = "level3_assignment")
      return(if (!is.na(r[[level_col]])) 1 else 0)
    }
    ncol(pd$dapc)
  }

  plot_height <- function(n_groups) paste0(max(150, n_groups * 80 + 70), "px")

  output$l1_plot_ui <- renderUI({
    req(input$detail_sample, rv$posterior_details)
    h <- plot_height(get_n_groups(input$detail_sample, "root"))
    plotlyOutput("l1_posteriors", height = h)
  })
  output$l2_plot_ui <- renderUI({
    req(input$detail_sample, rv$posterior_details)
    h <- plot_height(get_n_groups(input$detail_sample, "L2"))
    plotlyOutput("l2_posteriors", height = h)
  })
  output$l3_plot_ui <- renderUI({
    req(input$detail_sample, rv$posterior_details)
    h <- plot_height(get_n_groups(input$detail_sample, "L3"))
    plotlyOutput("l3_posteriors", height = h)
  })

  output$l1_posteriors <- renderPlotly({
    req(input$detail_sample)
    make_posterior_plot(sid = input$detail_sample, "root")
  })

  output$l2_posteriors <- renderPlotly({
    req(input$detail_sample)
    make_posterior_plot(sid = input$detail_sample, "L2")
  })

  output$l3_posteriors <- renderPlotly({
    req(input$detail_sample)
    make_posterior_plot(sid = input$detail_sample, "L3")
  })
}

shinyApp(ui, server)
