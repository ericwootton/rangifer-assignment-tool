# Rangifer Assignment Tool

A Shiny application for probabilistic assignment of caribou and reindeer
(*Rangifer tarandus*) samples to subspecies, ecotype, and herd using SNP
genotype data from the RangiferSNP63K chip.

## Overview

The tool classifies samples through a three-level hierarchy (subspecies,
ecotype, herd/population) using an ensemble of three statistical methods:

- **SVM** (Support Vector Machine) -- primary classifier, operating on the top
  5,000 F_ST-ranked SNPs per node
- **DAPC** (Discriminant Analysis of Principal Components) -- validator, using
  all available SNPs per node
- **Neural network** (popfinder architecture) -- validator, applied at nodes
  with sufficient group sizes

Confidence in each assignment is expressed as **Strongly Supported**,
**Supported**, or **Provisional**, based on agreement across methods.

The reference panel contains 362 quality-controlled samples spanning five
recognized subspecies, assembled from multiple chip versions and whole-genome
sequencing datasets.

## Requirements

- **R** >= 4.1.0 ([download](https://cran.r-project.org/))

The following R packages are required and will be installed automatically on
first launch if missing:

`shiny`, `data.table`, `adegenet`, `e1071`, `jsonlite`, `DT`, `plotly`,
`bslib`, `shinyjs`

## Quick start

### Windows

Double-click **launch.bat**, or from a terminal:

```
Rscript launch.R
```

### macOS

Double-click **launch.command**, or from a terminal:

```
Rscript launch.R
```

You may need to make the launcher executable first:

```
chmod +x launch.command
```

### From RStudio

Open the project folder in RStudio and run:

```r
shiny::runApp()
```

The app will open in your default web browser.

## Preparing input data

This tool accepts an Illumina FinalReport file generated from the
RangiferSNP63K chip.

To create the file from your sequencing data:

1. Open your project (`.bsc` file) in Illumina GenomeStudio using the
   **Genotyping Module**.
2. Go to **Analysis > Reports > Report Wizard**.
3. Select **Final Report** from the list of report types.
4. Click through using the default options to generate the report.

The resulting file should contain a `[Data]` section with `SNP Name`,
`Sample ID`, and `Allele1/2 - Top` columns.

## Output

The app provides:

- **Upload & QC** -- per-sample call rate, heterozygosity, quality flags, and
  genetic sex determination (SRY + X-chromosome heterozygosity)
- **All Results** -- assignment table with subspecies, ecotype, and herd calls,
  confidence tiers, and ensemble probabilities; downloadable as CSV
- **Individual Results** -- per-sample probability breakdowns and grouped bar
  charts for each classification level
- **Hierarchy** -- interactive sunburst plot of the classification tree with
  reference panel sizes and assignment counts

## Project structure

```
rangifer-assignment-tool/
├── app.R                  # Shiny UI and server
├── global.R               # Package loading and data initialization
├── launch.R               # Cross-platform launcher with dependency checks
├── launch.bat             # Windows double-click launcher
├── launch.command         # macOS double-click launcher
├── R/
│   ├── mod_classify.R     # Hierarchical classification engine
│   ├── utils_genotype.R   # FinalReport parsing and 012 conversion
│   └── utils_hierarchy.R  # Hierarchy traversal and routing
├── www/
│   └── styles.css         # Custom styles
├── data/
│   ├── hierarchy_definition.json
│   ├── chip_snp_map.rds
│   ├── snp_info.rds
│   ├── x_chromosome_snps.rds
│   ├── models/            # Pre-trained DAPC, SVM, and NN models
│   └── node_data/         # Per-node reference means and thresholds
└── prep/
    └── prepare_app_data.R # Bundles pipeline outputs into app data
```

## Reference

Carrier, A., et al. (2022). Design and validation of a 63K genome-wide
SNP-genotyping platform for caribou/reindeer (*Rangifer tarandus*). *BMC
Genomics*, 23(1), 687.
https://doi.org/10.1186/s12864-022-08899-6

## License

This project is licensed under the GNU General Public License v3.0. See
[LICENSE](LICENSE) for details.

## Contact

Eric Wootton -- eric.wootton@mail.mcgill.ca
