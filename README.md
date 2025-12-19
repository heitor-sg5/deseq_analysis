# Differential Expression Analysis

## Overview
This repository contains an R-based, modular workflow for performing differential gene expression analysis using the DESeq2 package. The project is reproducible and extensible, with clear separation between data input, analysis logic, and execution. It supports RNA-seq count data provided in .xlsx format and performs standard differential expression analysis along with exploratory and summary statistics. This repository serves both as a practical analysis tool and as a structured reference for RNA-seq differential expression processes in R.

---

## Table of Contents:
- [Overview](#differential-expression-analysis)
- [Content](#content)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Input](#input)

---

## Content
This analysis includes:

- Differential expression using DESeq2
- Summary statistics, including:
    - Number of significant genes
    - Upregulated versus downregulated genes
    - Strongest fold-change genes
    - Highest and lowest expressed genes
- Explanatory analysis, including:
    - Principal component analysis (PCA)
    - MA and Volcano plot visualizations
- Input validation for count and metadata tables
- Reproducible environment management using `renv`

---

## Project Structure

```deseq_analysis/
├── main.R
├── scripts/
│   ├── parser.R
│   ├── analysis.R
│   └── plots.R  
├── data/
│   └── example_data.xlsx
├── output/
│   ├── figures/
│   └── tables/ 
├── renv.lock
├── .Rprofile
├── deseq_analysis.Rproj
└── README.md
```

---

## Installation

Prerequisites:

- R (version ≥ 4.0)
- RStudio (optional, but recommended)

1. Clone the repository:

```bash
git clone https://github.com/heitor-sg5/deseq_analysis.git
cd deseq_analysis
```

2. Restore the R package environment:

```bash
install.packages("renv")
renv::restore()
```

This will install the exact package versions used for this project, including `DESeq2`, `ggplot2`, and `openxlsx`.

---

## Usage
1. Place your input `.xlsx` file in the `data/` directory.
2. Open the project in RStudio (`deseq_analysis.Rproj`).
3. Run the analysis:

```bash
source("main.R")
```

Alternatively, from the command line:

```bash
Rscript. main.R
```

4. Results, summaries, and plots will be generated and saved to the `output/` directory.

---

## Input

The input Excel file must contain two sheets:

1. `counts` sheet:

- First column: `gene` (unique gene identifiers)
- Remaining columns: raw integer counts for each sample

Example:

```bash
gene    sample1  sample2  sample3
GENE1   120      98       105
GENE2   0        3        1
```

2. `metadata` sheet:

- First column: `sample` (must match count column names exactly)
- Second column: `condition` (experimental condition)

Example:

```bash
sample      condition
sample1     control
sample2     treated
sample3     treated
```
