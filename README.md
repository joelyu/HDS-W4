# Healthcare Data Science Module 4 Assignment

## Cytotoxic immune infiltration effects on breast cancer survival and immunotherapy response

Submission repository for the Module 4 poster assignment (MSt Healthcare Data Science). Applies Danaher et al. (2017) immune gene scoring across three breast cancer cohorts (METABRIC, TCGA-BRCA, I-SPY2) to investigate the prognostic and predictive value of cytotoxic immune infiltration via Cox regression and logistic regression for immunotherapy response. A live build of the rendered report is [here](https://joelyu.github.io/HDS-W4/HDS_04_YuChungYan_2604.html)., and the poster [here](https://joelyu.github.io/HDS-W4/poster/poster.pdf).

### Environment Setup

Requires [Quarto](https://quarto.org/docs/get-started/) and [Miniforge](https://github.com/conda-forge/miniforge) (mamba).

```bash
# Create environment
mamba env create -f environment.yml

# Activate
mamba activate cytotoxcore

# Install Bioconductor packages (bioconda dependency chains broken on osx-arm64)
Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install(c("cBioPortalData", "GEOquery"), ask = FALSE, update = FALSE)'
```

Bioconductor packages are also auto-installed via `if(!require)` in the qmd setup chunk.

### Rendering

```bash
mamba activate cytotoxcore
quarto render HDS_04_YuChungYan_2604.qmd
```

The submission qmd sources scripts in order. No need to render intermediate files.

### File Structure

```
├── HDS_04_YuChungYan_2604.qmd       # Submission document — sources scripts/
├── scripts/
│   ├── _immune_markers.R             # Danaher gene lists (single source of truth)
│   ├── 00-fetch-data.R               # Cache-first: cBioPortal + GEO, 3-way gene harmonisation
│   ├── 01-immune-scoring.R           # 14 immune cell type scores per cohort
│   ├── 02-explore.R                  # Exploratory analysis
│   ├── 03-pca-cytotoxic.R            # PCA — cytotoxic axis identification
│   ├── 04-confirmatory-cox.R         # Cox regression: progressive adjustment, diagnostics
│   └── 05-ispy2.R                    # I-SPY2 logistic regression (pCR ~ cytotoxic)
├── poster/                           # Poster and assets (SVGs, QR code)
├── data/
│   └── processed/                    # Cleaned datasets (committed, cache-first pipeline)
├── environment.yml
├── references.bib
├── nature.csl
├── LICENSE.md
└── .github/workflows/render-qmd.yml  # GitHub Pages deploy
```

### Data

Three cohorts: METABRIC and TCGA-BRCA from [cBioPortal](https://www.cbioportal.org/), and I-SPY2 (GSE194040) from [GEO](https://www.ncbi.nlm.nih.gov/geo/). Immune scoring uses the Danaher et al. (2017) 14-cell-type gene panel (57/60 genes harmonised across platforms).

The data pipeline is cache-first: committed files in `data/processed/` are used if present; APIs are only called when the cache is missing.

### License

Code and analysis: CC BY 4.0. Data: see `LICENSE.md` for per-dataset terms.
