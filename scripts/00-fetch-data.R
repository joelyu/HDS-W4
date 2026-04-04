# =============================================================================
# 00-fetch-data.R
# Fetch METABRIC + TCGA-BRCA (cBioPortal) + I-SPY2 (GEO) → data/processed/ and data/raw/
#
# Two cohorts:
#   1. METABRIC (brca_metabric) — Illumina HT-12 v3, log2 intensity
#   2. TCGA-BRCA (brca_tcga_pan_can_atlas_2018) — RNA-Seq V2, RSEM
#
# Output per cohort:
#   data/processed/{cohort}_clinical.csv
#   data/processed/{cohort}_expression.csv   (gene × patient, log2 scale)
#
# Expects: proc_dir defined by parent qmd (falls back to data/processed/)
# =============================================================================

if (!exists("proc_dir")) {
  proc_dir <- file.path("data", "processed")
}
dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)

source("scripts/_immune_markers.R")
danaher_genes <- unique(unlist(immune_markers))

.log <- character()

# --- Load packages -----------------------------------------------------------
for (pkg in c("cBioPortalData", "GEOquery", "dplyr", "tidyr", "tibble")) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% c("cBioPortalData", "GEOquery")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

cbio <- cBioPortal()

# =============================================================================
# Helper: fetch clinical + expression for a study
# =============================================================================
fetch_cohort <- function(
  cbio,
  study_id,
  mrna_profile_id,
  label,
  proc_dir,
  danaher_genes,
  is_linear_scale = FALSE
) {
  log <- character()

  clin_file <- file.path(proc_dir, paste0(label, "_clinical.csv"))
  expr_file <- file.path(proc_dir, paste0(label, "_expression.csv"))

  if (all(file.exists(c(clin_file, expr_file)))) {
    log <- c(
      log,
      sprintf("[%s] Using cached files (delete to re-fetch)", label)
    )
    return(log)
  }

  log <- c(
    log,
    sprintf("[%s] Fetching from cBioPortal (study: %s)", label, study_id)
  )

  # ---- Clinical data --------------------------------------------------------
  clin_raw <- clinicalData(cbio, studyId = study_id)
  log <- c(
    log,
    sprintf("  Raw clinical: %d rows x %d cols", nrow(clin_raw), ncol(clin_raw))
  )

  # Standardise column names across cohorts
  # METABRIC uses CLAUDIN_SUBTYPE, TCGA uses SUBTYPE
  pam50_col <- if ("SUBTYPE" %in% names(clin_raw)) {
    "SUBTYPE"
  } else if ("CLAUDIN_SUBTYPE" %in% names(clin_raw)) {
    "CLAUDIN_SUBTYPE"
  } else {
    NA_character_
  }

  os_status_col <- if ("OS_STATUS" %in% names(clin_raw)) {
    "OS_STATUS"
  } else {
    NA_character_
  }
  os_months_col <- if ("OS_MONTHS" %in% names(clin_raw)) {
    "OS_MONTHS"
  } else {
    NA_character_
  }

  clinical_clean <- clin_raw %>%
    transmute(
      patient_id = patientId,
      cohort = label,
      age_at_diagnosis = if ("AGE_AT_DIAGNOSIS" %in% names(.)) {
        as.numeric(AGE_AT_DIAGNOSIS)
      } else if ("AGE" %in% names(.)) {
        as.numeric(AGE)
      } else {
        NA_real_
      },
      sex = if ("SEX" %in% names(.)) SEX else NA_character_,
      er_status = if ("ER_IHC" %in% names(.)) {
        ER_IHC
      } else if ("ER_STATUS_BY_IHC" %in% names(.)) {
        ER_STATUS_BY_IHC
      } else {
        NA_character_
      },
      her2_status = if ("HER2_SNP6" %in% names(.)) {
        HER2_SNP6
      } else if ("HER2_STATUS_BY_IHC" %in% names(.)) {
        HER2_STATUS_BY_IHC
      } else {
        NA_character_
      },
      pam50_subtype = if (!is.na(pam50_col)) .[[pam50_col]] else NA_character_,
      os_months = if (!is.na(os_months_col)) {
        as.numeric(.[[os_months_col]])
      } else {
        NA_real_
      },
      os_status = if (!is.na(os_status_col)) {
        .[[os_status_col]]
      } else {
        NA_character_
      }
    )

  write.csv(clinical_clean, clin_file, row.names = FALSE)
  log <- c(
    log,
    sprintf("  %s_clinical.csv: %d patients", label, nrow(clinical_clean))
  )

  # ---- Expression data (Danaher marker genes only) ---------------------------
  expr_raw <- getDataByGenes(
    cbio,
    studyId = study_id,
    genes = danaher_genes,
    by = "hugoGeneSymbol",
    molecularProfileIds = mrna_profile_id
  )

  df <- expr_raw[[mrna_profile_id]]
  if (is.null(df) || nrow(df) == 0) {
    stop(sprintf(
      "[%s] No expression data returned for profile %s",
      label,
      mrna_profile_id
    ))
  }

  genes_returned <- unique(df$hugoGeneSymbol)
  genes_missing <- setdiff(danaher_genes, genes_returned)
  log <- c(
    log,
    sprintf(
      "  Expression: %d genes x %d samples",
      length(genes_returned),
      length(unique(df$sampleId))
    )
  )
  if (length(genes_missing) > 0) {
    log <- c(
      log,
      sprintf("  Missing genes: %s", paste(genes_missing, collapse = ", "))
    )
  }

  expr_wide <- df %>%
    select(hugoGeneSymbol, sampleId, value) %>%
    pivot_wider(
      names_from = sampleId,
      values_from = value,
      values_fn = ~ mean(.x, na.rm = TRUE)
    )

  # Log2-transform RNA-seq RSEM counts (microarray data already log2)
  if (is_linear_scale) {
    gene_col <- colnames(expr_wide)[1]
    expr_mat <- expr_wide %>%
      column_to_rownames(var = gene_col) %>%
      as.matrix()
    # RSEM values can be 0; add pseudocount of 1
    expr_mat <- log2(expr_mat + 1)
    expr_wide <- expr_mat %>%
      as.data.frame() %>%
      rownames_to_column(var = gene_col)
    log <- c(log, "  Applied log2(x+1) transformation to RNA-seq counts")
  }

  write.csv(expr_wide, expr_file, row.names = FALSE)
  log <- c(log, sprintf("  %s_expression.csv written", label))

  log
}

# =============================================================================
# 1. METABRIC
# =============================================================================
.log <- c(
  .log,
  "",
  tryCatch(
    fetch_cohort(
      cbio,
      "brca_metabric",
      "brca_metabric_mrna",
      "metabric",
      proc_dir,
      danaher_genes,
      is_linear_scale = FALSE
    ),
    error = function(e) sprintf("[metabric] FAILED: %s", conditionMessage(e))
  )
)

# =============================================================================
# 2. TCGA-BRCA (Pan-Cancer Atlas)
# =============================================================================
# First, discover available mRNA profiles for TCGA-BRCA
tcga_study <- "brca_tcga_pan_can_atlas_2018"
tcga_profiles <- molecularProfiles(cbio, studyId = tcga_study)
mrna_profiles <- tcga_profiles %>%
  as.data.frame() %>%
  filter(grepl("mrna|rna", molecularProfileId, ignore.case = TRUE))

.log <- c(
  .log,
  "",
  sprintf("[tcga] Available mRNA profiles for %s:", tcga_study),
  sprintf("  %s (%s)", mrna_profiles$molecularProfileId, mrna_profiles$name)
)

# Use the base RNA-seq profile (not z-scores) for Danaher scoring
tcga_mrna_id <- mrna_profiles$molecularProfileId[
  grepl("rna_seq.*mrna$", mrna_profiles$molecularProfileId) |
    grepl("rna_seq_v2_mrna$", mrna_profiles$molecularProfileId)
]
if (length(tcga_mrna_id) == 0) {
  # Fallback: take first non-zscore profile
  tcga_mrna_id <- mrna_profiles$molecularProfileId[
    !grepl(
      "z_scores|zscores",
      mrna_profiles$molecularProfileId,
      ignore.case = TRUE
    )
  ][1]
}
tcga_mrna_id <- tcga_mrna_id[1] # take first match
.log <- c(.log, sprintf("[tcga] Selected profile: %s", tcga_mrna_id))

.log <- c(
  .log,
  "",
  tryCatch(
    fetch_cohort(
      cbio,
      tcga_study,
      tcga_mrna_id,
      "tcga",
      proc_dir,
      danaher_genes,
      is_linear_scale = TRUE
    ),
    error = function(e) sprintf("[tcga] FAILED: %s", conditionMessage(e))
  )
)

# =============================================================================
# 3. I-SPY2 (from GEO)
# =============================================================================
.log <- c(.log, "", "--- I-SPY2 (from GEO) ---")

raw_dir <- file.path("data", "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

# The gene-level file will be inside a GSE194040 subdirectory
ispy_dir <- file.path(raw_dir, "GSE194040")
ispy_file <- file.path(
  ispy_dir,
  "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz"
)

if (!file.exists(ispy_file)) {
  .log <- c(.log, "[ispy2] Downloading gene-level expression from GEO...")

  # Download only the gene-level file (47 MB, not 670 MB!)
  getGEOSuppFiles(
    "GSE194040",
    baseDir = raw_dir,
    filter_regex = "meanCol_geneLevel",
    makeDirectory = TRUE
  )

  .log <- c(
    .log,
    sprintf(
      "  Downloaded: %s (%.1f MB)",
      basename(ispy_file),
      file.size(ispy_file) / 1024^2
    )
  )
} else {
  .log <- c(
    .log,
    sprintf(
      "[ispy2] Using cached file (%.1f MB)",
      file.size(ispy_file) / 1024^2
    )
  )
}


# =============================================================================
# Verify outputs
# =============================================================================
.log <- c(.log, "", "--- File check ---")
for (f in list.files(proc_dir, pattern = "\\.csv$")) {
  fp <- file.path(proc_dir, f)
  .log <- c(
    .log,
    sprintf("  %s (%s)", f, format(file.size(fp), big.mark = ","))
  )
}

message(paste(.log, collapse = "\n"))
