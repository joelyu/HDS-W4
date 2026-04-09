# =============================================================================
# 00-fetch-data.R
# Fetch METABRIC + TCGA-BRCA (cBioPortal) + I-SPY2 (GEO) → data/processed/
#
# Three cohorts:
#   1. METABRIC (brca_metabric) — Illumina HT-12 v3, log2 intensity
#   2. TCGA-BRCA (brca_tcga_pan_can_atlas_2018) — RNA-Seq V2, RSEM
#   3. I-SPY2 (GSE194040) — Agilent microarray, log2 intensity
#
# Output per cohort:
#   data/processed/{cohort}_clinical.csv
#   data/processed/{cohort}_expression.csv   (gene × patient, log2 scale)
#
# After processing, all three expression files are harmonized to the
# intersection of Danaher genes present on all platforms (57/60 genes).
# Excluded: TPSB2, XCL2, KIR2DL3 (platform-specific gaps).
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

# Skip cBioPortal API initialisation if all cached files already exist —
# avoids hard failure when the cBioPortal service is down (HTTP 503).
required_files <- c(
  file.path(proc_dir, "metabric_clinical.csv"),
  file.path(proc_dir, "metabric_clinical_extended.csv"),
  file.path(proc_dir, "metabric_expression.csv"),
  file.path(proc_dir, "tcga_clinical.csv"),
  file.path(proc_dir, "tcga_clinical_extended.csv"),
  file.path(proc_dir, "tcga_expression.csv")
)

if (all(file.exists(required_files))) {
  .log <- c(.log, "All cBioPortal cache files present — skipping API init")
  cbio <- NULL
} else {
  cbio <- cBioPortal()
}

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
  ext_file  <- file.path(proc_dir, paste0(label, "_clinical_extended.csv"))

  if (all(file.exists(c(clin_file, expr_file, ext_file)))) {
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

  # ---- Extended clinical (study-specific confounders) ------------------------
  if (label == "metabric") {
    ext_clean <- clin_raw %>%
      transmute(
        patient_id      = patientId,
        tumor_stage     = as.integer(suppressWarnings(as.numeric(TUMOR_STAGE))),
        grade           = as.integer(suppressWarnings(as.numeric(GRADE))),
        tumor_size      = as.numeric(TUMOR_SIZE),
        lymph_nodes_pos = as.integer(suppressWarnings(as.numeric(LYMPH_NODES_EXAMINED_POSITIVE))),
        npi             = as.numeric(NPI),
        chemotherapy    = if ("CHEMOTHERAPY"    %in% names(.)) CHEMOTHERAPY    else NA_character_,
        hormone_therapy = if ("HORMONE_THERAPY" %in% names(.)) HORMONE_THERAPY else NA_character_
      )
  } else if (label == "tcga") {
    ext_clean <- clin_raw %>%
      transmute(
        patient_id      = patientId,
        ajcc_stage      = AJCC_PATHOLOGIC_TUMOR_STAGE,
        # Order matters: check most specific first to avoid "STAGE I" matching "STAGE II/III"
        tumor_stage_num = case_when(
          grepl("^STAGE IV$",        AJCC_PATHOLOGIC_TUMOR_STAGE) ~ 4L,
          grepl("^STAGE III[ABC]?$", AJCC_PATHOLOGIC_TUMOR_STAGE) ~ 3L,
          grepl("^STAGE II[AB]?$",   AJCC_PATHOLOGIC_TUMOR_STAGE) ~ 2L,
          grepl("^STAGE I[AB]?$",    AJCC_PATHOLOGIC_TUMOR_STAGE) ~ 1L,
          TRUE ~ NA_integer_
        ),
        path_t_stage = PATH_T_STAGE,
        radiation    = if ("RADIATION_THERAPY" %in% names(.)) RADIATION_THERAPY else NA_character_
      )
  } else {
    ext_clean <- NULL
  }
  if (!is.null(ext_clean)) {
    write.csv(ext_clean, ext_file, row.names = FALSE)
    log <- c(log, sprintf("  %s_clinical_extended.csv: %d rows", label, nrow(ext_clean)))
  }

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

# Skip the molecular profile discovery API call if cache files already exist
# (cbio is NULL when all required files were already on disk)
if (!is.null(cbio)) {
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
} else {
  # All cached — use the known profile id
  tcga_mrna_id <- "brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
}

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

# ---- Process I-SPY2 expression (Danaher genes only) -------------------------
ispy_expr_file <- file.path(proc_dir, "ispy2_expression.csv")

if (!file.exists(ispy_expr_file)) {
  .log <- c(.log, "[ispy2] Processing expression data...")

  ispy_raw <- read.delim(gzfile(ispy_file), header = TRUE, row.names = 1,
                         check.names = FALSE)
  .log <- c(.log, sprintf("  Raw: %d genes x %d samples",
                           nrow(ispy_raw), ncol(ispy_raw)))

  # Rename KIAA0125 → FAM30A (Agilent uses old symbol)
  if ("KIAA0125" %in% rownames(ispy_raw) && !"FAM30A" %in% rownames(ispy_raw)) {
    rownames(ispy_raw)[rownames(ispy_raw) == "KIAA0125"] <- "FAM30A"
    .log <- c(.log, "  Renamed KIAA0125 → FAM30A")
  }

  # Filter to Danaher genes
  available <- danaher_genes[danaher_genes %in% rownames(ispy_raw)]
  missing   <- setdiff(danaher_genes, available)
  .log <- c(.log, sprintf("  Danaher genes: %d/%d found", length(available),
                           length(danaher_genes)))
  if (length(missing) > 0) {
    .log <- c(.log, sprintf("  Missing genes: %s", paste(missing, collapse = ", ")))
  }

  # Extract and format like METABRIC/TCGA (hugoGeneSymbol as first column)
  ispy_expr <- ispy_raw[available, , drop = FALSE]
  ispy_expr <- ispy_expr %>%
    as.data.frame() %>%
    tibble::rownames_to_column("hugoGeneSymbol")

  write.csv(ispy_expr, ispy_expr_file, row.names = FALSE)
  .log <- c(.log, sprintf("  ispy2_expression.csv: %d genes x %d patients",
                           nrow(ispy_expr), ncol(ispy_expr) - 1))
} else {
  .log <- c(.log, "[ispy2] Using cached ispy2_expression.csv")
}

# ---- I-SPY2 clinical/phenotype (from GEO series matrix) ---------------------
ispy_clin_file <- file.path(proc_dir, "ispy2_clinical.csv")

if (!file.exists(ispy_clin_file)) {
  .log <- c(.log, "[ispy2] Extracting clinical data from GEO series matrix...")

  raw_dir <- file.path("data", "raw")
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

  gse <- getGEO("GSE194040", destdir = raw_dir, GSEMatrix = TRUE)
  pd_all <- bind_rows(pData(gse[[1]]), pData(gse[[2]])) %>%
    transmute(
      patient_id = `patient id:ch1`,
      arm        = `arm:ch1`,
      pcr        = as.integer(`pcr:ch1`),
      hr         = as.integer(`hr:ch1`),
      her2       = as.integer(`her2:ch1`)
    )

  write.csv(pd_all, ispy_clin_file, row.names = FALSE)
  .log <- c(.log, sprintf("  ispy2_clinical.csv: %d patients", nrow(pd_all)))
} else {
  .log <- c(.log, "[ispy2] Using cached ispy2_clinical.csv")
}

# =============================================================================
# 4. Three-way gene harmonization
# =============================================================================
# Restrict all three expression datasets to genes present on all platforms.
# Excludes 3 genes: TPSB2 (absent METABRIC + I-SPY2), XCL2 (absent METABRIC),
# KIR2DL3 (absent I-SPY2). Yields 57/60 Danaher genes.
# Cytotoxic cells (focal predictor): 10/10 markers on all platforms.
# =============================================================================
.log <- c(.log, "", "--- Three-way gene harmonization ---")

expr_files <- c(
  metabric = file.path(proc_dir, "metabric_expression.csv"),
  tcga     = file.path(proc_dir, "tcga_expression.csv"),
  ispy2    = file.path(proc_dir, "ispy2_expression.csv")
)

expr_list <- lapply(expr_files, read.csv, check.names = FALSE)
gene_col  <- "hugoGeneSymbol"

common_genes <- Reduce(intersect, lapply(expr_list, function(df) df[[gene_col]]))
excluded     <- setdiff(danaher_genes, common_genes)

.log <- c(.log, sprintf("  Common genes: %d/%d", length(common_genes),
                         length(danaher_genes)))
.log <- c(.log, sprintf("  Excluded: %s", paste(excluded, collapse = ", ")))

for (nm in names(expr_list)) {
  n_before <- nrow(expr_list[[nm]])
  expr_list[[nm]] <- expr_list[[nm]][expr_list[[nm]][[gene_col]] %in% common_genes, ]
  write.csv(expr_list[[nm]], expr_files[[nm]], row.names = FALSE)
  .log <- c(.log, sprintf("  %s: %d → %d genes", nm, n_before,
                           nrow(expr_list[[nm]])))
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
