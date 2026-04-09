# =============================================================================
# 05-ispy2.R
# Test: does cytotoxic effector score predict immunotherapy response (pCR)
# in the I-SPY2 pembrolizumab arm?
#
# Data: GSE194040 (I-SPY2 TRIAL, 988 breast cancer patients)
#
# Structure:
#   A. Primary analysis — pre-specified cytotoxic score, HR-adjusted, pembro arm
#   B. Discovery screen — all 14 cell types, HR-adjusted, FDR corrected
#   C. Control arm — cytotoxic in paclitaxel-only arm (specificity check)
#   D. Interaction test — cytotoxic × pembro arm
#
# Adjustment: HR status only. HER2 dropped because both pembro and control
# arms are HER2-negative by I-SPY2 trial design (all her2=0), so it cannot
# vary as a covariate. PAM50, age, and stage are not available in I-SPY2.
#
# Expects: scores_ispy2 from 01-immune-scoring.R
#          proc_dir defined by parent qmd
# =============================================================================

library(dplyr)

if (!exists("proc_dir")) proc_dir <- file.path("data", "processed")

# --- Setup: load scores and merge with phenotype ------------------------------
if (!exists("scores_ispy2")) {
  scores_ispy2 <- read.csv(file.path(proc_dir, "ispy2_immune_scores.csv"))
}
message(sprintf("=== I-SPY2 immune scores: %d patients ===", nrow(scores_ispy2)))

pd_all <- read.csv(file.path(proc_dir, "ispy2_clinical.csv"))
pd_all$patient_id <- as.character(pd_all$patient_id)

merged <- inner_join(scores_ispy2, pd_all, by = "patient_id")

message("\nArm label distribution (sanity check):")
print(sort(table(merged$arm), decreasing = TRUE))

n_pembro <- sum(grepl("Pembrolizumab", merged$arm, ignore.case = TRUE))
n_pcr_pembro <- sum(merged$pcr[grepl("Pembrolizumab", merged$arm, ignore.case = TRUE)])
message(sprintf("\n=== Pembrolizumab arm: %d patients ===", n_pembro))
message(sprintf("pCR rate: %d/%d = %.1f%%", n_pcr_pembro, n_pembro,
                100 * n_pcr_pembro / n_pembro))

# --- Build pooled (pembro + control) sample with shared scaling --------------
# All cell type scores are scaled ONCE on the pooled sample so that "1 SD"
# means the same thing in every model below. This makes per-arm ORs (Sections
# A and C) directly comparable AND consistent with the interaction model
# (Section D), avoiding the apples-to-oranges problem of arm-specific scaling.
all_cell_types <- setdiff(
  colnames(scores_ispy2),
  c("patient_id", "cohort", "T_cells")  # exclude T_cells (parent of CD8)
)

pooled <- merged %>%
  filter(grepl("Pembrolizumab", arm, ignore.case = TRUE) | arm == "Paclitaxel") %>%
  mutate(is_pembro = ifelse(grepl("Pembrolizumab", arm), 1, 0))

# Scale every cell type once on the pooled sample
for (ct in all_cell_types) {
  pooled[[paste0(ct, "_z")]] <- scale(pooled[[ct]])[, 1]
}

pembro_pooled <- pooled %>% filter(is_pembro == 1)
control_pooled <- pooled %>% filter(is_pembro == 0)

message(sprintf("\nPooled sample: %d (pembro=%d, control=%d) — scaled once",
                nrow(pooled), nrow(pembro_pooled), nrow(control_pooled)))

# --- Helper for logistic OR extraction ----------------------------------------
# Operates on the pooled-scaled `_z` columns so the OR is "per 1 pooled SD".
fit_logistic <- function(df, cell_type, formula_extra = "") {
  z_col <- paste0(cell_type, "_z")
  formula <- as.formula(paste("pcr ~", z_col, formula_extra))
  fit <- glm(formula, data = df, family = binomial)
  s <- summary(fit)
  or <- exp(coef(fit)[2])
  ci <- exp(confint.default(fit)[2, ])
  data.frame(
    cell_type = cell_type,
    OR = or,
    OR_lower = ci[1],
    OR_upper = ci[2],
    p_value = s$coefficients[2, 4],
    n = nrow(df),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# A. PRIMARY ANALYSIS — pre-specified cytotoxic score (HR-adjusted)
# =============================================================================
# Cytotoxic was pre-specified as the focal predictor based on the survival
# analysis (METABRIC + TCGA). One hypothesis, one test, no FDR needed.
# =============================================================================
message("\n=== A. PRIMARY: Cytotoxic → pCR (HR-adjusted, pembro arm) ===")

primary <- fit_logistic(pembro_pooled, "Cytotoxic_cells", "+ factor(hr)")
message(sprintf("  Cytotoxic_cells: OR=%.2f [%.2f-%.2f] p=%.4f  n=%d",
                primary$OR, primary$OR_lower, primary$OR_upper,
                primary$p_value, primary$n))

# =============================================================================
# B. DISCOVERY SCREEN — all 14 cell types with FDR correction
# =============================================================================
# Exploratory: which immune cell types predict pCR in the pembro arm?
# Tests all 14 Danaher cell types simultaneously, controls FDR via BH.
# =============================================================================
message("\n=== B. DISCOVERY SCREEN: all 14 cell types (HR-adjusted, FDR corrected) ===")

discovery <- lapply(all_cell_types, function(ct) {
  fit_logistic(pembro_pooled, ct, "+ factor(hr)")
}) %>%
  bind_rows() %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

message(sprintf("  %d/%d cell types nominally significant (p<0.05)",
                sum(discovery$p_value < 0.05), nrow(discovery)))
message(sprintf("  %d/%d cell types FDR-significant (FDR<0.05)",
                sum(discovery$fdr < 0.05), nrow(discovery)))

for (i in seq_len(nrow(discovery))) {
  r <- discovery[i, ]
  sig <- ifelse(r$fdr < 0.05, "***",
         ifelse(r$p_value < 0.05, "*",
         ifelse(r$p_value < 0.1, ".", "")))
  message(sprintf("  %-20s OR=%.2f [%.2f-%.2f] p=%.4f FDR=%.4f %s",
                  r$cell_type, r$OR, r$OR_lower, r$OR_upper,
                  r$p_value, r$fdr, sig))
}

# =============================================================================
# C. CONTROL ARM — cytotoxic in paclitaxel-only (specificity check)
# =============================================================================
# If cytotoxic predicts pCR specifically because of pembrolizumab, the effect
# should be weaker in the chemo-only control arm. If it's equally strong in
# both, the signal is general prognosis, not immunotherapy-specific.
# =============================================================================
message("\n=== C. CONTROL: Cytotoxic → pCR (HR-adjusted, paclitaxel-only arm) ===")
message(sprintf("  Control: %d patients, pCR rate: %.1f%%",
                nrow(control_pooled), 100 * mean(control_pooled$pcr)))

control_cyto <- fit_logistic(control_pooled, "Cytotoxic_cells", "+ factor(hr)")
message(sprintf("  Cytotoxic_cells: OR=%.2f [%.2f-%.2f] p=%.4f  n=%d",
                control_cyto$OR, control_cyto$OR_lower, control_cyto$OR_upper,
                control_cyto$p_value, control_cyto$n))

# =============================================================================
# D. INTERACTION — cytotoxic × pembrolizumab arm
# =============================================================================
# Tests whether cytotoxic's effect on pCR DIFFERS between pembro and control.
# Significant interaction = treatment-specific effect (pembro amplifies it).
# Uses the same pooled-scaled cyto_z as Sections A and C — all three are now
# directly comparable.
# =============================================================================
message("\n=== D. INTERACTION: Cytotoxic × Pembrolizumab arm ===")

fi <- glm(pcr ~ Cytotoxic_cells_z * factor(is_pembro) + factor(hr),
          data = pooled, family = binomial)
si <- summary(fi)
message(sprintf("  Pooled n=%d (%d pembro, %d control)",
                nrow(pooled), sum(pooled$is_pembro), sum(!pooled$is_pembro)))
print(round(si$coefficients, 4))

# Extract interaction term explicitly
interaction_term <- "Cytotoxic_cells_z:factor(is_pembro)1"
interaction_p <- si$coefficients[interaction_term, "Pr(>|z|)"]
interaction_est <- si$coefficients[interaction_term, "Estimate"]
message(sprintf("\n  Interaction term: estimate=%.3f, p=%.4f",
                interaction_est, interaction_p))
if (interaction_p < 0.05) {
  message("  → Cytotoxic effect significantly differs between arms")
} else {
  message("  → No significant interaction (underpowered with n=69 vs n=control)")
}

# =============================================================================
# 6. Save outputs
# =============================================================================
# Pembrolizumab arm patient-level scores (for plotting) — uses pooled-scaled
# subset so plotted scores match the model inputs above
pembro_out <- pembro_pooled %>%
  select(patient_id, arm, pcr, hr, her2, all_of(all_cell_types))
write.csv(pembro_out, file.path(proc_dir, "ispy2_pembro_scores.csv"), row.names = FALSE)

# Primary result (cytotoxic, pre-specified)
primary_out <- bind_rows(
  primary %>% mutate(arm = "Pembrolizumab", model = "Primary (HR-adjusted)"),
  control_cyto %>% mutate(arm = "Paclitaxel (control)", model = "Control (HR-adjusted)")
)
write.csv(primary_out, file.path(proc_dir, "ispy2_pembro_logistic.csv"), row.names = FALSE)

# Discovery screen (all 14 cell types with FDR)
write.csv(discovery, file.path(proc_dir, "ispy2_discovery_screen.csv"), row.names = FALSE)

# Interaction model coefficients
write.csv(as.data.frame(si$coefficients) %>% tibble::rownames_to_column("term"),
          file.path(proc_dir, "ispy2_interaction.csv"), row.names = FALSE)

message("\n=== Done ===")
message(sprintf("Outputs: %s/ispy2_pembro_scores.csv", proc_dir))
message(sprintf("         %s/ispy2_pembro_logistic.csv (primary + control)", proc_dir))
message(sprintf("         %s/ispy2_discovery_screen.csv (14 cell types + FDR)", proc_dir))
message(sprintf("         %s/ispy2_interaction.csv", proc_dir))
