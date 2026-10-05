#======================================================================================================================================
#======================================================================================================================================
# Figure S15 (2.0) - Simple model
#======================================================================================================================================
# Sensitivity analysis with Delta_FEV1 with Delta_%ofpredicted(PCT) - WITH adjustment for FEV1 % of predicted
#======================================================================================================================================
# === CHECK MISSING VALUES ===
sum(is.na(ORACLE_Placebo_Lung_imputed$FEV1_preBD_PCT_0W))
sum(is.na(ORACLE_Placebo_Lung_imputed$FEV1_preBD_PCT_52W))

# === CREATE DATASET WITH DELTA FEV1 ===
ORACLE_Placebo_Lung_imputed_PCT <- ORACLE_Placebo_Lung_imputed |>
  mutate(
    delta_FEV1_preBD_PCT = FEV1_preBD_PCT_52W - FEV1_preBD_PCT_0W
  )

# === VERIFY ===
summary(ORACLE_Placebo_Lung_imputed_PCT$delta_FEV1_preBD_PCT)
sum(is.na(ORACLE_Placebo_Lung_imputed_PCT$delta_FEV1_preBD_PCT))

# === UNIVARIATE MODEL - PCT ===

Univariate_deltalung_model_PCT <- MI_estimates(
  data        = ORACLE_Placebo_Lung_imputed_PCT,
  outcome_var = "delta_FEV1_preBD_PCT",
  predictor_vars = c(
    "Age_per_10",                          # 1  - Demographic
    "Sex",                                 # 2  - Demographic
    "BMI_per_5",                           # 3  - Demographic
    "Smoking_history_yes_no",              # 4  - Demographic
    "Airborne_allergen_sensitisation",     # 5  - Comorbidities
    "Allergic_Rhinitis",                   # 6  - Comorbidities
    "Eczema",                              # 7  - Comorbidities
    "CRSwNP",                              # 8  - Comorbidities
    "CRSsNP",                              # 9  - Comorbidities
    "Psychiatric_disease",                 # 10 - Comorbidities
    "GINA_step_numeric",                   # 11 - Asthma history
    "ACQ_score_0W",                        # 12 - Asthma history
    "Attack_history_num",                  # 13 - Asthma history
    "FEV1_preBD_PCT_0W_per_10",             # 14 - Baseline lung function
    "FEV1_reversibility_0W_per_10",        # 15 - Baseline lung function
    "BEC_log10",                           # 16 - Inflammatory biomarkers
    "FeNO_log10",                          # 17 - Inflammatory biomarkers
    "IgE_log10"                            # 18 - Inflammatory biomarkers
  ),
  covariables          = c("FEV1_preBD_PCT_0W_per_10"),
  imp_col              = ".imp",
  followup_offset      = "No",
  followup_col         = "Follow_up_duration_days",
  random_intercept_var = "Enrolled_Trial_name",
  model_type           = "lm"
)

# === EXTRACT & LABEL RESULTS ===
Univariate_deltalung_results_PCT <- attr(Univariate_deltalung_model_PCT, "combined_results")

row.names(Univariate_deltalung_results_PCT) <- c(
  # Demographic (4)
  "Age (per 10 years)",
  "Sex (Male vs. Female)",
  "BMI (per 5 kg/m\u00b2)",
  "Smoking history",
  # Comorbidities (6)
  "Airborne allergen sensitisation",
  "Allergic rhinitis",
  "Eczema",
  "CRS with nasal polyps",
  "CRS without nasal polyps",
  "Psychiatric disease",
  # Asthma history (3)
  "GINA step",
  "ACQ score (baseline)",
  "Exacerbation history (count)",
  # Baseline lung function (2)
  "FEV1 pre-BD (% of predicted, per 10% decrease)",
  "FEV1 reversibility (per 10% increase)",
  # Inflammatory biomarkers (3)
  "Blood eosinophils (log10)",
  "FeNO (log10)",
  "Total IgE (log10)"
)

# === FOREST PLOT - UNIVARIATE PCT ===
pdf(
  "/Users/joel/ORACLE - Placebo - Exacerbations/Project/Sup_Forplot_Lung_PCT_Univar_final.pdf",
  width = 18, height = 12
)
par(mar = c(4, 2, 2, 2))

Sup_Forplot_Lung_PCT_Univar_final <- forplo(
  as.data.frame(Univariate_deltalung_results_PCT[, c("estimate", "2.5 %", "97.5 %")]),
  xlim            = c(-10, 10),
  em              = "aRC",
  linreg          = TRUE,
  row.labels      = row.names(Univariate_deltalung_results_PCT),
  left.align      = FALSE,
  add.arrow.right = FALSE,
  arrow.right.length = 20,
  add.arrow.left  = FALSE,
  arrow.left.length  = 20,
  left.bar        = FALSE,
  shade.every     = 1,
  shade.col       = "grey",
  shade.alpha     = 0.2,
  margin.left  = 28,   # ↑ more room for long left labels
  margin.right = 14,   # balanced — 20% less space on right
  groups          = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 3)),
  grouplabs       = c(
    "Demographic",
    "Comorbidities",
    "Asthma history",
    "Baseline lung function",
    "Inflammatory biomarkers"
  ),
  char = 20,
  size = 1.5,
  col  = c(
    rep("darkgreen",  4),
    rep("purple",     6),
    rep("darkred",    3),
    rep("darkblue",   2),
    rep("darkorange", 3)
  )
)
dev.off()

#======================================================================================================================================

# === MULTIVARIATE MODEL - PCT ===

Multivariate_deltalung_model_PCT <- MI_estimates(
  data        = ORACLE_Placebo_Lung_imputed_PCT,
  outcome_var = "delta_FEV1_preBD_PCT",
  predictor_vars = c(
    "Age_per_10",                          # 1  - Demographic
    "Sex",                                 # 2  - Demographic
    "BMI_per_5",                           # 3  - Demographic
    "Smoking_history_yes_no",              # 4  - Demographic
    "Airborne_allergen_sensitisation",     # 5  - Comorbidities
    "Allergic_Rhinitis",                   # 6  - Comorbidities
    "Eczema",                              # 7  - Comorbidities
    "CRSwNP",                              # 8  - Comorbidities
    "CRSsNP",                              # 9  - Comorbidities
    "Psychiatric_disease",                 # 10 - Comorbidities
    "GINA_step_numeric",                   # 11 - Asthma history
    "ACQ_score_0W",                        # 12 - Asthma history
    "Attack_history_num",                  # 13 - Asthma history
    "FEV1_preBD_PCT_0W_per_10",             # 14 - Baseline lung function
    "FEV1_reversibility_0W_per_10",        # 15 - Baseline lung function
    "BEC_log10",                           # 16 - Inflammatory biomarkers
    "FeNO_log10",                          # 17 - Inflammatory biomarkers
    "IgE_log10"                            # 18 - Inflammatory biomarkers
  ),
  covariables = c(
    "Age_per_10",
    "Sex",
    "BMI_per_5",
    "FEV1_preBD_PCT_0W_per_10",
    "FEV1_reversibility_0W_per_10",
    "FeNO_log10",
    "BEC_log10",
    "IgE_log10"
  ),
  imp_col              = ".imp",
  followup_offset      = "No",
  followup_col         = "Follow_up_duration_days",
  random_intercept_var = "Enrolled_Trial_name",
  model_type           = "lm"
)

# === EXTRACT & LABEL RESULTS ===
Multivariate_deltalung_results_PCT <- attr(Multivariate_deltalung_model_PCT, "combined_results")

row.names(Multivariate_deltalung_results_PCT) <- c(
  # Demographic (4)
  "Age (per 10 years)",
  "Sex (Male vs. Female)",
  "BMI (per 5 kg/m\u00b2)",
  "Smoking history",
  # Comorbidities (6)
  "Airborne allergen sensitisation",
  "Allergic rhinitis",
  "Eczema",
  "CRS with nasal polyps",
  "CRS without nasal polyps",
  "Psychiatric disease",
  # Asthma history (3)
  "GINA step",
  "ACQ score (baseline)",
  "Exacerbation history (count)",
  # Baseline lung function (2)
  "FEV1 pre-BD (% of predicted, per 10% decrease)",
  "FEV1 reversibility (per 10% increase)",
  # Inflammatory biomarkers (3)
  "Blood eosinophils (log10)",
  "FeNO (log10)",
  "Total IgE (log10)"
)

# === FOREST PLOT - MULTIVARIATE PCT ===
pdf(
  "/Users/joel/ORACLE - Placebo - Exacerbations/Project/Sup_Forplot_Lung_PCT_Multivar_final.pdf",
  width = 18, height = 12
)
par(mar = c(4, 2, 2, 2))

Sup_Forplot_Lung_PCT_Multivar_final <- forplo(
  as.data.frame(Multivariate_deltalung_results_PCT[, c("estimate", "2.5 %", "97.5 %")]),
  xlim            = c(-10, 10),
  em              = "aRC",
  linreg          = TRUE,
  row.labels      = row.names(Multivariate_deltalung_results_PCT),
  left.align      = FALSE,
  add.arrow.right = FALSE,
  arrow.right.length = 20,
  add.arrow.left  = FALSE,
  arrow.left.length  = 20,
  left.bar        = FALSE,
  shade.every     = 1,
  shade.col       = "grey",
  shade.alpha     = 0.2,
  margin.left  = 28,   # ↑ more room for long left labels
  margin.right = 14,   # balanced — 20% less space on right
  groups          = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 3)),
  grouplabs       = c(
    "Demographic",
    "Comorbidities",
    "Asthma history",
    "Baseline lung function",
    "Inflammatory biomarkers"
  ),
  char = 20,
  size = 1.5,
  col  = c(
    rep("darkgreen",  4),
    rep("purple",     6),
    rep("darkred",    3),
    rep("darkblue",   2),
    rep("darkorange", 3)
  )
)
dev.off()