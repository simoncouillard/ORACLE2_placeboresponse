#######################################################################################################################################
#######################################################################################################################################

### SCRIPT FOR THE MAIN ANALYSIS (IN THE MANUSCRIPT - TABLE 1 + FIGURE 1,2,3) FOR THE PLACEBO-ASTHMA ANALYSIS
#JSP and SML et al.

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# PART A: DATAFRAME PREPARATION AND TABLE 1
#######################################################################################################################################
#######################################################################################################################################
# Remove everything in the environment 
rm(list = ls())

# After removing large objects, prompt R to release memory back to the OS
gc()
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 1.	Loading the packages
#--------------------------------------------------------------------------------------------------------------------------------------
# 1.1.	From library

# Data manipulation
library(tidyverse)   # includes dplyr, tidyr, purrr, stringr, forcats, tibble, readr, ggplot2, lubridate
library(magrittr)
library(readxl)

# Statistical analysis
library(ggsci)
library(ggVennDiagram)
library(MASS)
library(mice)
library(mitools)
library(rms)
library(MuMIn)
library(boot)
library(AICcmodavg)
library(caret)
library(gam)
library(rstatix)
library(metafor)
library(lme4)
library(glmmTMB)
library(leaps)

# Survival analysis
library(survival)
library(ggsurvfit)
library(tidycmprsk)

# Visualization
library(ggpubr)
library(ggvenn)
library(cowplot)
library(pheatmap)
library(circlize)
library(forplo)
library(patchwork)
library(gridExtra)
library(grid)
library(DiagrammeR)
library(rsvg)

# Reporting & tables
library(gtsummary)
library(tableone)
library(table1)
library(Gmisc, quietly = TRUE)
library(glue)
library(knitr)
library(flextable)

# Utilities
library(writexl)
library(rspiro)
library(remotes)
library(coxme)
#--------------------------------------------------------------------------------------------------------------------------------------
# 1.2.	Install and load the homemade package MIAnalysis
#install.packages("devtools")
library(devtools)
required_packages <- c("dplyr", "MASS", "Hmisc", "mice", "rlang", "survival")
devtools::install_github("SamLar27/MIAnalysis", force = TRUE)

install.packages("doParallel")
devtools::install_github("SamLar27/MIAnalysis", force = TRUE)

library(MIAnalysis)

#VALIDATE : "MIAnalysis" appears in the Packages (section in R right menu) ?
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 2.	Defining the path
#--------------------------------------------------------------------------------------------------------------------------------------
# 2.1. Loading path
#Load_path <- file.path(
#"D:/Dossier CRCHUS/ORACLE - Validation Code/Validation/PlaceboAsthma_Data_and_Script_clean/For_MainAnalysis_Script/Main_FCS_single_pmm_logreg_polyreg_pmm"
#)

Load_path <- file.path(
  "Q:/Simon_Couillard/ORACLE - Validation Code/Validation/PlaceboAsthma_Data_and_Script_clean/For_MainAnalysis_Script/Main_FCS_single_pmm_logreg_polyreg_pmm"
)

#--------------------------------------------------------------------------------------------------------------------------------------
# 2.2. Saving data path
Saving_path <- file.path(
 "Q:/Simon_Couillard/ORACLE - Validation Code/Validation/PlaceboAsthma_Data_and_Script_clean/New_Dataset"
)

#--------------------------------------------------------------------------------------------------------------------------------------
# 2.3. Saving figure path
Figure_path <- file.path(
  "Q:/Simon_Couillard/ORACLE - Validation Code/Validation/PlaceboAsthma_Data_and_Script_clean/Figures_Final"
)
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 3.	Importing the dataset 

#--------------------------------------------------------------------------------------------------------------------------------------
# 3.1.	Read the dataframe
Load_file <- "ORACLE_MI_FULL_long_imp0_to_m.rds"
full_path  <- file.path(Load_path, Load_file)
file.exists(full_path)

data_long <- readRDS(full_path)

#--------------------------------------------------------------------------------------------------------------------------------------
# 3.2.	Convert to data.frame
ORACLE_after_imputation <- as.data.frame(data_long)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 4.	Variables formatting

## 4.1.	Creating of categorical variable

# 4.1.1.	Age_by_group
ORACLE_after_imputation$Age_by_group <- cut(
  ORACLE_after_imputation$Age,
  breaks = c(-10, 40, 50, 60, 100000),
  labels = c("<40", "≥40 to 50", "≥50 to <60", "≥60")
)

# 4.1.2.	BMI_by_group
ORACLE_after_imputation$BMI_by_group <- cut(
  ORACLE_after_imputation$BMI,
  breaks = c(-10, 25, 30, 35, 100000),
  labels = c("<25", "≥25 to 30", "≥30 to <35", "≥35")
)

# 4.1.3.	BEC_by_group
ORACLE_after_imputation$BEC_by_group <- cut(
  ORACLE_after_imputation$BEC,
  breaks = c(0, 0.15, 0.3, Inf),
  labels = c("<1.5", "≥1.5 to <0.3", "≥0.3"),
  right  = FALSE
)

# 4.1.4.	FeNO_by_group
ORACLE_after_imputation$FeNO_by_group <- cut(
  ORACLE_after_imputation$FeNO,
  breaks = c(0, 25, 50, Inf),
  labels = c("<25", "≥25 to <50", "≥50"),
  right  = TRUE
)

# 4.1.5.	IgE_by_group
ORACLE_after_imputation$IgE_by_group <- cut(
  ORACLE_after_imputation$IgE,
  breaks = c(0, 150, 600, 100000),
  labels = c("<150", "≥150 to <600", "≥600")
)

# 4.1.6.	ACQ_score_0W_by_group
ORACLE_after_imputation$ACQ_score_0W_by_group <- cut(
  ORACLE_after_imputation$ACQ_score_0W,
  breaks = c(-10, 1.5, 3, 100000),
  labels = c("<1.5", "≥1.5 to <3", "≥3")
)

# 4.1.7.	Attack_history_last12mo_yes_no
ORACLE_after_imputation <- ORACLE_after_imputation |>
  dplyr::mutate(
    Attack_history_last12mo_yes_no = factor(
      Attack_history_last12mo_yes_no,
      levels = c("No", "Yes")
    )
  )

# 4.1.8.	Attack_history_cat
ORACLE_after_imputation <- ORACLE_after_imputation |>
  dplyr::mutate(
    Attack_history_cat = dplyr::case_when(
      Attack_history_num == 0   ~ "0",
      Attack_history_num == 1   ~ "1",
      Attack_history_num >= 2   ~ "≥2",
      is.na(Attack_history_num) ~ NA_character_,
      TRUE                      ~ NA_character_
    ),
    Attack_history_cat = factor(Attack_history_cat, levels = c("0", "1", "≥2"))
  )

# 4.1.9.	Attack_number_during_followup_yes_no
ORACLE_after_imputation <- ORACLE_after_imputation |>
  dplyr::mutate(
    Attack_number_during_followup_yes_no = dplyr::case_when(
      Attack_number_during_followup == 0   ~ "No",
      Attack_number_during_followup >= 1   ~ "Yes",
      is.na(Attack_number_during_followup) ~ NA_character_,
      TRUE                                 ~ NA_character_
    ),
    Attack_number_during_followup_yes_no = factor(
      Attack_number_during_followup_yes_no,
      levels = c("No", "Yes")
    )
  )

# 4.1.10.	Attack_number_during_followup_0_1_2
ORACLE_after_imputation <- ORACLE_after_imputation |>
  dplyr::mutate(
    Attack_number_during_followup_0_1_2 = dplyr::case_when(
      Attack_number_during_followup == 0  ~ "0",
      Attack_number_during_followup == 1  ~ "1",
      Attack_number_during_followup >= 2  ~ "\u22652",
      is.na(Attack_number_during_followup) ~ NA_character_
    ),
    Attack_number_during_followup_0_1_2 = factor(
      Attack_number_during_followup_0_1_2,
      levels = c("0", "1", "\u22652")
    )
  )

## 4.2.	Creating of variable combining others

# 4.2.1.	LTRA_or_LAMA_or_Theophylline
ORACLE_after_imputation$LTRA_or_LAMA_or_Theophylline <- dplyr::case_when(
  ORACLE_after_imputation$LTRA          == "Yes" |
    ORACLE_after_imputation$LAMA        == "Yes" |
    ORACLE_after_imputation$Theophylline == "Yes" ~ "Yes",
  TRUE ~ "No"
) |>
  factor(levels = c("No", "Yes"))

# 4.2.2.	mOCS_dose_supraphysiologic
ORACLE_after_imputation <- ORACLE_after_imputation |>
  dplyr::mutate(
    mOCS_dose_supraphysiologic = dplyr::case_when(
      is.na(mOCS_dose) ~ NA_character_,
      mOCS_dose > 5    ~ "Yes",
      TRUE             ~ "No"
    ),
    mOCS_dose_supraphysiologic = factor(
      mOCS_dose_supraphysiologic,
      levels = c("No", "Yes")
    )
  )

# 4.2.3.	GINA_step
ORACLE_after_imputation <- ORACLE_after_imputation |>
  dplyr::mutate(
    GINA_step = dplyr::case_when(
      # Step 1: No ICS and No LTRA
      ICS_Dose_Category == "No" & LTRA == "No"                                                ~ "Step 1",
      # Step 2: Low ICS + No LTRA + No LABA OR No ICS + LTRA + No LABA
      (ICS_Dose_Category == "Low"  & LABA == "No") |
        (ICS_Dose_Category == "No" & LTRA == "Yes" & LABA == "No")                           ~ "Step 2",
      # Step 3: Low ICS + Additional controller
      ICS_Dose_Category == "Low" & (LABA == "Yes" | LTRA_or_LAMA_or_Theophylline == "Yes")  ~ "Step 3",
      # Step 4: Medium ICS
      ICS_Dose_Category == "Medium"                                                           ~ "Step 4",
      # Step 5: High ICS or supraphysiologic mOCS
      ICS_Dose_Category == "High" | mOCS_dose_supraphysiologic == "Yes"                      ~ "Step 5",
      TRUE                                                                                    ~ NA_character_
    ),
    GINA_step = factor(
      GINA_step,
      levels  = c("Step 1", "Step 2", "Step 3", "Step 4", "Step 5"),
      ordered = TRUE
    )
  )

# 4.2.4.	GINA_step_numeric
ORACLE_after_imputation$GINA_step_numeric <- as.numeric(ORACLE_after_imputation$GINA_step)
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 5.	Calculation of the variable by definite change

## 5.1.	Per-unit variables
ORACLE_after_imputation <- ORACLE_after_imputation %>%
  # Age per 10 year increase
  dplyr::mutate(Age_per_10 = Age / 10) %>%
  # BMI per 5 increase
  dplyr::mutate(BMI_per_5 = BMI / 5) %>%
  # --- Pre-BD ---
  # FEV1 pre-BD per 10% decrease
  dplyr::mutate(FEV1_preBD_PCT_0W_per_10 = -FEV1_preBD_PCT_0W / 10) %>%
  # FEV1 pre-BD L per 0.1 decrease in absolute units
  dplyr::mutate(FEV1_preBD_L_0W_per_10 = -FEV1_preBD_L_0W / 0.1) %>%
  # FEV1 reversibility 10% decrease
  dplyr::mutate(FEV1_reversibility_0W_per_10 = FEV1_reversibility_0W / 10) %>%
  # FVC pre-BD per 10% decrease
  dplyr::mutate(FVC_preBD_PCT_0W_per_10 = -FVC_preBD_PCT_0W / 10) %>%
  # Tif (%) pre-BD per 10% decrease
  dplyr::mutate(Tif_preBD_PCT_0W_per_10 = -Tif_preBD_PCT_0W / 10) %>%
  # Tif (raw) pre-BD per 0.1 decrease in absolute units
  dplyr::mutate(Tif_preBD_0W_per_10 = -Tif_preBD_0W / 0.1) %>%
  # --- Post-BD ---
  # FEV1 post-BD per 10% decrease
  dplyr::mutate(FEV1_postBD_PCT_0W_per_10 = -FEV1_postBD_PCT_0W / 10) %>%
  # FVC post-BD per 10% decrease
  dplyr::mutate(FVC_postBD_PCT_0W_per_10 = -FVC_postBD_PCT_0W / 10) %>%
  # Tif (%) post-BD per 10% decrease
  dplyr::mutate(Tif_postBD_PCT_0W_per_10 = -Tif_postBD_PCT_0W / 10) %>%
  # Tif (raw) post-BD per 0.1 decrease in absolute units
  dplyr::mutate(Tif_postBD_0W_per_10 = -Tif_postBD_0W / 0.1)
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# Creation of the variable of ouctome placebo change

ORACLE_after_imputation <- ORACLE_after_imputation %>%
  ## Placebo_response_Attack:
  dplyr::mutate(Placebo_response_Attack = Attack_history_num - Attack_number_during_followup) %>%
  ## Placebo_response_Lung:
  dplyr::mutate(Placebo_response_Lung = (FEV1_preBD_L_52W - FEV1_preBD_L_0W) * 1000) %>%
  ## ACQ_score_26_25_24: follow-up ACQ (26W preferred, then 25W, then 24W)
  dplyr::mutate(ACQ_score_26_25_24 = dplyr::coalesce(ACQ_score_26W, ACQ_score_25W, ACQ_score_24W)) %>%
  ## Placebo_response_ACQ: baseline minus follow-up ACQ
  dplyr::mutate(Placebo_response_ACQ = ACQ_score_0W - ACQ_score_26_25_24)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
## Through Flowchart for the 3 analyses

# 6. Selection of the subjects to include in the analysis
#    All eligibility is determined on the non-imputed data (.imp == 0),
#    then applied across all imputations.
#    Common criteria: Trial_Design == "blinded" & Trial_control_type == "Placebo"

imp0 <- ORACLE_after_imputation |> dplyr::filter(.imp == 0)

# --- Counts at each selection step (non-imputed data) ---
total_n    <- nrow(imp0)
total_rcts <- dplyr::n_distinct(imp0$Enrolled_Trial_name)

imp0_placebo <- imp0 |>
  dplyr::filter(
    Trial_Design       == "blinded",
    Trial_control_type == "Placebo"
  )
n_placebo      <- nrow(imp0_placebo)
rcts_placebo   <- dplyr::n_distinct(imp0_placebo$Enrolled_Trial_name)
n_excl_placebo <- total_n - n_placebo

#--------------------------------------------------------------------------------------------------------------------------------------
# 6.1  ORACLE_placebo_attack
#      Requires Attack_history_num and Attack_number_during_followup both observed
#      Also excludes LUTE, MILLY, VERSE (Attack_history_num singly imputed)
unique(imp0_placebo$Enrolled_Trial_name)

trials_imputed_attack <- c("LUTE", "MILLY", "VERSE","COSTA")

ids_attack <- imp0_placebo |>
  dplyr::filter(
    !is.na(Attack_history_num),
    !is.na(Attack_number_during_followup),
    !Enrolled_Trial_name %in% trials_imputed_attack
  ) |>
  dplyr::pull(Sequential_number)

n_attack      <- length(ids_attack)
rcts_attack   <- dplyr::n_distinct(
  imp0_placebo$Enrolled_Trial_name[imp0_placebo$Sequential_number %in% ids_attack]
)
n_excl_attack <- n_placebo - n_attack

ORACLE_placebo_attack <- ORACLE_after_imputation |>
  dplyr::filter(
    Trial_Design       == "blinded",
    Trial_control_type == "Placebo",
    Sequential_number  %in% ids_attack
  )

ORACLE_Placebo_Attack_Raw <- ORACLE_placebo_attack |> dplyr::filter(.imp == 0)%>% 
  mutate(
    Groupe = 'Asthma attacks'
  )

ORACLE_Placebo_Attack_imputed <- ORACLE_placebo_attack |>
  dplyr::filter(.imp %in% 1:10)%>% 
  group_by(.imp) %>% 
  mutate(
    BEC_IQR = quantile(BEC, 0.75) - quantile(BEC, 0.25),
    FeNO_IQR = quantile(FeNO, 0.75) - quantile(FeNO, 0.25),
    BEC_per_IQR = BEC / BEC_IQR,
    FeNO_per_IQR = FeNO / FeNO_IQR,
  ) %>% 
  ungroup()

#--------------------------------------------------------------------------------------------------------------------------------------
# 6.2  ORACLE_placebo_lung
#      Requires FEV1_preBD_L_0W and FEV1_preBD_L_52W both observed

ids_lung <- imp0_placebo |>
  dplyr::filter(
    !is.na(FEV1_preBD_L_0W),
    !is.na(FEV1_preBD_L_52W)
  ) |>
  dplyr::pull(Sequential_number)

n_lung      <- length(ids_lung)
rcts_lung   <- dplyr::n_distinct(
  imp0_placebo$Enrolled_Trial_name[imp0_placebo$Sequential_number %in% ids_lung]
)
n_excl_lung <- n_placebo - n_lung

ORACLE_placebo_lung <- ORACLE_after_imputation |>
  dplyr::filter(
    Trial_Design       == "blinded",
    Trial_control_type == "Placebo",
    Sequential_number  %in% ids_lung
  )

ORACLE_Placebo_Lung_Raw <- ORACLE_placebo_lung |> dplyr::filter(.imp == 0)%>% 
  mutate(
    Groupe = 'Lung function'
  )

ORACLE_Placebo_Lung_imputed <- ORACLE_placebo_lung |>
  dplyr::filter(.imp %in% 1:10) %>% 
  group_by(.imp) %>% 
  mutate(
    BEC_IQR = quantile(BEC, 0.75) - quantile(BEC, 0.25),
    FeNO_IQR = quantile(FeNO, 0.75) - quantile(FeNO, 0.25),
    BEC_per_IQR = BEC / BEC_IQR,
    FeNO_per_IQR = FeNO / FeNO_IQR,
  ) %>% 
  ungroup()

#--------------------------------------------------------------------------------------------------------------------------------------
# 6.3  ORACLE_placebo_symptom
#      Requires ACQ_score_0W and at least one of ACQ_score_24W / _25W / _26W observed

ids_symptom <- imp0_placebo |>
  dplyr::filter(
    !is.na(ACQ_score_0W),
    (!is.na(ACQ_score_26W) | !is.na(ACQ_score_25W) | !is.na(ACQ_score_24W))
  ) |>
  dplyr::pull(Sequential_number)

n_symptom      <- length(ids_symptom)
rcts_symptom   <- dplyr::n_distinct(
  imp0_placebo$Enrolled_Trial_name[imp0_placebo$Sequential_number %in% ids_symptom]
)
n_excl_symptom <- n_placebo - n_symptom

ORACLE_placebo_symptom <- ORACLE_after_imputation |>
  dplyr::filter(
    Trial_Design       == "blinded",
    Trial_control_type == "Placebo",
    Sequential_number  %in% ids_symptom
  )

ORACLE_Placebo_ACQ_Raw <- ORACLE_placebo_symptom |> dplyr::filter(.imp == 0)%>% 
  mutate(
    Groupe = 'ACQ-5 score'
  )

ORACLE_Placebo_ACQ_imputed <- ORACLE_placebo_symptom |>
  dplyr::filter(.imp %in% 1:10) %>% 
  group_by(.imp) %>% 
  mutate(
    BEC_IQR = quantile(BEC, 0.75) - quantile(BEC, 0.25),
    FeNO_IQR = quantile(FeNO, 0.75) - quantile(FeNO, 0.25),
    BEC_per_IQR = BEC / BEC_IQR,
    FeNO_per_IQR = FeNO / FeNO_IQR,
  ) %>% 
  ungroup()

# TABLE 1

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 8.	Creation of Tables with p-values

Data_table1 <- rbind(ORACLE_Placebo_Attack_Raw, ORACLE_Placebo_Lung_Raw, ORACLE_Placebo_ACQ_Raw) %>% 
  dplyr::select(c(
    Groupe,
    
    ## Demographic
    Age,
    Sex,
    BMI,
    Smoking_history_yes_no,
    
    ## Comorbidities
    Airborne_allergen_sensitisation,
    Eczema,
    Allergic_Rhinitis,
    CRSsNP,
    CRSwNP,
    
    ## Baseline asthma
    GINA_step,
    mOCS,
    
    ## Asthma symptoms
    ACQ_score_0W,
    Attack_history_last12mo_yes_no,
    Attack_history_cat,
    Attack_number_during_followup_0_1_2,
    
    ## Baseline lung function
    FEV1_preBD_L_0W,
    FEV1_preBD_PCT_0W,
    Tif_preBD_0W,
    FEV1_reversibility_0W,
    
    ## Biomarkers
    BEC,
    FeNO,
    IgE
    
  )) %>%
  mutate(
    Groupe = factor(Groupe, levels = c("Asthma attacks", "Lung function", "ACQ-5 score"))
  )

# ============================================================
# P-value formatting function
# ============================================================

my_pvalue <- function(x) {
  
  case_when(
    is.na(x) ~ NA_character_,
    x < 0.001 ~ "<0.001",
    .default = formatC(
      signif(x, 2),
      format = "fg",
      digits = 2
    )
  )
}

# ============================================================
# Create Table 1
# ============================================================

Table1 <- Data_table1 %>% 
  
  tbl_summary(
    
    # --------------------------------------------------------
    # Grouping variable
    # --------------------------------------------------------
    
    by = Groupe,
    
    
    # --------------------------------------------------------
    # Summary statistics
    # --------------------------------------------------------
    
    statistic = list(
      
      # Continuous variables:
      # Median [25th percentile - 75th percentile]
      
      all_continuous() ~ "{median} [{p25} - {p75}]",
      
      
      # Categorical variables:
      # n (%)
      
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    
    # --------------------------------------------------------
    # Missing data
    # --------------------------------------------------------
    
    missing = "ifany",
    
    missing_text = "Missing",
    
    missing_stat = "{N_miss} ({p_miss}%)",
    
    
    # --------------------------------------------------------
    # Number of decimal places
    # --------------------------------------------------------
    
    digits = list(
      BMI ~ 0,
      FEV1_preBD_PCT_0W ~ 0,
      BEC ~ 2,
      FeNO ~ 0,
      IgE ~ 0
    )
  ) %>%
  
  
  # ==========================================================
# Add p-values
# ==========================================================

add_p(
  
  test = list(
    
    # ------------------------------------------------------
    # Continuous variables
    # ------------------------------------------------------
    
    # One-way ANOVA
    all_continuous() ~ "oneway.test",
    
    
    # ------------------------------------------------------
    # Categorical variables
    # ------------------------------------------------------
    
    # Fisher's exact test
    all_categorical() ~ "fisher.test",
    
    
    # ------------------------------------------------------
    # Variables specifically analyzed with Kruskal-Wallis
    # ------------------------------------------------------
    
    BEC ~ "kruskal.test",
    FeNO ~ "kruskal.test",
    IgE ~ "kruskal.test"
  ),
  
  
  # --------------------------------------------------------
  # Test arguments
  # --------------------------------------------------------
  
  test.args = list(
    
    # One-way ANOVA assuming equal variances
    all_tests("oneway.test") ~ list(
      var.equal = TRUE
    ),
    
    
    # Fisher's exact test using Monte Carlo simulation
    #
    # This avoids FEXACT errors for large contingency tables
    # such as GINA_step, Sex, Allergic_Rhinitis, etc.
    
    all_tests("fisher.test") ~ list(
      simulate.p.value = TRUE,
      B = 10000
    )
  ),
  
  
  # --------------------------------------------------------
  # P-value formatting
  # --------------------------------------------------------
  
  pvalue_fun = my_pvalue
)

# ============================================================
# Display table and extract 
# ============================================================

Table1

Table1_word <- as_flex_table(Table1)

# Export to Word
save_as_docx(
  "Table 1" = Table1_word,
  path = "Table1withP.docx"
)

#FIGURE 2

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# PART A : Models for ∆ Attack
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

# SET UP THE CPU PROCESSING FOR THE ANALYSES
parallel::detectCores()

## In the function, parallelism in the process or not

#######################################################################################################################################

# 4.1 Simple model with only Attack_history_num as covariable

# Predictor → random-slope variable mapping
# (scaled/log vars use their raw counterpart; binary/categorical use themselves)
.predictors_simple <- list(
  list(pred = "Age_per_10",                      slope = "Age"),
  list(pred = "Sex",                             slope = "Sex"),
  list(pred = "BMI_per_5",                       slope = "BMI"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no"),
  list(pred = "Airborne_allergen_sensitisation", slope = "Airborne_allergen_sensitisation"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis"),
  list(pred = "Eczema",                          slope = "Eczema"),
  list(pred = "CRSwNP",                          slope = "CRSwNP"),
  list(pred = "CRSsNP",                          slope = "CRSsNP"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num"),  # no covariable for this one
  list(pred = "FEV1_preBD_L_0W_per_10",          slope = "FEV1_preBD_L_0W"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_per_IQR",                     slope = "BEC_per_IQR"),
  list(pred = "FeNO_per_IQR",                    slope = "FeNO_per_IQR"),
  list(pred = "IgE_log10",                       slope = "IgE_log10")
)

.tables_simple <- vector("list", length(.predictors_simple))

for (i in seq_along(.predictors_simple)) {
  pred  <- .predictors_simple[[i]]$pred
  slope <- .predictors_simple[[i]]$slope
  
  # When the predictor IS Attack_history_num, drop it from covariables
  cov       <- if (pred == "Attack_history_num") NULL else "Attack_history_num"
  cov_slope <- if (pred == "Attack_history_num") NULL else "Attack_history_num"
  
  .model <- IPD_one_stage(
    data                        = ORACLE_Placebo_Attack_imputed,
    outcome_var                 = "Placebo_response_Attack",
    predictor_vars              = pred,
    covariables                 = "Attack_history_num",
    imp_col                     = ".imp",
    followup_offset             = "Yes",
    followup_col                = "Follow_up_duration_days",
    random_intercept_var        = "Enrolled_Trial_name",
    stratified_intercept_var    = "Enrolled_Trial_name",
    predictor_vars_random_slope = slope,
    covariables_random_slope    = "Attack_history_num",
    model_type                  = "lm",
    model_performance           = FALSE,
    weighted_intercept          = FALSE,
    parallel                    = TRUE,
    n_cores                     = 2
  )
  
  # Keep only the row for the predictor of interest
  # startsWith() handles dichotomic variables whose term gains a level suffix (e.g. Sex -> Sex1)
  .tables_simple[[i]] <- .model$table[startsWith(.model$table$term, pred), ]
  rm(.model); gc()
}

Coefficients_Placebo_Attack_simple_model <- dplyr::bind_rows(.tables_simple)
rm(.predictors_simple, .tables_simple)
Coefficients_Placebo_Attack_simple_model

#######################################################################################################################################
# 4.2.  Forest plot – Simple model (∆ Asthma Attack Rate)
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Label and category mapping ────────────────────────────────────────────────
df_plot <- Coefficients_Placebo_Attack_simple_model |>
  mutate(
    label = case_match(
      term,
      "Age_per_10"                         ~ "Age (per 10 years)",
      "SexMale"                            ~ "Sex (Male vs. Female)",
      "BMI_per_5"                          ~ "BMI (per 5 kg/m\u00b2)",
      "Smoking_history_yes_noYes"          ~ "Smoking history",
      "Airborne_allergen_sensitisationYes" ~ "Airborne allergen sensitisation",
      "Allergic_RhinitisYes"              ~ "Allergic rhinitis",
      "EczemaYes"                          ~ "Eczema",
      "CRSwNPYes"                          ~ "CRS with nasal polyps",
      "CRSsNPYes"                          ~ "CRS without nasal polyps",
      "Psychiatric_diseaseYes"             ~ "Psychiatric disease",
      "GINA_step_numeric"                  ~ "GINA step",
      "ACQ_score_0W"                       ~ "ACQ score (baseline)",
      "Attack_history_num"                 ~ "Exacerbation history (count)",
      "FEV1_preBD_L_0W_per_10"             ~ "FEV1 pre-BD (L, per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"       ~ "FEV1 reversibility (per 10% increase)",
      "BEC_per_IQR"                        ~ "Blood eosinophils (per 1 IQR)",
      "FeNO_per_IQR"                       ~ "FeNO (per 1 IQR)",
      "IgE_log10"                          ~ "Total IgE (log10)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5","Smoking_history_yes_noYes")
      ~ "Demographic",
      c("Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline lung function",
      c("BEC_per_IQR", "FeNO_per_IQR", "IgE_log10")
      ~ "Inflammatory biomarkers"
    ),
    category    = factor(category, levels = c(
      "Demographic", "Comorbidities", "Asthma history",
      "Baseline lung function", "Inflammatory biomarkers"
    )),
    label       = factor(label, levels = rev(label)),
    p_label     = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p.value)),
    significant = p.value < 0.05
  )

# ── Build flat data frame: spacer + header rows inserted per category ─────────
.cats     <- levels(df_plot$category)
.all_rows <- list()

for (.cat in .cats) {
  if (.cat != .cats[1]) {
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = paste0(".spacer.", gsub(" ", "_", .cat)), type = "spacer",
      row_label = "", estimate = NA_real_, lo = NA_real_, hi = NA_real_,
      p_label = "", significant = FALSE
    )
  }
  .all_rows[[length(.all_rows) + 1]] <- tibble(
    y_id = paste0(".header.", gsub(" ", "_", .cat)), type = "header",
    row_label = .cat, estimate = NA_real_, lo = NA_real_, hi = NA_real_,
    p_label = "", significant = FALSE
  )
  .cat_vars <- df_plot |> filter(category == .cat)
  for (.k in seq_len(nrow(.cat_vars))) {
    .v <- .cat_vars[.k, ]
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = as.character(.v$label), type = "variable",
      row_label = as.character(.v$label),
      estimate = .v$estimate, lo = .v$`2.5 %`, hi = .v$`97.5 %`,
      p_label = .v$p_label, significant = .v$significant
    )
  }
}

.df_flat <- bind_rows(.all_rows) |>
  mutate(y_id = factor(y_id, levels = rev(y_id)))

.df_vars <- .df_flat |> filter(type == "variable")

# Bold headers via ggtext markdown, plain variables, blank spacers
.axis_labels <- setNames(
  case_when(
    .df_flat$type == "header"   ~ paste0("**", .df_flat$row_label, "**"),
    .df_flat$type == "variable" ~ .df_flat$row_label,
    TRUE                         ~ ""
  ),
  as.character(.df_flat$y_id)
)

# ── Forest plot (p-values placed outside panel via clip = "off") ──────────────
library(ggtext)
.x_pval <- 1.28   # data-unit x for p-value column (outside xlim of -1 to 1)
.n_rows  <- nrow(.df_flat)   # total factor levels (variables + headers + spacers)

.p_forest <- ggplot(.df_flat, aes(y = y_id)) +
  geom_blank(aes(x = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_errorbarh(
    data = .df_vars,
    aes(xmin = lo, xmax = hi, x = estimate),
    height = 0.3, linewidth = 0.55, color = "gray30"
  ) +
  geom_point(
    data = .df_vars,
    aes(x = estimate, fill = significant),
    shape = 22, size = 3, color = "gray20"
  ) +
  # P-values outside panel (right side)
  geom_text(
    data = .df_vars,
    aes(x = .x_pval, label = p_label),
    size = 3, hjust = 0.5
  ) +
  annotate("text",
           x = .x_pval, y = .n_rows + 1,
           label = "p-value", fontface = "bold", size = 3.2, hjust = 0.5
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "steelblue4"), guide = "none") +
  scale_x_continuous(breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) +
  scale_y_discrete(labels = .axis_labels) +
  coord_cartesian(xlim = c(-1, 1), clip = "off") +
  labs(
    x        = expression(Delta ~ "Asthma Attack Rate"),
    y        = NULL,
    title    = "Predictors of placebo response \u2014 Exacerbations",
    subtitle = "Simple models (each predictor adjusted for exacerbation history only)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y   = element_markdown(size = 9, hjust = 1, color = "black"),
    plot.margin   = margin(20, 85, 0, 5),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ── Directional arrow strip (below x-axis title) ──────────────────────────────
.p_arrows <- ggplot() +
  annotate("segment",
           x = -0.02, xend = -1, y = 0.7, yend = 0.7,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = -0.5, y = 0.15, label = "Reduced placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  annotate("segment",
           x = 0.02, xend = 1, y = 0.7, yend = 0.7,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = 0.5, y = 0.15, label = "Greater placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(0, 85, 5, 5))

# ── Combine and save ──────────────────────────────────────────────────────────
Figure_Placebo_Attack_simple_forest <- (.p_forest / .p_arrows) + plot_layout(heights = c(20, 1))

#######################################################################################################################################
#######################################################################################################################################
# FIGURE 2.A - In the manuscript
#######################################################################################################################################
#######################################################################################################################################
# 4.3 Predictors ∆Asthma attack in the final model
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Map every term back to its predictor / slope pair ────────────────────────
.predictor_map_final <- list(
  list(pred = "Age_per_10",                      slope = "Age",                             term = "Age_per_10"),
  list(pred = "Sex",                             slope = "Sex",                             term = "SexMale"),
  list(pred = "BMI_per_5",                       slope = "BMI",                             term = "BMI_per_5"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no",          term = "Smoking_history_yes_noYes"),
  list(pred = "Airborne_allergen_sensitisation", slope = "Airborne_allergen_sensitisation", term = "Airborne_allergen_sensitisationYes"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis",               term = "Allergic_RhinitisYes"),
  list(pred = "Eczema",                          slope = "Eczema",                          term = "EczemaYes"),
  list(pred = "CRSwNP",                          slope = "CRSwNP",                          term = "CRSwNPYes"),
  list(pred = "CRSsNP",                          slope = "CRSsNP",                          term = "CRSsNPYes"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease",             term = "Psychiatric_diseaseYes"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric",               term = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W",                    term = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num",              term = "Attack_history_num"),
  list(pred = "FEV1_preBD_L_0W_per_10",          slope = "FEV1_preBD_L_0W",                 term = "FEV1_preBD_L_0W_per_10"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10",    term = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_per_IQR",                     slope = "BEC_per_IQR",                     term = "BEC_per_IQR"),
  list(pred = "FeNO_per_IQR",                    slope = "FeNO_per_IQR",                    term = "FeNO_per_IQR"),
  list(pred = "IgE_log10",                       slope = "IgE_log10",                       term = "IgE_log10")
)

# ── Build covariate pool from predictors significant in the simple models ─────
# Each predictor is tested individually, adjusted for all significant predictors
# except itself (to avoid self-adjustment).
.sig_terms <- Coefficients_Placebo_Attack_simple_model |>
  filter(p.value < 0.05) |>
  pull(term)

.sig_entries    <- purrr::keep(.predictor_map_final, ~ .x$term %in% .sig_terms)
.sig_cov_preds  <- purrr::map_chr(.sig_entries, "pred")
.sig_cov_slopes <- purrr::map_chr(.sig_entries, "slope")

# ── Loop: test each predictor adjusted for significant covariables ─────────────
.tables_final <- vector("list", length(.predictor_map_final))
.cov_used     <- vector("list", length(.predictor_map_final))

for (i in seq_along(.predictor_map_final)) {
  pred  <- .predictor_map_final[[i]]$pred
  slope <- .predictor_map_final[[i]]$slope
  
  # Remove current predictor from covariate pool (no self-adjustment)
  .excl     <- which(.sig_cov_preds == pred)
  cov       <- if (length(.excl) > 0) .sig_cov_preds[-.excl]  else .sig_cov_preds
  #cov_slope <- if (length(.excl) > 0) .sig_cov_slopes[-.excl] else .sig_cov_slopes
  
  cov       <- if (length(cov) == 0) NULL else cov
  .cov_used[[i]] <- list(predictor = pred, covariates = cov) 
  #cov_slope <- if (length(cov_slope) == 0) NULL else cov_slope
  
  .model <- IPD_one_stage(
    data                        = ORACLE_Placebo_Attack_imputed,
    outcome_var                 = "Placebo_response_Attack",
    predictor_vars              = pred,
    covariables                 = cov,
    imp_col                     = ".imp",
    followup_offset             = "Yes",
    followup_col                = "Follow_up_duration_days",
    random_intercept_var        = "Enrolled_Trial_name",
    stratified_intercept_var    = "Enrolled_Trial_name",
    predictor_vars_random_slope = slope,
    covariables_random_slope    = "Attack_history_num",
    model_type                  = "lm",
    model_performance           = FALSE,
    weighted_intercept          = FALSE,
    parallel                    = FALSE,
    n_cores                     = 1
  )
  
  .tables_final[[i]] <- .model$table[startsWith(.model$table$term, pred), ]
  rm(.model); gc()
}

Coefficients_Placebo_Attack_final_model <- dplyr::bind_rows(.tables_final)
rm(.tables_final, .sig_entries, .predictor_map_final)
Coefficients_Placebo_Attack_final_model

# Check for the correct Covariates (in 'cov') - Should significant covariates from Model 1 (Simple model)
cov_table <- purrr::map_dfr(.cov_used, function(x) {
  tibble::tibble(
    predictor  = x$predictor,
    covariates = if (is.null(x$covariates)) "none" else paste(x$covariates, collapse = ", ")
  )
})

#######################################################################################################################################
# 4.4.  Forest plot – Final model (∆ Asthma Attack Rate)
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Label and category mapping ────────────────────────────────────────────────
df_plot <- Coefficients_Placebo_Attack_final_model |>
  mutate(
    label = case_match(
      term,
      "Age_per_10"                         ~ "Age (per 10 years)",
      "SexMale"                            ~ "Sex (Male vs. Female)",
      "BMI_per_5"                          ~ "BMI (per 5 kg/m\u00b2)",
      "Smoking_history_yes_noYes"          ~ "Smoking history",
      "Airborne_allergen_sensitisationYes" ~ "Airborne allergen sensitisation",
      "Allergic_RhinitisYes"              ~ "Allergic rhinitis",
      "EczemaYes"                          ~ "Eczema",
      "CRSwNPYes"                          ~ "CRS with nasal polyps",
      "CRSsNPYes"                          ~ "CRS without nasal polyps",
      "Psychiatric_diseaseYes"             ~ "Psychiatric disease",
      "GINA_step_numeric"                  ~ "GINA step",
      "ACQ_score_0W"                       ~ "ACQ score (baseline)",
      "Attack_history_num"                 ~ "Exacerbation history (count)",
      "FEV1_preBD_L_0W_per_10"             ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"       ~ "FEV1 reversibility (per 10% increase)",
      "BEC_per_IQR"                          ~ "Blood eosinophils (per 1 IQR)",
      "FeNO_per_IQR"                         ~ "FeNO (per 1 IQR)",
      "IgE_log10"                          ~ "Total IgE (log10)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5","Smoking_history_yes_noYes")
      ~ "Demographic",
      c("Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline lung function",
      c("BEC_per_IQR", "FeNO_per_IQR", "IgE_log10")
      ~ "Inflammatory biomarkers"
    ),
    category    = factor(category, levels = c(
      "Demographic", "Comorbidities", "Asthma history",
      "Baseline lung function", "Inflammatory biomarkers"
    )),
    # Drop unused category levels (categories absent from final model)
    category    = droplevels(category),
    label       = factor(label, levels = rev(label)),
    p_label     = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p.value)),
    rc_label    = sprintf("%.2f [%.2f; %.2f]", estimate, `2.5 %`, `97.5 %`),
    significant = p.value < 0.05
  )

# ── Build flat data frame: spacer + header rows inserted per category ─────────
.cats     <- levels(df_plot$category)
.all_rows <- list()

for (.cat in .cats) {
  if (.cat != .cats[1]) {
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = paste0(".spacer.", gsub(" ", "_", .cat)), type = "spacer",
      row_label = "", estimate = NA_real_, lo = NA_real_, hi = NA_real_,
      p_label = "", rc_label = "", significant = FALSE
    )
  }
  .all_rows[[length(.all_rows) + 1]] <- tibble(
    y_id = paste0(".header.", gsub(" ", "_", .cat)), type = "header",
    row_label = .cat, estimate = NA_real_, lo = NA_real_, hi = NA_real_,
    p_label = "", rc_label = "", significant = FALSE
  )
  .cat_vars <- df_plot |> filter(category == .cat)
  for (.k in seq_len(nrow(.cat_vars))) {
    .v <- .cat_vars[.k, ]
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = as.character(.v$label), type = "variable",
      row_label = as.character(.v$label),
      estimate = .v$estimate, lo = .v$`2.5 %`, hi = .v$`97.5 %`,
      p_label = .v$p_label, rc_label = .v$rc_label, significant = .v$significant
    )
  }
}

.df_flat <- bind_rows(.all_rows) |>
  mutate(y_id = factor(y_id, levels = rev(y_id)))

.df_vars <- .df_flat |> filter(type == "variable")

# Bold headers via ggtext markdown, plain variables, blank spacers
.axis_labels <- setNames(
  case_when(
    .df_flat$type == "header"   ~ paste0("**", .df_flat$row_label, "**"),
    .df_flat$type == "variable" ~ .df_flat$row_label,
    TRUE                         ~ ""
  ),
  as.character(.df_flat$y_id)
)

# ── Forest plot ───────────────────────────────────────────────────────────────
.x_pval <- 1.35
.n_rows  <- nrow(.df_flat)

.p_forest <- ggplot(.df_flat, aes(y = y_id)) +
  geom_blank(aes(x = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_errorbarh(
    data = .df_vars,
    aes(xmin = lo, xmax = hi, x = estimate),
    height = 0.3, linewidth = 0.55, color = "gray30"
  ) +
  geom_point(
    data = .df_vars,
    aes(x = estimate, fill = significant),
    shape = 22, size = 3, color = "gray20"
  ) +
  geom_text(
    data = .df_vars,
    aes(x = .x_pval, label = rc_label),
    size = 3, hjust = 0.5
  ) +
  annotate("text",
           x = .x_pval, y = .n_rows + 1,
           label = "aRC [CI95]", fontface = "bold", size = 3.2, hjust = 0.5
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "steelblue4"), guide = "none") +
  scale_x_continuous(breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) +
  scale_y_discrete(labels = .axis_labels) +
  coord_cartesian(xlim = c(-1, 1), clip = "off") +
  labs(
    x = expression(Delta ~ "Asthma Attack Rate"),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y   = element_markdown(size = 9, hjust = 1, color = "black"),
    plot.margin   = margin(20, 115, 0, 5),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ── Directional arrow strip (below x-axis title) ──────────────────────────────
.p_arrows <- ggplot() +
  annotate("segment",
           x = -0.02, xend = -1, y = 0.85, yend = 0.85,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = -0.5, y = 0.40, label = "Reduced placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  annotate("segment",
           x = 0.02, xend = 1, y = 0.85, yend = 0.85,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = 0.5, y = 0.40, label = "Greater placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(0, 85, 5, 5))

# ── Combine and save ──────────────────────────────────────────────────────────
Figure_Placebo_Attack_final_forest <- (.p_forest / .p_arrows) + plot_layout(heights = c(20, 1))
Figure_Placebo_Attack_final_forest

ggsave(
  filename = file.path(Figure_path, "Figure_Placebo_Attack_final_forest_IQR.pdf"),
  plot     = Figure_Placebo_Attack_final_forest,
  width    = 10, height = 6, units = "in"
)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# PART B : Models for ∆ FEV1
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

# 4.1 Simple model with only FEV1_preBD_L_0W_per_10 as covariable

# Predictor → random-slope variable mapping
# (scaled/log vars use their raw counterpart; binary/categorical use themselves)
.predictors_simple <- list(
  list(pred = "Age_per_10",                      slope = "Age"),
  list(pred = "Sex",                             slope = "Sex"),
  list(pred = "BMI_per_5",                       slope = "BMI"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no"),
  list(pred = "Airborne_allergen_sensitisation", slope = "Airborne_allergen_sensitisation"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis"),
  list(pred = "Eczema",                          slope = "Eczema"),
  list(pred = "CRSwNP",                          slope = "CRSwNP"),
  list(pred = "CRSsNP",                          slope = "CRSsNP"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num"),  # no covariable for this one
  list(pred = "FEV1_preBD_L_0W_per_10",          slope = "FEV1_preBD_L_0W"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_per_IQR",                     slope = "BEC_per_IQR"),
  list(pred = "FeNO_per_IQR",                    slope = "FeNO_per_IQR"),
  list(pred = "IgE_log10",                       slope = "IgE_log10")
)

.tables_simple <- vector("list", length(.predictors_simple))

for (i in seq_along(.predictors_simple)) {
  pred  <- .predictors_simple[[i]]$pred
  slope <- .predictors_simple[[i]]$slope
  
  # When the predictor IS FEV1_preBD_L_0W_per_10, drop it from covariables
  cov       <- if (pred == "FEV1_preBD_L_0W_per_10") NULL else "FEV1_preBD_L_0W_per_10"
 #cov_slope <- if (pred == "FEV1_preBD_L_0W_per_10") NULL else "FEV1_preBD_L_0W_per_10"
  
  .model <- IPD_one_stage(
    data                        = ORACLE_Placebo_Lung_imputed,
    outcome_var                 = "Placebo_response_Lung",
    predictor_vars              = pred,
    covariables                 = "FEV1_preBD_L_0W_per_10",
    imp_col                     = ".imp",
    #followup_offset = "Yes",
    #followup_col    = "Follow_up_duration_days",
    random_intercept_var        = "Enrolled_Trial_name",
    stratified_intercept_var    = "Enrolled_Trial_name",
    predictor_vars_random_slope = slope,
    covariables_random_slope    = "FEV1_preBD_L_0W_per_10",
    model_type                  = "lm",
    model_performance           = FALSE,
    weighted_intercept          = FALSE,
    parallel                    = TRUE,
    n_cores                     = 3
  )
  
  # Keep only the row for the predictor of interest
  # startsWith() handles dichotomic variables whose term gains a level suffix (e.g. Sex -> Sex1)
  .tables_simple[[i]] <- .model$table[startsWith(.model$table$term, pred), ]
  rm(.model); gc()
}

Coefficients_Placebo_Lung_simple_model <- dplyr::bind_rows(.tables_simple)
rm(.predictors_simple, .tables_simple)
Coefficients_Placebo_Lung_simple_model

#######################################################################################################################################
# 4.2.  Forest plot – Simple model (∆ lung)
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Label and category mapping ────────────────────────────────────────────────
df_plot <- Coefficients_Placebo_Lung_simple_model |>
  mutate(
    label = case_match(
      term,
      "Age_per_10"                         ~ "Age (per 10 years)",
      "SexMale"                            ~ "Sex (Male vs. Female)",
      "BMI_per_5"                          ~ "BMI (per 5 kg/m\u00b2)",
      "Smoking_history_yes_noYes"          ~ "Smoking history",
      "Airborne_allergen_sensitisationYes" ~ "Airborne allergen sensitisation",
      "Allergic_RhinitisYes"              ~ "Allergic rhinitis",
      "EczemaYes"                          ~ "Eczema",
      "CRSwNPYes"                          ~ "CRS with nasal polyps",
      "CRSsNPYes"                          ~ "CRS without nasal polyps",
      "Psychiatric_diseaseYes"             ~ "Psychiatric disease",
      "GINA_step_numeric"                  ~ "GINA step",
      "ACQ_score_0W"                       ~ "ACQ score (baseline)",
      "Attack_history_num"                 ~ "Exacerbation history (count)",
      "FEV1_preBD_L_0W_per_10"             ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"       ~ "FEV1 reversibility (per 10% increase)",
      "BEC_per_IQR"                        ~ "Blood eosinophils (per 1 IQR)",
      "FeNO_per_IQR"                       ~ "FeNO (per 1 IQR)",
      "IgE_log10"                          ~ "Total IgE (log10)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5","Smoking_history_yes_noYes")
      ~ "Demographic",
      c("Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline lung function",
      c("BEC_per_IQR", "FeNO_per_IQR", "IgE_log10")
      ~ "Inflammatory biomarkers"
    ),
    category    = factor(category, levels = c(
      "Demographic", "Comorbidities", "Asthma history",
      "Baseline lung function", "Inflammatory biomarkers"
    )),
    label       = factor(label, levels = rev(label)),
    p_label     = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p.value)),
    significant = p.value < 0.05
  )

# ── Build flat data frame: spacer + header rows inserted per category ─────────
.cats     <- levels(df_plot$category)
.all_rows <- list()

for (.cat in .cats) {
  if (.cat != .cats[1]) {
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = paste0(".spacer.", gsub(" ", "_", .cat)), type = "spacer",
      row_label = "", estimate = NA_real_, lo = NA_real_, hi = NA_real_,
      p_label = "", significant = FALSE
    )
  }
  .all_rows[[length(.all_rows) + 1]] <- tibble(
    y_id = paste0(".header.", gsub(" ", "_", .cat)), type = "header",
    row_label = .cat, estimate = NA_real_, lo = NA_real_, hi = NA_real_,
    p_label = "", significant = FALSE
  )
  .cat_vars <- df_plot |> filter(category == .cat)
  for (.k in seq_len(nrow(.cat_vars))) {
    .v <- .cat_vars[.k, ]
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = as.character(.v$label), type = "variable",
      row_label = as.character(.v$label),
      estimate = .v$estimate, lo = .v$`2.5 %`, hi = .v$`97.5 %`,
      p_label = .v$p_label, significant = .v$significant
    )
  }
}

.df_flat <- bind_rows(.all_rows) |>
  mutate(y_id = factor(y_id, levels = rev(y_id)))

.df_vars <- .df_flat |> filter(type == "variable")

# Bold headers via ggtext markdown, plain variables, blank spacers
.axis_labels <- setNames(
  case_when(
    .df_flat$type == "header"   ~ paste0("**", .df_flat$row_label, "**"),
    .df_flat$type == "variable" ~ .df_flat$row_label,
    TRUE                         ~ ""
  ),
  as.character(.df_flat$y_id)
)

# ── Forest plot (p-values placed outside panel via clip = "off") ──────────────
.x_pval <- 204   # data-unit x for p-value column (outside xlim of -160 to 160)
.n_rows  <- nrow(.df_flat)   # total factor levels (variables + headers + spacers)

.p_forest <- ggplot(.df_flat, aes(y = y_id)) +
  geom_blank(aes(x = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_errorbarh(
    data = .df_vars,
    aes(xmin = lo, xmax = hi, x = estimate),
    height = 0.3, linewidth = 0.55, color = "gray30"
  ) +
  geom_point(
    data = .df_vars,
    aes(x = estimate, fill = significant),
    shape = 22, size = 3, color = "gray20"
  ) +
  # P-values outside panel (right side)
  geom_text(
    data = .df_vars,
    aes(x = .x_pval, label = p_label),
    size = 3, hjust = 0.5
  ) +
  annotate("text",
           x = .x_pval, y = .n_rows + 1,
           label = "p-value", fontface = "bold", size = 3.2, hjust = 0.5
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "steelblue4"), guide = "none") +
  scale_x_continuous(breaks = seq(-160, 160, by = 40)) +
  scale_y_discrete(labels = .axis_labels) +
  coord_cartesian(xlim = c(-160, 160), clip = "off") +
  labs(
    x        = expression(Delta ~ "FEV (mL)"),
    y        = NULL,
    title    = "Predictors of placebo response \u2014 Lung function",
    subtitle = "Simple models (each predictor adjusted for FEV1% pre-BD)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y   = element_markdown(size = 9, hjust = 1, color = "black"),
    plot.margin   = margin(20, 85, 0, 5),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ── Directional arrow strip (below x-axis title) ──────────────────────────────
.p_arrows <- ggplot() +
  annotate("segment",
           x = -4, xend = -160, y = 0.7, yend = 0.7,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = -80, y = 0.15, label = "Reduced placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  annotate("segment",
           x = 4, xend = 160, y = 0.7, yend = 0.7,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = 80, y = 0.15, label = "Greater placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  scale_x_continuous(limits = c(-170, 170)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(0, 85, 5, 5))

# ── Combine and save ──────────────────────────────────────────────────────────
Figure_Placebo_Lung_simple_forest <- (.p_forest / .p_arrows) + plot_layout(heights = c(20, 1))

#######################################################################################################################################
#######################################################################################################################################
# FIGURE 2.B - In the manuscript
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 4.3 Predictors ∆Lung in the final model
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Map every term back to its predictor / slope pair ────────────────────────
.predictor_map_final <- list(
  list(pred = "Age_per_10",                      slope = "Age",                             term = "Age_per_10"),
  list(pred = "Sex",                             slope = "Sex",                             term = "SexMale"),
  list(pred = "BMI_per_5",                       slope = "BMI",                             term = "BMI_per_5"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no",          term = "Smoking_history_yes_noYes"),
  list(pred = "Airborne_allergen_sensitisation",  slope = "Airborne_allergen_sensitisation", term = "Airborne_allergen_sensitisationYes"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis",               term = "Allergic_RhinitisYes"),
  list(pred = "Eczema",                          slope = "Eczema",                          term = "EczemaYes"),
  list(pred = "CRSwNP",                          slope = "CRSwNP",                          term = "CRSwNPYes"),
  list(pred = "CRSsNP",                          slope = "CRSsNP",                          term = "CRSsNPYes"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease",             term = "Psychiatric_diseaseYes"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric",               term = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W",                    term = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num",              term = "Attack_history_num"),
  list(pred = "FEV1_preBD_L_0W_per_10",          slope = "FEV1_preBD_L_0W",                   term = "FEV1_preBD_L_0W_per_10"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10",    term = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_per_IQR",                       slope = "BEC_per_IQR",                   term = "BEC_per_IQR"),
  list(pred = "FeNO_per_IQR",                      slope = "FeNO_per_IQR",                  term = "FeNO_per_IQR"),
  list(pred = "IgE_log10",                       slope = "IgE_log10",                       term = "IgE_log10")
)

# ── Build covariate pool from predictors significant in the simple models ─────
# Each predictor is tested individually, adjusted for all significant predictors
# except itself (to avoid self-adjustment).
.sig_terms <- Coefficients_Placebo_Lung_simple_model |>
  filter(p.value < 0.05) |>
  pull(term)

.sig_entries    <- purrr::keep(.predictor_map_final, ~ .x$term %in% .sig_terms)
.sig_cov_preds  <- purrr::map_chr(.sig_entries, "pred")
.sig_cov_slopes <- purrr::map_chr(.sig_entries, "slope")

# ── Loop: test each predictor adjusted for significant covariables ─────────────
.tables_final <- vector("list", length(.predictor_map_final))
.cov_used     <- vector("list", length(.predictor_map_final))

for (i in seq_along(.predictor_map_final)) {
  pred  <- .predictor_map_final[[i]]$pred
  slope <- .predictor_map_final[[i]]$slope
  
  # Remove current predictor from covariate pool (no self-adjustment)
  .excl     <- which(.sig_cov_preds == pred)
  cov       <- if (length(.excl) > 0) .sig_cov_preds[-.excl]  else .sig_cov_preds
  #cov_slope <- if (length(.excl) > 0) .sig_cov_slopes[-.excl] else .sig_cov_slopes
  
  cov       <- if (length(cov) == 0) NULL else cov
  #cov_slope <- if (length(cov_slope) == 0) NULL else cov_slope
  .cov_used[[i]] <- list(predictor = pred, covariates = cov)
  
  .model <- IPD_one_stage(
    data                        = ORACLE_Placebo_Lung_imputed,,
    outcome_var                 = "Placebo_response_Lung",
    predictor_vars              = pred,
    covariables                 = cov,
    imp_col                     = ".imp",
    #followup_offset             = "Yes",
    #followup_col                = "Follow_up_duration_days",
    random_intercept_var        = "Enrolled_Trial_name",
    stratified_intercept_var    = "Enrolled_Trial_name",
    predictor_vars_random_slope = slope,
    covariables_random_slope    = "FEV1_preBD_L_0W_per_10",
    model_type                  = "lm",
    model_performance           = FALSE,
    weighted_intercept          = FALSE,
    parallel                    = TRUE,
    n_cores                     = 2
  )
  
  .tables_final[[i]] <- .model$table[startsWith(.model$table$term, pred), ]
  rm(.model); gc()
}

Coefficients_Placebo_Lung_final_model <- dplyr::bind_rows(.tables_final)
rm(.tables_final, .sig_entries, .predictor_map_final)
Coefficients_Placebo_Lung_final_model


cov_table <- purrr::map_dfr(.cov_used, function(x) {
  tibble::tibble(
    predictor  = x$predictor,
    covariates = if (is.null(x$covariates)) "none" else paste(x$covariates, collapse = ", ")
  )
})

#######################################################################################################################################
# 4.4.  Forest plot – Final model (∆ Lung)
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Label and category mapping ────────────────────────────────────────────────
df_plot <- Coefficients_Placebo_Lung_final_model |>
  mutate(
    label = case_match(
      term,
      "Age_per_10"                         ~ "Age (per 10 years)",
      "SexMale"                            ~ "Sex (Male vs. Female)",
      "BMI_per_5"                          ~ "BMI (per 5 kg/m\u00b2)",
      "Smoking_history_yes_noYes"          ~ "Smoking history",
      "Airborne_allergen_sensitisationYes" ~ "Airborne allergen sensitisation",
      "Allergic_RhinitisYes"              ~ "Allergic rhinitis",
      "EczemaYes"                          ~ "Eczema",
      "CRSwNPYes"                          ~ "CRS with nasal polyps",
      "CRSsNPYes"                          ~ "CRS without nasal polyps",
      "Psychiatric_diseaseYes"             ~ "Psychiatric disease",
      "GINA_step_numeric"                  ~ "GINA step",
      "ACQ_score_0W"                       ~ "ACQ score (baseline)",
      "Attack_history_num"                 ~ "Exacerbation history (count)",
      "FEV1_preBD_L_0W_per_10"             ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"       ~ "FEV1 reversibility (per 10% increase)",
      "BEC_per_IQR"                        ~ "Blood eosinophils (per 1 IQR)",
      "FeNO_per_IQR"                       ~ "FeNO (per 1 IQR)",
      "IgE_log10"                          ~ "Total IgE (log10)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5","Smoking_history_yes_noYes")
      ~ "Demographic",
      c("Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline lung function",
      c("BEC_per_IQR", "FeNO_per_IQR", "IgE_log10")
      ~ "Inflammatory biomarkers"
    ),
    category    = factor(category, levels = c(
      "Demographic", "Comorbidities", "Asthma history",
      "Baseline lung function", "Inflammatory biomarkers"
    )),
    # Drop unused category levels (categories absent from final model)
    category    = droplevels(category),
    label       = factor(label, levels = rev(label)),
    p_label     = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p.value)),
    rc_label    = sprintf("%.2f [%.2f; %.2f]", estimate, `2.5 %`, `97.5 %`),
    significant = p.value < 0.05
  )

# ── Build flat data frame: spacer + header rows inserted per category ─────────
.cats     <- levels(df_plot$category)
.all_rows <- list()

for (.cat in .cats) {
  if (.cat != .cats[1]) {
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = paste0(".spacer.", gsub(" ", "_", .cat)), type = "spacer",
      row_label = "", estimate = NA_real_, lo = NA_real_, hi = NA_real_,
      p_label = "", rc_label = "", significant = FALSE
    )
  }
  .all_rows[[length(.all_rows) + 1]] <- tibble(
    y_id = paste0(".header.", gsub(" ", "_", .cat)), type = "header",
    row_label = .cat, estimate = NA_real_, lo = NA_real_, hi = NA_real_,
    p_label = "", rc_label = "", significant = FALSE
  )
  .cat_vars <- df_plot |> filter(category == .cat)
  for (.k in seq_len(nrow(.cat_vars))) {
    .v <- .cat_vars[.k, ]
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = as.character(.v$label), type = "variable",
      row_label = as.character(.v$label),
      estimate = .v$estimate, lo = .v$`2.5 %`, hi = .v$`97.5 %`,
      p_label = .v$p_label, rc_label = .v$rc_label, significant = .v$significant
    )
  }
}

.df_flat <- bind_rows(.all_rows) |>
  mutate(y_id = factor(y_id, levels = rev(y_id)))

.df_vars <- .df_flat |> filter(type == "variable")

# Bold headers via ggtext markdown, plain variables, blank spacers
.axis_labels <- setNames(
  case_when(
    .df_flat$type == "header"   ~ paste0("**", .df_flat$row_label, "**"),
    .df_flat$type == "variable" ~ .df_flat$row_label,
    TRUE                         ~ ""
  ),
  as.character(.df_flat$y_id)
)

# ── Forest plot ───────────────────────────────────────────────────────────────
.x_pval <- 540   # data-unit x for p-value column (outside xlim of -400 to 400)
.n_rows  <- nrow(.df_flat)

.p_forest <- ggplot(.df_flat, aes(y = y_id)) +
  geom_blank(aes(x = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_errorbarh(
    data = .df_vars,
    aes(xmin = lo, xmax = hi, x = estimate),
    height = 0.3, linewidth = 0.55, color = "gray30"
  ) +
  geom_point(
    data = .df_vars,
    aes(x = estimate, fill = significant),
    shape = 22, size = 3, color = "gray20"
  ) +
  geom_text(
    data = .df_vars,
    aes(x = .x_pval, label = rc_label),
    size = 3, hjust = 0.5
  ) +
  annotate("text",
           x = .x_pval, y = .n_rows + 1,
           label = "aRC [CI95]", fontface = "bold", size = 3.2, hjust = 0.5
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "steelblue4"), guide = "none") +
  scale_x_continuous(breaks = seq(-400, 400, by = 100)) +
  scale_y_discrete(labels = .axis_labels) +
  coord_cartesian(xlim = c(-400, 400), clip = "off") +
  labs(
    x = expression(Delta ~ "FEV1 (mL)"),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y   = element_markdown(size = 9, hjust = 1, color = "black"),
    plot.margin   = margin(20, 115, 0, 5),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ── Directional arrow strip (below x-axis title) ──────────────────────────────
.p_arrows <- ggplot() +
  annotate("segment",
           x = -10, xend = -400, y = 0.85, yend = 0.85,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = -200, y = 0.40, label = "Reduced placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  annotate("segment",
           x = 10, xend = 400, y = 0.85, yend = 0.85,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = 200, y = 0.40, label = "Greater placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  scale_x_continuous(limits = c(-440, 440)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(0, 85, 5, 5))

# ── Combine and save ──────────────────────────────────────────────────────────
Figure_Placebo_Lung_final_forest <- (.p_forest / .p_arrows) + plot_layout(heights = c(20, 1))
Figure_Placebo_Lung_final_forest

ggsave(
  filename = file.path(Figure_path, "Figure_Placebo_Lung_final_forest_IQR.pdf"),
  plot     = Figure_Placebo_Lung_final_forest,
  width    = 10, height = 6, units = "in"
)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# PART C : Models for ∆ ACQ
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

# 4.	Models for ∆ ACQ

# 4.1 Simple model with only ACQ_score_0W as covariable

# Predictor → random-slope variable mapping
# (scaled/log vars use their raw counterpart; binary/categorical use themselves)
.predictors_simple <- list(
  list(pred = "Age_per_10",                      slope = "Age"),
  list(pred = "Sex",                             slope = "Sex"),
  list(pred = "BMI_per_5",                       slope = "BMI"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no"),
  list(pred = "Airborne_allergen_sensitisation", slope = "Airborne_allergen_sensitisation"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis"),
  list(pred = "Eczema",                          slope = "Eczema"),
  list(pred = "CRSwNP",                          slope = "CRSwNP"),
  list(pred = "CRSsNP",                          slope = "CRSsNP"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num"),  # no covariable for this one
  list(pred = "FEV1_preBD_L_0W_per_10",          slope = "FEV1_preBD_L_0W"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_per_IQR",                     slope = "BEC_per_IQR"),
  list(pred = "FeNO_per_IQR",                    slope = "FeNO_per_IQR"),
  list(pred = "IgE_log10",                       slope = "IgE_log10")
)

.tables_simple <- vector("list", length(.predictors_simple))

for (i in seq_along(.predictors_simple)) {
  pred  <- .predictors_simple[[i]]$pred
  slope <- .predictors_simple[[i]]$slope
  
  # When the predictor IS ACQ_score_0W, drop it from covariables
  cov       <- if (pred == "ACQ_score_0W") NULL else "ACQ_score_0W"
  #cov_slope <- if (pred == "ACQ_score_0W") NULL else "ACQ_score_0W"
  
  .model <- IPD_one_stage(
    data                        = ORACLE_Placebo_ACQ_imputed,
    outcome_var                 = "Placebo_response_ACQ",
    predictor_vars              = pred,
    covariables                 = "ACQ_score_0W",
    imp_col                     = ".imp",
    #followup_offset = "Yes",
    #followup_col    = "Follow_up_duration_days",
    random_intercept_var        = "Enrolled_Trial_name",
    stratified_intercept_var    = "Enrolled_Trial_name",
    predictor_vars_random_slope = slope,
    covariables_random_slope    = "ACQ_score_0W",
    model_type                  = "lm",
    model_performance           = FALSE,
    weighted_intercept          = FALSE,
    parallel                    = TRUE,
    n_cores                     = 2
  )
  
  # Keep only the row for the predictor of interest
  # startsWith() handles dichotomic variables whose term gains a level suffix (e.g. Sex -> Sex1)
  .tables_simple[[i]] <- .model$table[startsWith(.model$table$term, pred), ]
  rm(.model); gc()
}

Coefficients_Placebo_ACQ_simple_model <- dplyr::bind_rows(.tables_simple)
rm(.predictors_simple, .tables_simple)
Coefficients_Placebo_ACQ_simple_model

#######################################################################################################################################
# 4.2.  Forest plot – Simple model (∆ ACQ)
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Label and category mapping ────────────────────────────────────────────────
df_plot <- Coefficients_Placebo_ACQ_simple_model |>
  mutate(
    label = case_match(
      term,
      "Age_per_10"                         ~ "Age (per 10 years)",
      "SexMale"                            ~ "Sex (Male vs. Female)",
      "BMI_per_5"                          ~ "BMI (per 5 kg/m\u00b2)",
      "Smoking_history_yes_noYes"          ~ "Smoking history",
      "Airborne_allergen_sensitisationYes" ~ "Airborne allergen sensitisation",
      "Allergic_RhinitisYes"               ~ "Allergic rhinitis",
      "EczemaYes"                          ~ "Eczema",
      "CRSwNPYes"                          ~ "CRS with nasal polyps",
      "CRSsNPYes"                          ~ "CRS without nasal polyps",
      "Psychiatric_diseaseYes"             ~ "Psychiatric disease",
      "GINA_step_numeric"                  ~ "GINA step",
      "ACQ_score_0W"                       ~ "ACQ score (baseline)",
      "Attack_history_num"                 ~ "Exacerbation history (count)",
      "FEV1_preBD_L_0W_per_10"             ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"       ~ "FEV1 reversibility (per 10% increase)",
      "BEC_per_IQR"                        ~ "Blood eosinophils (per 1 IQR)",
      "FeNO_per_IQR"                       ~ "FeNO (per 1 IQR)",
      "IgE_log10"                          ~ "Total IgE (log10)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5","Smoking_history_yes_noYes")
      ~ "Demographic",
      c("Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline ACQ function",
      c("BEC_per_IQR", "FeNO_per_IQR", "IgE_log10")
      ~ "Inflammatory biomarkers"
    ),
    category    = factor(category, levels = c(
      "Demographic", "Comorbidities", "Asthma history",
      "Baseline ACQ function", "Inflammatory biomarkers"
    )),
    label       = factor(label, levels = rev(label)),
    p_label     = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p.value)),
    significant = p.value < 0.05
  )

# ── Build flat data frame: spacer + header rows inserted per category ─────────
.cats     <- levels(df_plot$category)
.all_rows <- list()

for (.cat in .cats) {
  if (.cat != .cats[1]) {
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = paste0(".spacer.", gsub(" ", "_", .cat)), type = "spacer",
      row_label = "", estimate = NA_real_, lo = NA_real_, hi = NA_real_,
      p_label = "", significant = FALSE
    )
  }
  .all_rows[[length(.all_rows) + 1]] <- tibble(
    y_id = paste0(".header.", gsub(" ", "_", .cat)), type = "header",
    row_label = .cat, estimate = NA_real_, lo = NA_real_, hi = NA_real_,
    p_label = "", significant = FALSE
  )
  .cat_vars <- df_plot |> filter(category == .cat)
  for (.k in seq_len(nrow(.cat_vars))) {
    .v <- .cat_vars[.k, ]
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = as.character(.v$label), type = "variable",
      row_label = as.character(.v$label),
      estimate = .v$estimate, lo = .v$`2.5 %`, hi = .v$`97.5 %`,
      p_label = .v$p_label, significant = .v$significant
    )
  }
}

.df_flat <- bind_rows(.all_rows) |>
  mutate(y_id = factor(y_id, levels = rev(y_id)))

.df_vars <- .df_flat |> filter(type == "variable")

# Bold headers via ggtext markdown, plain variables, blank spacers
.axis_labels <- setNames(
  case_when(
    .df_flat$type == "header"   ~ paste0("**", .df_flat$row_label, "**"),
    .df_flat$type == "variable" ~ .df_flat$row_label,
    TRUE                         ~ ""
  ),
  as.character(.df_flat$y_id)
)

# ── Forest plot (p-values placed outside panel via clip = "off") ──────────────
.x_pval <- 0.675   # data-unit x for p-value column (outside xlim of -0.5 to 0.5)
.n_rows  <- nrow(.df_flat)   # total factor levels (variables + headers + spacers)

.p_forest <- ggplot(.df_flat, aes(y = y_id)) +
  geom_blank(aes(x = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_errorbarh(
    data = .df_vars,
    aes(xmin = lo, xmax = hi, x = estimate),
    height = 0.3, linewidth = 0.55, color = "gray30"
  ) +
  geom_point(
    data = .df_vars,
    aes(x = estimate, fill = significant),
    shape = 22, size = 3, color = "gray20"
  ) +
  # P-values outside panel (right side)
  geom_text(
    data = .df_vars,
    aes(x = .x_pval, label = p_label),
    size = 3, hjust = 0.5
  ) +
  annotate("text",
           x = .x_pval, y = .n_rows + 1,
           label = "p-value", fontface = "bold", size = 3.2, hjust = 0.5
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "steelblue4"), guide = "none") +
  scale_x_continuous(breaks = seq(-0.5, 0.5, by = 0.25)) +
  scale_y_discrete(labels = .axis_labels) +
  coord_cartesian(xlim = c(-0.5, 0.5), clip = "off") +
  labs(
    x        = expression(Delta ~ "ACQ-5 Score"),
    y        = NULL,
    title    = "Predictors of placebo response \u2014 ACQ",
    subtitle = "Simple models (each predictor adjusted for ACQ score (baseline))"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y   = element_markdown(size = 9, hjust = 1, color = "black"),
    plot.margin   = margin(20, 85, 0, 5),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ── Directional arrow strip (below x-axis title) ──────────────────────────────
.p_arrows <- ggplot() +
  annotate("segment",
           x = -0.01, xend = -0.5, y = 0.7, yend = 0.7,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = -0.25, y = 0.15, label = "Reduced placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  annotate("segment",
           x = 0.01, xend = 0.5, y = 0.7, yend = 0.7,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = 0.25, y = 0.15, label = "Greater placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  scale_x_continuous(limits = c(-0.55, 0.55)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(0, 85, 5, 5))

# ── Combine and save ──────────────────────────────────────────────────────────
Figure_Placebo_ACQ_simple_forest <- (.p_forest / .p_arrows) + plot_layout(heights = c(20, 1))

#######################################################################################################################################
#######################################################################################################################################
# FIGURE 2.C - In the manuscript
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 4.3 Predictors ∆ACQ in the final model
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Map every term back to its predictor / slope pair ────────────────────────
.predictor_map_final <- list(
  list(pred = "Age_per_10",                      slope = "Age",                             term = "Age_per_10"),
  list(pred = "Sex",                             slope = "Sex",                             term = "SexMale"),
  list(pred = "BMI_per_5",                       slope = "BMI",                             term = "BMI_per_5"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no",          term = "Smoking_history_yes_noYes"),
  list(pred = "Airborne_allergen_sensitisation",  slope = "Airborne_allergen_sensitisation", term = "Airborne_allergen_sensitisationYes"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis",               term = "Allergic_RhinitisYes"),
  list(pred = "Eczema",                          slope = "Eczema",                          term = "EczemaYes"),
  list(pred = "CRSwNP",                          slope = "CRSwNP",                          term = "CRSwNPYes"),
  list(pred = "CRSsNP",                          slope = "CRSsNP",                          term = "CRSsNPYes"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease",             term = "Psychiatric_diseaseYes"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric",               term = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W",                    term = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num",              term = "Attack_history_num"),
  list(pred = "FEV1_preBD_L_0W_per_10",        slope = "FEV1_preBD_L_0W",               term = "FEV1_preBD_L_0W_per_10"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10",    term = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_per_IQR",                       slope = "BEC_per_IQR",                       term = "BEC_per_IQR"),
  list(pred = "FeNO_per_IQR",                      slope = "FeNO_per_IQR",                      term = "FeNO_per_IQR"),
  list(pred = "IgE_log10",                       slope = "IgE_log10",                       term = "IgE_log10")
)

# ── Build covariate pool from predictors significant in the simple models ─────
# Each predictor is tested individually, adjusted for all significant predictors
# except itself (to avoid self-adjustment).
.sig_terms <- Coefficients_Placebo_ACQ_simple_model |>
  filter(p.value < 0.05) |>
  pull(term)

.sig_entries    <- purrr::keep(.predictor_map_final, ~ .x$term %in% .sig_terms)
.sig_cov_preds  <- purrr::map_chr(.sig_entries, "pred")
.sig_cov_slopes <- purrr::map_chr(.sig_entries, "slope")

# ── Loop: test each predictor adjusted for significant covariables ─────────────
.tables_final <- vector("list", length(.predictor_map_final))
.cov_used     <- vector("list", length(.predictor_map_final))

for (i in seq_along(.predictor_map_final)) {
  pred  <- .predictor_map_final[[i]]$pred
  slope <- .predictor_map_final[[i]]$slope
  
  # Remove current predictor from covariate pool (no self-adjustment)
  .excl     <- which(.sig_cov_preds == pred)
  cov       <- if (length(.excl) > 0) .sig_cov_preds[-.excl]  else .sig_cov_preds
  #cov_slope <- if (length(.excl) > 0) .sig_cov_slopes[-.excl] else .sig_cov_slopes
  
  cov       <- if (length(cov) == 0) NULL else cov
  #cov_slope <- if (length(cov_slope) == 0) NULL else cov_slope
  .cov_used[[i]] <- list(predictor = pred, covariates = cov)
  
  .model <- IPD_one_stage(
    data                        = ORACLE_Placebo_ACQ_imputed,,
    outcome_var                 = "Placebo_response_ACQ",
    predictor_vars              = pred,
    covariables                 = cov,
    imp_col                     = ".imp",
    #followup_offset             = "Yes",
    #followup_col                = "Follow_up_duration_days",
    random_intercept_var        = "Enrolled_Trial_name",
    stratified_intercept_var    = "Enrolled_Trial_name",
    predictor_vars_random_slope = slope,
    covariables_random_slope    = "ACQ_score_0W",
    model_type                  = "lm",
    model_performance           = FALSE,
    weighted_intercept          = FALSE,
    parallel                    = TRUE,
    n_cores                     = 2
  )
  
  .tables_final[[i]] <- .model$table[startsWith(.model$table$term, pred), ]
  rm(.model); gc()
}

Coefficients_Placebo_ACQ_final_model <- dplyr::bind_rows(.tables_final)
rm(.tables_final, .sig_entries, .predictor_map_final)
Coefficients_Placebo_ACQ_final_model

cov_table <- purrr::map_dfr(.cov_used, function(x) {
  tibble::tibble(
    predictor  = x$predictor,
    covariates = if (is.null(x$covariates)) "none" else paste(x$covariates, collapse = ", ")
  )
})

#######################################################################################################################################
# 4.4.  Forest plot – Final model (∆ ACQ)
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Label and category mapping ────────────────────────────────────────────────
df_plot <- Coefficients_Placebo_ACQ_final_model |>
  mutate(
    label = case_match(
      term,
      "Age_per_10"                         ~ "Age (per 10 years)",
      "SexMale"                            ~ "Sex (Male vs. Female)",
      "BMI_per_5"                          ~ "BMI (per 5 kg/m\u00b2)",
      "Smoking_history_yes_noYes"          ~ "Smoking history",
      "Airborne_allergen_sensitisationYes" ~ "Airborne allergen sensitisation",
      "Allergic_RhinitisYes"               ~ "Allergic rhinitis",
      "EczemaYes"                          ~ "Eczema",
      "CRSwNPYes"                          ~ "CRS with nasal polyps",
      "CRSsNPYes"                          ~ "CRS without nasal polyps",
      "Psychiatric_diseaseYes"             ~ "Psychiatric disease",
      "GINA_step_numeric"                  ~ "GINA step",
      "ACQ_score_0W"                       ~ "ACQ score (baseline)",
      "Attack_history_num"                 ~ "Exacerbation history (count)",
      "FEV1_preBD_L_0W_per_10"             ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"       ~ "FEV1 reversibility (per 10% increase)",
      "BEC_per_IQR"                        ~ "Blood eosinophils (per 1 IQR)",
      "FeNO_per_IQR"                       ~ "FeNO (per 1 IQR)",
      "IgE_log10"                          ~ "Total IgE (log10)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5","Smoking_history_yes_noYes")
      ~ "Demographic",
      c("Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline ACQ function",
      c("BEC_per_IQR", "FeNO_per_IQR", "IgE_log10")
      ~ "Inflammatory biomarkers"
    ),
    category    = factor(category, levels = c(
      "Demographic", "Comorbidities", "Asthma history",
      "Baseline ACQ function", "Inflammatory biomarkers"
    )),
    # Drop unused category levels (categories absent from final model)
    category    = droplevels(category),
    label       = factor(label, levels = rev(label)),
    p_label     = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p.value)),
    rc_label    = sprintf("%.2f [%.2f; %.2f]", estimate, `2.5 %`, `97.5 %`),
    significant = p.value < 0.05
  )

# ── Build flat data frame: spacer + header rows inserted per category ─────────
.cats     <- levels(df_plot$category)
.all_rows <- list()

for (.cat in .cats) {
  if (.cat != .cats[1]) {
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = paste0(".spacer.", gsub(" ", "_", .cat)), type = "spacer",
      row_label = "", estimate = NA_real_, lo = NA_real_, hi = NA_real_,
      p_label = "", rc_label = "", significant = FALSE
    )
  }
  .all_rows[[length(.all_rows) + 1]] <- tibble(
    y_id = paste0(".header.", gsub(" ", "_", .cat)), type = "header",
    row_label = .cat, estimate = NA_real_, lo = NA_real_, hi = NA_real_,
    p_label = "", rc_label = "", significant = FALSE
  )
  .cat_vars <- df_plot |> filter(category == .cat)
  for (.k in seq_len(nrow(.cat_vars))) {
    .v <- .cat_vars[.k, ]
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = as.character(.v$label), type = "variable",
      row_label = as.character(.v$label),
      estimate = .v$estimate, lo = .v$`2.5 %`, hi = .v$`97.5 %`,
      p_label = .v$p_label, rc_label = .v$rc_label, significant = .v$significant
    )
  }
}

.df_flat <- bind_rows(.all_rows) |>
  mutate(y_id = factor(y_id, levels = rev(y_id)))

.df_vars <- .df_flat |> filter(type == "variable")

# Bold headers via ggtext markdown, plain variables, blank spacers
.axis_labels <- setNames(
  case_when(
    .df_flat$type == "header"   ~ paste0("**", .df_flat$row_label, "**"),
    .df_flat$type == "variable" ~ .df_flat$row_label,
    TRUE                         ~ ""
  ),
  as.character(.df_flat$y_id)
)

# ── Forest plot ───────────────────────────────────────────────────────────────
.x_pval <- 1.35   # data-unit x for p-value column (outside xlim of -1 to 1)
.n_rows  <- nrow(.df_flat)

.p_forest <- ggplot(.df_flat, aes(y = y_id)) +
  geom_blank(aes(x = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_errorbarh(
    data = .df_vars,
    aes(xmin = lo, xmax = hi, x = estimate),
    height = 0.3, linewidth = 0.55, color = "gray30"
  ) +
  geom_point(
    data = .df_vars,
    aes(x = estimate, fill = significant),
    shape = 22, size = 3, color = "gray20"
  ) +
  geom_text(
    data = .df_vars,
    aes(x = .x_pval, label = rc_label),
    size = 3, hjust = 0.5
  ) +
  annotate("text",
           x = .x_pval, y = .n_rows + 1,
           label = "aRC [CI95]", fontface = "bold", size = 3.2, hjust = 0.5
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "steelblue4"), guide = "none") +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.25)) +
  scale_y_discrete(labels = .axis_labels) +
  coord_cartesian(xlim = c(-1, 1), clip = "off") +
  labs(
    x = expression(Delta ~ "ACQ-5 Score"),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y   = element_markdown(size = 9, hjust = 1, color = "black"),
    plot.margin   = margin(20, 115, 0, 5),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ── Directional arrow strip (below x-axis title) ──────────────────────────────
.p_arrows <- ggplot() +
  annotate("segment",
           x = -0.02, xend = -1, y = 0.85, yend = 0.85,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = -0.5, y = 0.40, label = "Reduced placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  annotate("segment",
           x = 0.02, xend = 1, y = 0.85, yend = 0.85,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = 0.5, y = 0.40, label = "Greater placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  scale_x_continuous(limits = c(-1.1, 1.1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(0, 85, 5, 5))

# ── Combine and save ──────────────────────────────────────────────────────────
Figure_Placebo_ACQ_final_forest <- (.p_forest / .p_arrows) + plot_layout(heights = c(20, 1))
Figure_Placebo_ACQ_final_forest

ggsave(
  filename = file.path(Figure_path, "Figure_Placebo_ACQ_final_forest_IQR.pdf"),
  plot     = Figure_Placebo_ACQ_final_forest,
  width    = 10, height = 6, units = "in"
)

#######################################################################################################################################
#######################################################################################################################################
#     Exploratory analysis : Responder
#######################################################################################################################################
#######################################################################################################################################

# Prepare common dataset (Only include patients with all three outcomes)

ORACLE_all3 <- ORACLE_Placebo_Attack_Raw %>%
  dplyr::select(
    -Placebo_response_Lung,
    -Placebo_response_ACQ
  ) %>%
  dplyr::inner_join(
    ORACLE_Placebo_Lung_Raw %>%
      dplyr::select(
        Subject_ID,
        Placebo_response_Lung
      ),
    by = "Subject_ID"
  ) %>%
  dplyr::inner_join(
    ORACLE_Placebo_ACQ_Raw %>%
      dplyr::select(
        Subject_ID,
        Placebo_response_ACQ
      ),
    by = "Subject_ID"
  )

# Create responder variables 
# no asthma attacks during follow-up among patients with ≥120 days of follow-up; ≥100 mL improvement in FEV1; or ≥0.5-point improvement in ACQ-5.

ORACLE_all3 <- ORACLE_all3 %>%
  dplyr::filter(
    Follow_up_duration_days >= 120
  ) %>%
  dplyr::mutate(
    Responder_Attack = dplyr::case_when(
      Attack_number_during_followup == 0 ~ TRUE,
      Attack_number_during_followup > 0 ~ FALSE,
      TRUE ~ NA
    ),
    
    Responder_Lung = dplyr::case_when(
      Placebo_response_Lung >= 100 ~ TRUE,
      Placebo_response_Lung < 100 ~ FALSE,
      TRUE ~ NA
    ),
    
    Responder_ACQ = dplyr::case_when(
      Placebo_response_ACQ >= 0.5 ~ TRUE,
      Placebo_response_ACQ < 0.5 ~ FALSE,
      TRUE ~ NA
    )
  )

# ============================================================
# Table 1 — Characteristics of the common cohort (n = 1453)
# ============================================================

Table1_Placebo <- ORACLE_all3 %>%
  dplyr::mutate(
    GINA_step = droplevels(GINA_step)
  ) %>%
  gtsummary::tbl_summary(
    include = c(
      Age,
      Sex,
      BMI,
      Smoking_history_yes_no,
      Psychiatric_disease,
      Eczema,
      Allergic_Rhinitis,
      Airborne_allergen_sensitisation,
      CRSsNP,
      CRSwNP,
      GINA_step,
      mOCS,
      Attack_history_cat,
      ICU_or_ETI_history,
      ACQ_score_0W,
      FEV1_preBD_PCT_0W,
      FEV1_preBD_L_0W,
      Tif_preBD_0W,
      FEV1_reversibility_0W,
      BEC,
      FeNO,
      Follow_up_duration_days,
      Attack_number_during_followup_0_1_2,
      Placebo_response_Attack,
      FEV1_preBD_L_52W,
      Placebo_response_Lung,
      ACQ_score_26_25_24,
      Placebo_response_ACQ
    ),
    
    type = list(
      Sex ~ "dichotomous"
    ),
    
    value = list(
      Sex ~ "Female"
    ),
    
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    
    missing_text = "(Missing)",
    
    label = list(
      Age                             ~ "Age (y)",
      Sex                             ~ "Female sex",
      BMI                             ~ "BMI (kg/m²)",
      Smoking_history_yes_no          ~ "Former smokers",
      Psychiatric_disease             ~ "Psychiatric disease",
      Eczema                          ~ "Eczema",
      Allergic_Rhinitis               ~ "Allergic Rhinitis",
      Airborne_allergen_sensitisation ~ "Airborne allergen sensitisation",
      CRSsNP                          ~ "CRSsNP",
      CRSwNP                          ~ "CRSwNP",
      GINA_step                       ~ "GINA treatment step",
      mOCS                            ~ "mOCS",
      Attack_history_cat              ~ "Severe attack in the past 12 mo",
      ICU_or_ETI_history              ~ "History of ICU admission or intubation",
      ACQ_score_0W                    ~ "ACQ-5 score",
      FEV1_preBD_PCT_0W               ~ "FEV1 (% of predicted)",
      FEV1_preBD_L_0W                ~ "FEV1 (L)",
      Tif_preBD_0W                     ~ "FEV1/FVC",
      FEV1_reversibility_0W           ~ "FEV1 reversibility (%)",
      BEC                             ~ "Blood eosinophil count (×10⁹ cells per L)",
      FeNO                            ~ "FeNO (ppb)",
      Follow_up_duration_days         ~ "Follow-up duration (days)",
      Attack_number_during_followup_0_1_2 ~ "Number of severe attacks during follow-up",
      Placebo_response_Attack         ~ "Change in asthma attack rate by year",
      FEV1_preBD_L_52W                ~ "FEV1 (L) pre-BD at 52 weeks of follow-up",
      Placebo_response_Lung           ~ "FEV1 (mL) change at 52 weeks",
      ACQ_score_26_25_24              ~ "ACQ-5 at 24–26 weeks of follow-up",
      Placebo_response_ACQ            ~ "ACQ-5 score change at 24–26 weeks"
    )
  ) %>%
  
  # Remove missing rows for categorical/dichotomous variables
  gtsummary::modify_table_body(function(tbl) {
    tbl %>%
      dplyr::filter(
        !(row_type == "missing" &
            var_type %in% c("categorical", "dichotomous"))
      )
  }) %>%
  
  # Add section headers
  gtsummary::modify_table_body(function(tbl) {
    
    headers <- list(
      list(before = "Age",                             label = "Demographic"),
      list(before = "Psychiatric_disease",             label = "Comorbidities"),
      list(before = "GINA_step",                       label = "Medication"),
      list(before = "Attack_history_cat",              label = "Exacerbation history"),
      list(before = "ACQ_score_0W",                    label = "Symptoms"),
      list(before = "FEV1_preBD_L_0W",                 label = "Lung function pre-BD"),
      list(before = "BEC",                             label = "Inflammatory biomarkers"),
      list(before = "Follow_up_duration_days",         label = "Outcomes/follow-up data")
    )
    
    for (h in rev(headers)) {
      
      idx <- min(which(tbl$variable == h$before))
      
      new_row <- tibble::tibble(
        variable  = paste0("header_", gsub("[^a-zA-Z]", "_", h$label)),
        var_label = h$label,
        label     = h$label,
        row_type  = "label",
        var_type  = NA_character_
      )
      
      tbl <- dplyr::bind_rows(
        tbl[seq_len(idx - 1), ],
        new_row,
        tbl[idx:nrow(tbl), ]
      )
    }
    
    tbl %>%
      dplyr::mutate(
        label = dplyr::case_when(
          startsWith(variable, "header_") ~ label,
          row_type %in% c("level", "missing") ~
            paste0("\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0", label),
          TRUE ~
            paste0("\u00a0\u00a0\u00a0\u00a0", label)
        )
      )
  }) %>%
  
  gtsummary::modify_table_styling(
    columns = label,
    rows = startsWith(variable, "header_"),
    text_format = "bold"
  ) %>%
  
  gtsummary::modify_header(
    label ~ "",
    stat_0 ~ paste0(
      "**Common cohort (N = ",
      dplyr::n_distinct(ORACLE_all3$Subject_ID),
      ")**"
    )
  ) %>%
  
  gtsummary::modify_footnote(
    all_stat_cols() ~ NA
  )

Table1_Placebo


## Aucun critère vs au moins un 


ORACLE_responder <- ORACLE_all3 %>%
  dplyr::mutate(
    
    # Responder to at least one component
    Responder_any = dplyr::case_when(
      
      # At least one component meets the responder criterion
      Responder_Attack == TRUE |
        Responder_Lung == TRUE |
        Responder_ACQ == TRUE ~ TRUE,
      
      # All three available and none meets the criterion
      Responder_Attack == FALSE &
        Responder_Lung == FALSE &
        Responder_ACQ == FALSE ~ FALSE,
      
      # Otherwise insufficient information
      TRUE ~ NA
    ),
    
    Responder_group = factor(
      Responder_any,
      levels = c(FALSE, TRUE),
      labels = c(
        "Non-responder",
        "Responder to ≥1 component"
      )
    )
  ) %>%
  
  # Exclude patients without enough information to classify
  dplyr::filter(!is.na(Responder_group))




Table_Responder <- ORACLE_responder %>%
  
  dplyr::mutate(
    GINA_step = droplevels(GINA_step)
  ) %>%
  
  gtsummary::tbl_summary(
    
    # Same variables as Table 1
    by = Responder_group,
    
    include = c(
      Age,
      Sex,
      BMI,
      Smoking_history_yes_no,
      Psychiatric_disease,
      Eczema,
      Allergic_Rhinitis,
      Airborne_allergen_sensitisation,
      CRSsNP,
      CRSwNP,
      GINA_step,
      mOCS,
      Attack_history_cat,
      ICU_or_ETI_history,
      ACQ_score_0W,
      FEV1_preBD_PCT_0W,
      FEV1_preBD_L_0W,
      Tif_preBD_0W,
      FEV1_reversibility_0W,
      BEC,
      FeNO,
      Follow_up_duration_days,
      Attack_number_during_followup_0_1_2,
      Placebo_response_Attack,
      FEV1_preBD_L_52W,
      Placebo_response_Lung,
      ACQ_score_26_25_24,
      Placebo_response_ACQ
    ),
    
    type = list(
      Sex ~ "dichotomous"
    ),
    
    value = list(
      Sex ~ "Female"
    ),
    
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    
    missing_text = "(Missing)",
    
    label = list(
      Age                             ~ "Age (y)",
      Sex                             ~ "Female sex",
      BMI                             ~ "BMI (kg/m²)",
      Smoking_history_yes_no          ~ "Former smokers",
      Psychiatric_disease             ~ "Psychiatric disease",
      Eczema                          ~ "Eczema",
      Allergic_Rhinitis               ~ "Allergic Rhinitis",
      Airborne_allergen_sensitisation ~ "Airborne allergen sensitisation",
      CRSsNP                          ~ "CRSsNP",
      CRSwNP                          ~ "CRSwNP",
      GINA_step                       ~ "GINA treatment step",
      mOCS                            ~ "mOCS",
      Attack_history_cat              ~ "Severe attack in the past 12 mo",
      ICU_or_ETI_history              ~ "History of ICU admission or intubation",
      ACQ_score_0W                    ~ "ACQ-5 score",
      FEV1_preBD_PCT_0W               ~ "FEV1 (% of predicted)",
      FEV1_preBD_L_0W                 ~ "FEV1 (L)",
      Tif_preBD_0W                    ~ "FEV1/FVC",
      FEV1_reversibility_0W           ~ "FEV1 reversibility (%)",
      BEC                             ~ "Blood eosinophil count (×10⁹ cells per L)",
      FeNO                            ~ "FeNO (ppb)",
      Follow_up_duration_days         ~ "Follow-up duration (days)",
      Attack_number_during_followup_0_1_2 ~
        "Number of severe attacks during follow-up",
      Placebo_response_Attack         ~
        "Change in asthma attack rate by year",
      FEV1_preBD_L_52W                ~
        "FEV1 (L) pre-BD at 52 weeks of follow-up",
      Placebo_response_Lung           ~
        "FEV1 (mL) change at 52 weeks",
      ACQ_score_26_25_24              ~
        "ACQ-5 at 24–26 weeks of follow-up",
      Placebo_response_ACQ            ~
        "ACQ-5 score change at 24–26 weeks"
    )
  ) %>%
  
  # ----------------------------------------------------------
# P-values
# ----------------------------------------------------------

gtsummary::add_p(
  test = list(
    all_continuous() ~ "wilcox.test",
    all_categorical() ~ "chisq.test"
  ),
  
  # Use exact tests when expected cell counts are small
  test.args = list(
    all_categorical() ~ list(correct = FALSE)
  ),
  
  pvalue_fun = function(x) {
    gtsummary::style_pvalue(
      x,
      digits = 3
    )
  }
) %>%
  
  # ----------------------------------------------------------
# Remove missing rows for categorical/dichotomous variables
# ----------------------------------------------------------

gtsummary::modify_table_body(function(tbl) {
  tbl %>%
    dplyr::filter(
      !(row_type == "missing" &
          var_type %in% c("categorical", "dichotomous"))
    )
}) %>%
  
  # ----------------------------------------------------------
# Add section headers
# ----------------------------------------------------------

gtsummary::modify_table_body(function(tbl) {
  
  headers <- list(
    list(before = "Age",                             label = "Demographic"),
    list(before = "Psychiatric_disease",             label = "Comorbidities"),
    list(before = "GINA_step",                       label = "Medication"),
    list(before = "Attack_history_cat",              label = "Exacerbation history"),
    list(before = "ACQ_score_0W",                    label = "Symptoms"),
    list(before = "FEV1_preBD_L_0W",                 label = "Lung function pre-BD"),
    list(before = "BEC",                             label = "Inflammatory biomarkers"),
    list(before = "Follow_up_duration_days",         label = "Outcomes/follow-up data")
  )
  
  for (h in rev(headers)) {
    
    idx <- min(which(tbl$variable == h$before))
    
    new_row <- tibble::tibble(
      variable  = paste0(
        "header_",
        gsub("[^a-zA-Z]", "_", h$label)
      ),
      var_label = h$label,
      label     = h$label,
      row_type  = "label",
      var_type  = NA_character_
    )
    
    tbl <- dplyr::bind_rows(
      tbl[seq_len(idx - 1), ],
      new_row,
      tbl[idx:nrow(tbl), ]
    )
  }
  
  tbl %>%
    dplyr::mutate(
      label = dplyr::case_when(
        startsWith(variable, "header_") ~ label,
        
        row_type %in% c("level", "missing") ~
          paste0(
            "\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0",
            label
          ),
        
        TRUE ~
          paste0(
            "\u00a0\u00a0\u00a0\u00a0",
            label
          )
      )
    )
}) %>%
  
  # ----------------------------------------------------------
# Bold section headers
# ----------------------------------------------------------

gtsummary::modify_table_styling(
  columns = label,
  rows = startsWith(variable, "header_"),
  text_format = "bold"
) %>%
  
  # ----------------------------------------------------------
# Column headers
# ----------------------------------------------------------

gtsummary::modify_header(
  label ~ "",
  stat_1 ~ paste0(
    "**Non-responder**  \n",
    "*N = ", sum(ORACLE_responder$Responder_group == "Non-responder", na.rm = TRUE), "*"
  ),
  stat_2 ~ paste0(
    "**Responder to ≥1 component**  \n",
    "*N = ", sum(ORACLE_responder$Responder_group == "Responder to ≥1 component", na.rm = TRUE), "*"
  ),
  p.value ~ "**P-value**"
) %>%
  
  # Remove footnotes
  gtsummary::modify_footnote(
    all_stat_cols() ~ NA
  )


Table_Responder

Table_Responder %>%
  as_flex_table() %>%
  flextable::save_as_docx(
    path = "Table_Responder.docx"
  )

##### Responder analysis (Venn diagram)

ORACLE_all3 %>%
  dplyr::summarise(
    N = dplyr::n(),
    
    Attack_responders = sum(Responder_Attack),
    Lung_responders   = sum(Responder_Lung),
    ACQ_responders    = sum(Responder_ACQ),
    
    Attack_percent = mean(Responder_Attack) * 100,
    Lung_percent   = mean(Responder_Lung) * 100,
    ACQ_percent    = mean(Responder_ACQ) * 100
  )

ORACLE_all3 <- ORACLE_all3 %>%
  dplyr::mutate(
    Response_pattern = dplyr::case_when(
      Responder_Attack & !Responder_Lung & !Responder_ACQ ~ "Attack only",
      !Responder_Attack & Responder_Lung & !Responder_ACQ ~ "Lung only",
      !Responder_Attack & !Responder_Lung & Responder_ACQ ~ "ACQ only",
      Responder_Attack & Responder_Lung & !Responder_ACQ ~ "Attack + Lung",
      Responder_Attack & !Responder_Lung & Responder_ACQ ~ "Attack + ACQ",
      !Responder_Attack & Responder_Lung & Responder_ACQ ~ "Lung + ACQ",
      Responder_Attack & Responder_Lung & Responder_ACQ ~ "All 3",
      !Responder_Attack & !Responder_Lung & !Responder_ACQ ~ "None"
    )
  )

Response_pattern_summary <- ORACLE_all3 %>%
  dplyr::count(Response_pattern) %>%
  dplyr::mutate(
    Percent = 100 * n / sum(n)
  ) %>%
  dplyr::arrange(desc(n))

Response_pattern_summary

ORACLE_all3 <- ORACLE_all3 %>%
  dplyr::mutate(
    N_domains_response =
      Responder_Attack +
      Responder_Lung +
      Responder_ACQ
  )

ORACLE_all3 %>%
  dplyr::count(N_domains_response) %>%
  dplyr::mutate(
    Percent = 100 * n / sum(n)
  ) %>%
  dplyr::arrange(N_domains_response)

# ============================================================
# Venn diagram
# ============================================================

# ------------------------------------------------------------
# 1. Define the common cohort
# ------------------------------------------------------------

N_common <- nrow(ORACLE_all3)

# ------------------------------------------------------------
# 2. Define responder sets
# ------------------------------------------------------------

venn_data <- list(
  `No severe asthma attack` =
    ORACLE_all3$Subject_ID[
      ORACLE_all3$Responder_Attack
    ],
  
  `FEV1 change ≥ 0.1 L` =
    ORACLE_all3$Subject_ID[
      ORACLE_all3$Responder_Lung
    ],
  
  `ACQ-5 improvement ≥ 0.5` =
    ORACLE_all3$Subject_ID[
      ORACLE_all3$Responder_ACQ
    ]
)

# ------------------------------------------------------------
# 3. JAMA colour palette
# ------------------------------------------------------------

jama_colors <- ggsci::pal_jama()(3)

# ------------------------------------------------------------
# 4. Create Venn diagram
# ------------------------------------------------------------

Venn_plot <- ggvenn(
  venn_data,
  
  # JAMA colours
  fill_color = jama_colors,
  fill_alpha = 0.55,
  
  # Circle borders
  stroke_color = "black",
  stroke_alpha = 1,
  stroke_size = 0.8,
  
  # Set names
  set_name_color = "black",
  set_name_size = 5,
  
  # We add our own n (%) labels
  show_percentage = FALSE,
  show_elements = FALSE,
  text_size = 0,
  padding = 0
) +
  
  # ----------------------------------------------------------
# 5. Custom n (%) labels
# ----------------------------------------------------------

# Attack only
annotate(
  "text",
  x = -0.90,
  y = 0.75,
  label = "152\n(10.5%)",
  size = 4.5
) +
  
  # Lung only
  annotate(
    "text",
    x = 0.90,
    y = 0.75,
    label = "72\n(5.0%)",
    size = 4.5
  ) +
  
  # ACQ only
  annotate(
    "text",
    x = 0.015,
    y = -0.80,
    label = "171\n(11.8%)",
    size = 4.5
  ) +
  
  # Attack + Lung
  annotate(
    "text",
    x = 0,
    y = 0.75,
    label = "126\n(8.7%)",
    size = 4.5
  ) +
  
  # Attack + ACQ
  annotate(
    "text",
    x = -0.48,
    y = -0.07,
    label = "321\n(22.1%)",
    size = 4.5
  ) +
  
  # Lung + ACQ
  annotate(
    "text",
    x = 0.48,
    y = -0.07,
    label = "148\n(10.2%)",
    size = 4.5
  ) +
  
  # All 3
  annotate(
    "text",
    x = 0.005,
    y = 0.24,
    label = "326\n(22.4%)",
    size = 4.2,
    fontface = "bold"
  ) +
  
  # ----------------------------------------------------------
# 6. Title and subtitle
# ----------------------------------------------------------

labs(
  title = "Overlap of response across clinical domains",
  subtitle = paste0(
    "Common cohort (N = ", N_common, ")"
  ),
  caption = paste0(
     "Non-responders in all domains: 137 (9.4%)."
  )
) +
  
theme_void() +
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 4)
    ),
    
    plot.subtitle = element_text(
      size = 13,
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    
    plot.caption = element_text(
      size = 10,
      hjust = 0.5,
      margin = margin(t = 15)
    ),
    
    plot.margin = margin(
      t = 5,
      r = 0,
      b = 5,
      l = 0
    )
  )

# Display
Venn_plot

# Export 

ggsave(
  filename = "Venn_response_ORACLE.png",
  plot = Venn_plot,
  width = 9,
  height = 8,
  units = "in",
  dpi = 600
)


# Provide a table for values

Table_response <- ORACLE_all3 %>%
  dplyr::mutate(
    Response_pattern = dplyr::case_when(
      Responder_Attack & !Responder_Lung & !Responder_ACQ ~
        "Asthma attacks only",
      
      !Responder_Attack & Responder_Lung & !Responder_ACQ ~
        "Lung function only",
      
      !Responder_Attack & !Responder_Lung & Responder_ACQ ~
        "ACQ-5 score only",
      
      Responder_Attack & Responder_Lung & !Responder_ACQ ~
        "Asthma attacks + Lung function",
      
      Responder_Attack & !Responder_Lung & Responder_ACQ ~
        "Asthma attacks + ACQ-5 score",
      
      !Responder_Attack & Responder_Lung & Responder_ACQ ~
        "Lung function + ACQ-5 score",
      
      Responder_Attack & Responder_Lung & Responder_ACQ ~
        "All three domains",
      
      !Responder_Attack & !Responder_Lung & !Responder_ACQ ~
        "No response in any domain"
    )
  ) %>%
  dplyr::count(Response_pattern, name = "n") %>%
  dplyr::mutate(
    Percent = 100 * n / N_common
  ) %>%
  dplyr::mutate(
    Response_pattern = factor(
      Response_pattern,
      levels = c(
        "No response in any domain",
        "Asthma attacks only",
        "Lung function only",
        "ACQ-5 score only",
        "Asthma attacks + Lung function",
        "Asthma attacks + ACQ-5 score",
        "Lung function + ACQ-5 score",
        "All three domains"
      )
    )
  ) %>%
  dplyr::arrange(Response_pattern) %>%
  dplyr::mutate(
    `n (%)` = paste0(
      n,
      " (",
      sprintf("%.1f", Percent),
      "%)"
    )
  ) %>%
  dplyr::select(
    `Response pattern` = Response_pattern,
    `n (%)`
  )

Table_response


ft_response <- Table_response %>%
  flextable() %>%
  theme_booktabs() %>%
  bold(
    part = "header"
  ) %>%
  align(
    j = 2,
    align = "center",
    part = "all"
  ) %>%
  align(
    j = 1,
    align = "left",
    part = "all"
  ) %>%
  autofit() %>%
  set_caption(
    caption = paste0(
      "Table. Response patterns across clinical domains ",
      "(common cohort, N = ", N_common, ")"
    )
  )

# ------------------------------------------------------------
# Export to Word
# ------------------------------------------------------------

doc <- officer::read_docx() %>%
  flextable::body_add_flextable(
    value = ft_response
  )

print(
  doc,
  target = "ORACLE_multidomain_response_patterns.docx"
)


## Among patients who respond in 2 domains, are some combinations more common than others?

two_domain <- ORACLE_all3 %>%
  dplyr::filter(N_domains_response == 2)

two_domain %>%
  dplyr::count(Response_pattern) %>%
  dplyr::mutate(
    Percent = 100 * n / sum(n)
  )

## if you responded in 2 domains, were you more likely to respond in 3 domains?

# ============================================================
# Create combined responder variables
# ============================================================

ORACLE_all3 <- ORACLE_all3 %>%
  dplyr::mutate(
    Responder_Attack_Lung =
      Responder_Attack & Responder_Lung,
    
    Responder_Attack_ACQ =
      Responder_Attack & Responder_ACQ,
    
    Responder_Lung_ACQ =
      Responder_Lung & Responder_ACQ
  )

# ============================================================
# Function to create one 2 x 2 table
# ============================================================

.make_crosstable <- function(
    data,
    row_var,
    col_var,
    row_label,
    col_label) {
  
  # 2 x 2 contingency table
  tab <- table(
    data[[row_var]],
    data[[col_var]]
  )
  
  # Row percentages
  prop <- prop.table(tab, margin = 1) * 100
  
  # Fisher's exact test
  p_value <- fisher.test(tab)$p.value
  
  # Create formatted table
  tibble::tibble(
    
    `Response status` = c(
      paste0("No ", row_label),
      paste0(row_label)
    ),
    
    `No response` = c(
      paste0(
        tab[1, 1],
        " (",
        sprintf("%.1f", prop[1, 1]),
        "%)"
      ),
      paste0(
        tab[2, 1],
        " (",
        sprintf("%.1f", prop[2, 1]),
        "%)"
      )
    ),
    
    `Response` = c(
      paste0(
        tab[1, 2],
        " (",
        sprintf("%.1f", prop[1, 2]),
        "%)"
      ),
      paste0(
        tab[2, 2],
        " (",
        sprintf("%.1f", prop[2, 2]),
        "%)"
      )
    ),
    
    `P-value` = c(
      sprintf("%.3f", p_value),
      ""
    )
  )
}

# ============================================================
# 1. ACQ-5 response vs Attack + Lung response
# ============================================================

ct1 <- .make_crosstable(
  data = ORACLE_all3,
  row_var = "Responder_ACQ",
  col_var = "Responder_Attack_Lung",
  row_label = "ACQ-5 response",
  col_label = "Attack + Lung response"
) %>%
  dplyr::mutate(
    Comparison =
      "ACQ-5 response vs Asthma attacks + Lung function"
  )

# ============================================================
# 2. Attack response vs ACQ + Lung response
# ============================================================

ct2 <- .make_crosstable(
  data = ORACLE_all3,
  row_var = "Responder_Attack",
  col_var = "Responder_Lung_ACQ",
  row_label = "Asthma attack response",
  col_label = "ACQ-5 + Lung response"
) %>%
  dplyr::mutate(
    Comparison =
      "Asthma attack response vs ACQ-5 + Lung function"
  )


# ============================================================
# 3. Lung response vs Attack + ACQ response
# ============================================================

ct3 <- .make_crosstable(
  data = ORACLE_all3,
  row_var = "Responder_Lung",
  col_var = "Responder_Attack_ACQ",
  row_label = "Lung function response",
  col_label = "Asthma attacks + ACQ-5 response"
) %>%
  dplyr::mutate(
    Comparison =
      "Lung function response vs Asthma attacks + ACQ-5"
  )

# ============================================================
# Combine the three tables
# ============================================================

Table_3_crosstabs <- dplyr::bind_rows(
  ct1,
  ct2,
  ct3
) %>%
  dplyr::select(
    Comparison,
    `Response status`,
    `No response`,
    Response,
    `P-value`
  )


# View
Table_3_crosstabs

ft_3_crosstabs <- Table_3_crosstabs %>%
  flextable::flextable() %>%
  
  # Clean table style
  flextable::theme_booktabs() %>%
  
  # Bold header
  flextable::bold(
    part = "header"
  ) %>%
  
  # Alignment
  flextable::align(
    j = c(
      "No response",
      "Response",
      "P-value"
    ),
    align = "center",
    part = "all"
  ) %>%
  
  flextable::align(
    j = c(
      "Comparison",
      "Response status"
    ),
    align = "left",
    part = "all"
  ) %>%
  
  # Merge comparison labels
  flextable::merge_v(
    j = "Comparison"
  ) %>%
  
  flextable::valign(
    j = "Comparison",
    valign = "top",
    part = "body"
  ) %>%
  
  # Adjust widths
  flextable::autofit()

ft_3_crosstabs

doc <- officer::read_docx() %>%
  flextable::body_add_flextable(
    value = ft_3_crosstabs
  )

doc <- doc %>%
  officer::body_add_par(
    paste0(
      "Note: Percentages are row percentages. P-values are from ",
      "Fisher's exact tests. The common cohort includes patients ",
      "with ≥120 days of follow-up (N = ",
      nrow(ORACLE_all3),
      ")."
    )
  )

print(
  doc,
  target = "ORACLE_three_domain_response_combinations.docx"
)

###############################################################################
###############################################################################
# Add outcomes in Multivariable forest plot
###############################################################################
###############################################################################

# SET UP THE CPU PROCESSING FOR THE ANALYSES
parallel::detectCores()

## In the function, parallelism in the process or not

#######################################################################################################################################

# Rescale Placebo_response_Lung

ORACLE_Placebo_Attack_imputed$Placebo_response_Lung_per_100mL = ORACLE_Placebo_Attack_imputed$Placebo_response_Lung / 100

# 4.1 Simple model with only Attack_history_num as covariable

# Predictor → random-slope variable mapping
# (scaled/log vars use their raw counterpart; binary/categorical use themselves)
.predictors_simple <- list(
  list(pred = "Age_per_10",                      slope = "Age"),
  list(pred = "Sex",                             slope = "Sex"),
  list(pred = "BMI_per_5",                       slope = "BMI"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no"),
  list(pred = "Airborne_allergen_sensitisation", slope = "Airborne_allergen_sensitisation"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis"),
  list(pred = "Eczema",                          slope = "Eczema"),
  list(pred = "CRSwNP",                          slope = "CRSwNP"),
  list(pred = "CRSsNP",                          slope = "CRSsNP"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num"),  # no covariable for this one
  list(pred = "FEV1_preBD_L_0W_per_10",          slope = "FEV1_preBD_L_0W"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_log10",                     slope = "BEC_log10"),
  list(pred = "FeNO_log10",                    slope = "FeNO_log10"),
  list(pred = "IgE_log10",                       slope = "IgE_log10"),
  list(pred = "Placebo_response_Lung_per_100mL",                       slope = "Placebo_response_Lung_per_100mL"),
  list(pred = "Placebo_response_ACQ",                       slope = "Placebo_response_ACQ")
)

.tables_simple <- vector("list", length(.predictors_simple))

for (i in seq_along(.predictors_simple)) {
  pred  <- .predictors_simple[[i]]$pred
  slope <- .predictors_simple[[i]]$slope
  
  # When the predictor IS Attack_history_num, drop it from covariables
  cov       <- if (pred == "Attack_history_num") NULL else "Attack_history_num"
  cov_slope <- if (pred == "Attack_history_num") NULL else "Attack_history_num"
  
  .model <- IPD_one_stage(
    data                        = ORACLE_Placebo_Attack_imputed,
    outcome_var                 = "Placebo_response_Attack",
    predictor_vars              = pred,
    covariables                 = "Attack_history_num",
    imp_col                     = ".imp",
    followup_offset             = "Yes",
    followup_col                = "Follow_up_duration_days",
    random_intercept_var        = "Enrolled_Trial_name",
    stratified_intercept_var    = "Enrolled_Trial_name",
    predictor_vars_random_slope = slope,
    covariables_random_slope    = "Attack_history_num",
    model_type                  = "lm",
    model_performance           = FALSE,
    weighted_intercept          = FALSE,
    parallel                    = TRUE,
    n_cores                     = 2
  )
  
  # Keep only the row for the predictor of interest
  # startsWith() handles dichotomic variables whose term gains a level suffix (e.g. Sex -> Sex1)
  .tables_simple[[i]] <- .model$table[startsWith(.model$table$term, pred), ]
  rm(.model); gc()
}

Coefficients_Placebo_Attack_simple_model <- dplyr::bind_rows(.tables_simple)

#######################################################################################################################################
#######################################################################################################################################
# FIGURE 2.A - In the manuscript
#######################################################################################################################################
#######################################################################################################################################
# 4.3 Predictors ∆Asthma attack in the final model
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Map every term back to its predictor / slope pair ────────────────────────
.predictor_map_final <- list(
  list(pred = "Age_per_10",                      slope = "Age",                             term = "Age_per_10"),
  list(pred = "Sex",                             slope = "Sex",                             term = "SexMale"),
  list(pred = "BMI_per_5",                       slope = "BMI",                             term = "BMI_per_5"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no",          term = "Smoking_history_yes_noYes"),
  list(pred = "Airborne_allergen_sensitisation", slope = "Airborne_allergen_sensitisation", term = "Airborne_allergen_sensitisationYes"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis",               term = "Allergic_RhinitisYes"),
  list(pred = "Eczema",                          slope = "Eczema",                          term = "EczemaYes"),
  list(pred = "CRSwNP",                          slope = "CRSwNP",                          term = "CRSwNPYes"),
  list(pred = "CRSsNP",                          slope = "CRSsNP",                          term = "CRSsNPYes"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease",             term = "Psychiatric_diseaseYes"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric",               term = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W",                    term = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num",              term = "Attack_history_num"),
  list(pred = "FEV1_preBD_L_0W_per_10",          slope = "FEV1_preBD_L_0W",                 term = "FEV1_preBD_L_0W_per_10"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10",    term = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_log10",                     slope = "BEC_log10",                     term = "BEC_log10"),
  list(pred = "FeNO_log10",                    slope = "FeNO_log10",                    term = "FeNO_log10"),
  list(pred = "IgE_log10",                       slope = "IgE_log10",                       term = "IgE_log10"),
  list(pred = "Placebo_response_Lung_per_100mL",                       slope = "Placebo_response_Lung_per_100mL",                       term = "Placebo_response_Lung_per_100mL"),
  list(pred = "Placebo_response_ACQ",                       slope = "Placebo_response_ACQ",                       term = "Placebo_response_ACQ")
)

# ── Build covariate pool from predictors significant in the simple models ─────
# Each predictor is tested individually, adjusted for all significant predictors
# except itself (to avoid self-adjustment).
.sig_terms <- Coefficients_Placebo_Attack_simple_model |>
  filter(p.value < 0.05) |>
  pull(term)

.sig_entries    <- purrr::keep(.predictor_map_final, ~ .x$term %in% .sig_terms)
.sig_cov_preds  <- purrr::map_chr(.sig_entries, "pred")
.sig_cov_slopes <- purrr::map_chr(.sig_entries, "slope")

# ── Loop: test each predictor adjusted for significant covariables ─────────────
.tables_final <- vector("list", length(.predictor_map_final))
.cov_used     <- vector("list", length(.predictor_map_final))


for (i in seq_along(.predictor_map_final)) {
  pred  <- .predictor_map_final[[i]]$pred
  slope <- .predictor_map_final[[i]]$slope
  
  # Remove current predictor from covariate pool (no self-adjustment)
  .excl     <- which(.sig_cov_preds == pred)
  cov       <- if (length(.excl) > 0) .sig_cov_preds[-.excl]  else .sig_cov_preds
  #cov_slope <- if (length(.excl) > 0) .sig_cov_slopes[-.excl] else .sig_cov_slopes
  
  cov       <- if (length(cov) == 0) NULL else cov
  .cov_used[[i]] <- list(predictor = pred, covariates = cov) 
  #cov_slope <- if (length(cov_slope) == 0) NULL else cov_slope
  
  .model <- IPD_one_stage(
    data                        = ORACLE_Placebo_Attack_imputed,
    outcome_var                 = "Placebo_response_Attack",
    predictor_vars              = pred,
    covariables                 = cov,
    imp_col                     = ".imp",
    followup_offset             = "Yes",
    followup_col                = "Follow_up_duration_days",
    random_intercept_var        = "Enrolled_Trial_name",
    stratified_intercept_var    = "Enrolled_Trial_name",
    predictor_vars_random_slope = slope,
    covariables_random_slope    = "Attack_history_num",
    model_type                  = "lm",
    model_performance           = FALSE,
    weighted_intercept          = FALSE,
    parallel                    = FALSE,
    n_cores                     = 1
  )
  
  .tables_final[[i]] <- .model$table[startsWith(.model$table$term, pred), ]
  rm(.model); gc()
}

Coefficients_Placebo_Attack_final_model <- dplyr::bind_rows(.tables_final)
rm(.tables_final, .sig_entries, .predictor_map_final)
Coefficients_Placebo_Attack_final_model

# Check for the correct Covariates (in 'cov') - Should significant covariates from Model 1 (Simple model)
cov_table <- purrr::map_dfr(.cov_used, function(x) {
  tibble::tibble(
    predictor  = x$predictor,
    covariates = if (is.null(x$covariates)) "none" else paste(x$covariates, collapse = ", ")
  )
})

#######################################################################################################################################
# 4.4.  Forest plot – Final model (∆ Asthma Attack Rate)
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Label and category mapping ────────────────────────────────────────────────
df_plot <- Coefficients_Placebo_Attack_final_model |>
  mutate(
    label = case_match(
      term,
      "Age_per_10"                         ~ "Age (per 10 years)",
      "SexMale"                            ~ "Sex (Male vs. Female)",
      "BMI_per_5"                          ~ "BMI (per 5 kg/m\u00b2)",
      "Smoking_history_yes_noYes"          ~ "Smoking history",
      "Airborne_allergen_sensitisationYes" ~ "Airborne allergen sensitisation",
      "Allergic_RhinitisYes"              ~ "Allergic rhinitis",
      "EczemaYes"                          ~ "Eczema",
      "CRSwNPYes"                          ~ "CRS with nasal polyps",
      "CRSsNPYes"                          ~ "CRS without nasal polyps",
      "Psychiatric_diseaseYes"             ~ "Psychiatric disease",
      "GINA_step_numeric"                  ~ "GINA step",
      "ACQ_score_0W"                       ~ "ACQ score (baseline)",
      "Attack_history_num"                 ~ "Exacerbation history (count)",
      "FEV1_preBD_L_0W_per_10"             ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"       ~ "FEV1 reversibility (per 10% increase)",
      "BEC_log10"                          ~ "Blood eosinophils (per 10 fold increase)",
      "FeNO_log10"                         ~ "FeNO (per 10 fold increase)",
      "IgE_log10"                          ~ "Total IgE (log10)",
      'Placebo_response_Lung_per_100mL' ~ "Delta FEV1 (per 100 mL)",
      'Placebo_response_ACQ' ~ "Delta ACQ-5"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5","Smoking_history_yes_noYes")
      ~ "Demographic",
      c("Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline lung function",
      c("BEC_log10", "FeNO_log10", "IgE_log10",'Placebo_response_Lung_per_100mL','Placebo_response_ACQ')
      ~ "Inflammatory biomarkers"
    ),
    category    = factor(category, levels = c(
      "Demographic", "Comorbidities", "Asthma history",
      "Baseline lung function", "Inflammatory biomarkers"
    )),
    # Drop unused category levels (categories absent from final model)
    category    = droplevels(category),
    label       = factor(label, levels = rev(label)),
    p_label     = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p.value)),
    rc_label    = sprintf("%.2f [%.2f; %.2f]", estimate, `2.5 %`, `97.5 %`),
    significant = p.value < 0.05
  )

# ── Build flat data frame: spacer + header rows inserted per category ─────────
.cats     <- levels(df_plot$category)
.all_rows <- list()

for (.cat in .cats) {
  if (.cat != .cats[1]) {
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = paste0(".spacer.", gsub(" ", "_", .cat)), type = "spacer",
      row_label = "", estimate = NA_real_, lo = NA_real_, hi = NA_real_,
      p_label = "", rc_label = "", significant = FALSE
    )
  }
  .all_rows[[length(.all_rows) + 1]] <- tibble(
    y_id = paste0(".header.", gsub(" ", "_", .cat)), type = "header",
    row_label = .cat, estimate = NA_real_, lo = NA_real_, hi = NA_real_,
    p_label = "", rc_label = "", significant = FALSE
  )
  .cat_vars <- df_plot |> filter(category == .cat)
  for (.k in seq_len(nrow(.cat_vars))) {
    .v <- .cat_vars[.k, ]
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = as.character(.v$label), type = "variable",
      row_label = as.character(.v$label),
      estimate = .v$estimate, lo = .v$`2.5 %`, hi = .v$`97.5 %`,
      p_label = .v$p_label, rc_label = .v$rc_label, significant = .v$significant
    )
  }
}

.df_flat <- bind_rows(.all_rows) |>
  mutate(y_id = factor(y_id, levels = rev(y_id)))

.df_vars <- .df_flat |> filter(type == "variable")

# Bold headers via ggtext markdown, plain variables, blank spacers
.axis_labels <- setNames(
  case_when(
    .df_flat$type == "header"   ~ paste0("**", .df_flat$row_label, "**"),
    .df_flat$type == "variable" ~ .df_flat$row_label,
    TRUE                         ~ ""
  ),
  as.character(.df_flat$y_id)
)

# ── Forest plot ───────────────────────────────────────────────────────────────
.x_pval <- 1.35
.n_rows  <- nrow(.df_flat)

.p_forest <- ggplot(.df_flat, aes(y = y_id)) +
  geom_blank(aes(x = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_errorbarh(
    data = .df_vars,
    aes(xmin = lo, xmax = hi, x = estimate),
    height = 0.3, linewidth = 0.55, color = "gray30"
  ) +
  geom_point(
    data = .df_vars,
    aes(x = estimate, fill = significant),
    shape = 22, size = 3, color = "gray20"
  ) +
  geom_text(
    data = .df_vars,
    aes(x = .x_pval, label = rc_label),
    size = 3, hjust = 0.5
  ) +
  annotate("text",
           x = .x_pval, y = .n_rows + 1,
           label = "aRC [CI95]", fontface = "bold", size = 3.2, hjust = 0.5
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "steelblue4"), guide = "none") +
  scale_x_continuous(breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)) +
  scale_y_discrete(labels = .axis_labels) +
  coord_cartesian(xlim = c(-1, 1), clip = "off") +
  labs(
    x = expression(Delta ~ "Asthma Attack Rate"),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y   = element_markdown(size = 9, hjust = 1, color = "black"),
    plot.margin   = margin(20, 115, 0, 5),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ── Directional arrow strip (below x-axis title) ──────────────────────────────
.p_arrows <- ggplot() +
  annotate("segment",
           x = -0.02, xend = -1, y = 0.85, yend = 0.85,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = -0.5, y = 0.40, label = "Reduced placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  annotate("segment",
           x = 0.02, xend = 1, y = 0.85, yend = 0.85,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = 0.5, y = 0.40, label = "Greater placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(0, 85, 5, 5))

# ── Combine and save ──────────────────────────────────────────────────────────
Figure_Placebo_Attack_final_forest <- (.p_forest / .p_arrows) + plot_layout(heights = c(20, 1))
Figure_Placebo_Attack_final_forest

ggsave(
  filename = file.path(Figure_path, "Figure_Placebo_Attack_final_forest_With_outcomes.pdf"),
  plot     = Figure_Placebo_Attack_final_forest,
  width    = 10, height = 6, units = "in"
)








































#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# PART B : Models for ∆ FEV1
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

# 4.1 Simple model with only FEV1_preBD_L_0W_per_10 as covariable

# Predictor → random-slope variable mapping
# (scaled/log vars use their raw counterpart; binary/categorical use themselves)
.predictors_simple <- list(
  list(pred = "Age_per_10",                      slope = "Age"),
  list(pred = "Sex",                             slope = "Sex"),
  list(pred = "BMI_per_5",                       slope = "BMI"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no"),
  list(pred = "Airborne_allergen_sensitisation", slope = "Airborne_allergen_sensitisation"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis"),
  list(pred = "Eczema",                          slope = "Eczema"),
  list(pred = "CRSwNP",                          slope = "CRSwNP"),
  list(pred = "CRSsNP",                          slope = "CRSsNP"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num"),  # no covariable for this one
  list(pred = "FEV1_preBD_L_0W_per_10",          slope = "FEV1_preBD_L_0W"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_log10",                     slope = "BEC_log10"),
  list(pred = "FeNO_log10",                    slope = "FeNO_log10"),
  list(pred = "IgE_log10",                       slope = "IgE_log10"),
  list(pred = "Placebo_response_Attack",                       slope = "Placebo_response_Attack"),
  list(pred = "Placebo_response_ACQ",                       slope = "Placebo_response_ACQ")
)

.tables_simple <- vector("list", length(.predictors_simple))

for (i in seq_along(.predictors_simple)) {
  pred  <- .predictors_simple[[i]]$pred
  slope <- .predictors_simple[[i]]$slope
  
  # When the predictor IS FEV1_preBD_L_0W_per_10, drop it from covariables
  cov       <- if (pred == "FEV1_preBD_L_0W_per_10") NULL else "FEV1_preBD_L_0W_per_10"
  #cov_slope <- if (pred == "FEV1_preBD_L_0W_per_10") NULL else "FEV1_preBD_L_0W_per_10"
  
  .model <- IPD_one_stage(
    data                        = ORACLE_Placebo_Lung_imputed,
    outcome_var                 = "Placebo_response_Lung",
    predictor_vars              = pred,
    covariables                 = "FEV1_preBD_L_0W_per_10",
    imp_col                     = ".imp",
    #followup_offset = "Yes",
    #followup_col    = "Follow_up_duration_days",
    random_intercept_var        = "Enrolled_Trial_name",
    stratified_intercept_var    = "Enrolled_Trial_name",
    predictor_vars_random_slope = slope,
    covariables_random_slope    = "FEV1_preBD_L_0W_per_10",
    model_type                  = "lm",
    model_performance           = FALSE,
    weighted_intercept          = FALSE,
    parallel                    = TRUE,
    n_cores                     = 3
  )
  
  # Keep only the row for the predictor of interest
  # startsWith() handles dichotomic variables whose term gains a level suffix (e.g. Sex -> Sex1)
  .tables_simple[[i]] <- .model$table[startsWith(.model$table$term, pred), ]
  rm(.model); gc()
}

Coefficients_Placebo_Lung_simple_model <- dplyr::bind_rows(.tables_simple)

#######################################################################################################################################
#######################################################################################################################################
# FIGURE 2.B - In the manuscript
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 4.3 Predictors ∆Lung in the final model
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Map every term back to its predictor / slope pair ────────────────────────
.predictor_map_final <- list(
  list(pred = "Age_per_10",                      slope = "Age",                             term = "Age_per_10"),
  list(pred = "Sex",                             slope = "Sex",                             term = "SexMale"),
  list(pred = "BMI_per_5",                       slope = "BMI",                             term = "BMI_per_5"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no",          term = "Smoking_history_yes_noYes"),
  list(pred = "Airborne_allergen_sensitisation",  slope = "Airborne_allergen_sensitisation", term = "Airborne_allergen_sensitisationYes"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis",               term = "Allergic_RhinitisYes"),
  list(pred = "Eczema",                          slope = "Eczema",                          term = "EczemaYes"),
  list(pred = "CRSwNP",                          slope = "CRSwNP",                          term = "CRSwNPYes"),
  list(pred = "CRSsNP",                          slope = "CRSsNP",                          term = "CRSsNPYes"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease",             term = "Psychiatric_diseaseYes"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric",               term = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W",                    term = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num",              term = "Attack_history_num"),
  list(pred = "FEV1_preBD_L_0W_per_10",          slope = "FEV1_preBD_L_0W",                   term = "FEV1_preBD_L_0W_per_10"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10",    term = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_log10",                       slope = "BEC_log10",                   term = "BEC_log10"),
  list(pred = "FeNO_log10",                      slope = "FeNO_log10",                  term = "FeNO_log10"),
  list(pred = "IgE_log10",                       slope = "IgE_log10",                       term = "IgE_log10"),
  list(pred = "Placebo_response_Attack",                      slope = "Placebo_response_Attack",                  term = "Placebo_response_Attack"),
  list(pred = "Placebo_response_ACQ",                       slope = "Placebo_response_ACQ",                       term = "Placebo_response_ACQ")
)

# ── Build covariate pool from predictors significant in the simple models ─────
# Each predictor is tested individually, adjusted for all significant predictors
# except itself (to avoid self-adjustment).
.sig_terms <- Coefficients_Placebo_Lung_simple_model |>
  filter(p.value < 0.05) |>
  pull(term)

.sig_entries    <- purrr::keep(.predictor_map_final, ~ .x$term %in% .sig_terms)
.sig_cov_preds  <- purrr::map_chr(.sig_entries, "pred")
.sig_cov_slopes <- purrr::map_chr(.sig_entries, "slope")

# ── Loop: test each predictor adjusted for significant covariables ─────────────
.tables_final <- vector("list", length(.predictor_map_final))
.cov_used     <- vector("list", length(.predictor_map_final))

for (i in seq_along(.predictor_map_final)) {
  pred  <- .predictor_map_final[[i]]$pred
  slope <- .predictor_map_final[[i]]$slope
  
  # Remove current predictor from covariate pool (no self-adjustment)
  .excl     <- which(.sig_cov_preds == pred)
  cov       <- if (length(.excl) > 0) .sig_cov_preds[-.excl]  else .sig_cov_preds
  #cov_slope <- if (length(.excl) > 0) .sig_cov_slopes[-.excl] else .sig_cov_slopes
  
  cov       <- if (length(cov) == 0) NULL else cov
  #cov_slope <- if (length(cov_slope) == 0) NULL else cov_slope
  .cov_used[[i]] <- list(predictor = pred, covariates = cov)
  
  .model <- IPD_one_stage(
    data                        = ORACLE_Placebo_Lung_imputed,,
    outcome_var                 = "Placebo_response_Lung",
    predictor_vars              = pred,
    covariables                 = cov,
    imp_col                     = ".imp",
    #followup_offset             = "Yes",
    #followup_col                = "Follow_up_duration_days",
    random_intercept_var        = "Enrolled_Trial_name",
    stratified_intercept_var    = "Enrolled_Trial_name",
    predictor_vars_random_slope = slope,
    covariables_random_slope    = "FEV1_preBD_L_0W_per_10",
    model_type                  = "lm",
    model_performance           = FALSE,
    weighted_intercept          = FALSE,
    parallel                    = TRUE,
    n_cores                     = 2
  )
  
  .tables_final[[i]] <- .model$table[startsWith(.model$table$term, pred), ]
  rm(.model); gc()
}

Coefficients_Placebo_Lung_final_model <- dplyr::bind_rows(.tables_final)
rm(.tables_final, .sig_entries, .predictor_map_final)
Coefficients_Placebo_Lung_final_model


cov_table <- purrr::map_dfr(.cov_used, function(x) {
  tibble::tibble(
    predictor  = x$predictor,
    covariates = if (is.null(x$covariates)) "none" else paste(x$covariates, collapse = ", ")
  )
})

#######################################################################################################################################
# 4.4.  Forest plot – Final model (∆ Lung)
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Label and category mapping ────────────────────────────────────────────────
df_plot <- Coefficients_Placebo_Lung_final_model |>
  mutate(
    label = case_match(
      term,
      "Age_per_10"                         ~ "Age (per 10 years)",
      "SexMale"                            ~ "Sex (Male vs. Female)",
      "BMI_per_5"                          ~ "BMI (per 5 kg/m\u00b2)",
      "Smoking_history_yes_noYes"          ~ "Smoking history",
      "Airborne_allergen_sensitisationYes" ~ "Airborne allergen sensitisation",
      "Allergic_RhinitisYes"              ~ "Allergic rhinitis",
      "EczemaYes"                          ~ "Eczema",
      "CRSwNPYes"                          ~ "CRS with nasal polyps",
      "CRSsNPYes"                          ~ "CRS without nasal polyps",
      "Psychiatric_diseaseYes"             ~ "Psychiatric disease",
      "GINA_step_numeric"                  ~ "GINA step",
      "ACQ_score_0W"                       ~ "ACQ score (baseline)",
      "Attack_history_num"                 ~ "Exacerbation history (count)",
      "FEV1_preBD_L_0W_per_10"             ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"       ~ "FEV1 reversibility (per 10% increase)",
      "BEC_log10"                        ~ "Blood eosinophils (per 10 fold increase)",
      "FeNO_log10"                       ~ "FeNO (per 10 fold increase)",
      "IgE_log10"                          ~ "Total IgE (log10)",
      "Placebo_response_Attack"                       ~ "Delta number of attack",
      "Placebo_response_ACQ"                          ~ "Delta ACQ-5"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5","Smoking_history_yes_noYes")
      ~ "Demographic",
      c("Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline lung function",
      c("BEC_log10", "FeNO_log10", "IgE_log10","Placebo_response_Attack", "Placebo_response_ACQ" )
      ~ "Inflammatory biomarkers"
    ),
    category    = factor(category, levels = c(
      "Demographic", "Comorbidities", "Asthma history",
      "Baseline lung function", "Inflammatory biomarkers"
    )),
    # Drop unused category levels (categories absent from final model)
    category    = droplevels(category),
    label       = factor(label, levels = rev(label)),
    p_label     = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p.value)),
    rc_label    = sprintf("%.2f [%.2f; %.2f]", estimate, `2.5 %`, `97.5 %`),
    significant = p.value < 0.05
  )

# ── Build flat data frame: spacer + header rows inserted per category ─────────
.cats     <- levels(df_plot$category)
.all_rows <- list()

for (.cat in .cats) {
  if (.cat != .cats[1]) {
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = paste0(".spacer.", gsub(" ", "_", .cat)), type = "spacer",
      row_label = "", estimate = NA_real_, lo = NA_real_, hi = NA_real_,
      p_label = "", rc_label = "", significant = FALSE
    )
  }
  .all_rows[[length(.all_rows) + 1]] <- tibble(
    y_id = paste0(".header.", gsub(" ", "_", .cat)), type = "header",
    row_label = .cat, estimate = NA_real_, lo = NA_real_, hi = NA_real_,
    p_label = "", rc_label = "", significant = FALSE
  )
  .cat_vars <- df_plot |> filter(category == .cat)
  for (.k in seq_len(nrow(.cat_vars))) {
    .v <- .cat_vars[.k, ]
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = as.character(.v$label), type = "variable",
      row_label = as.character(.v$label),
      estimate = .v$estimate, lo = .v$`2.5 %`, hi = .v$`97.5 %`,
      p_label = .v$p_label, rc_label = .v$rc_label, significant = .v$significant
    )
  }
}

.df_flat <- bind_rows(.all_rows) |>
  mutate(y_id = factor(y_id, levels = rev(y_id)))

.df_vars <- .df_flat |> filter(type == "variable")

# Bold headers via ggtext markdown, plain variables, blank spacers
.axis_labels <- setNames(
  case_when(
    .df_flat$type == "header"   ~ paste0("**", .df_flat$row_label, "**"),
    .df_flat$type == "variable" ~ .df_flat$row_label,
    TRUE                         ~ ""
  ),
  as.character(.df_flat$y_id)
)

# ── Forest plot ───────────────────────────────────────────────────────────────
.x_pval <- 540   # data-unit x for p-value column (outside xlim of -400 to 400)
.n_rows  <- nrow(.df_flat)

.p_forest <- ggplot(.df_flat, aes(y = y_id)) +
  geom_blank(aes(x = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_errorbarh(
    data = .df_vars,
    aes(xmin = lo, xmax = hi, x = estimate),
    height = 0.3, linewidth = 0.55, color = "gray30"
  ) +
  geom_point(
    data = .df_vars,
    aes(x = estimate, fill = significant),
    shape = 22, size = 3, color = "gray20"
  ) +
  geom_text(
    data = .df_vars,
    aes(x = .x_pval, label = rc_label),
    size = 3, hjust = 0.5
  ) +
  annotate("text",
           x = .x_pval, y = .n_rows + 1,
           label = "aRC [CI95]", fontface = "bold", size = 3.2, hjust = 0.5
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "steelblue4"), guide = "none") +
  scale_x_continuous(breaks = seq(-400, 400, by = 100)) +
  scale_y_discrete(labels = .axis_labels) +
  coord_cartesian(xlim = c(-400, 400), clip = "off") +
  labs(
    x = expression(Delta ~ "FEV1 (mL)"),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y   = element_markdown(size = 9, hjust = 1, color = "black"),
    plot.margin   = margin(20, 115, 0, 5),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ── Directional arrow strip (below x-axis title) ──────────────────────────────
.p_arrows <- ggplot() +
  annotate("segment",
           x = -10, xend = -400, y = 0.85, yend = 0.85,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = -200, y = 0.40, label = "Reduced placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  annotate("segment",
           x = 10, xend = 400, y = 0.85, yend = 0.85,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = 200, y = 0.40, label = "Greater placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  scale_x_continuous(limits = c(-440, 440)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(0, 85, 5, 5))

# ── Combine and save ──────────────────────────────────────────────────────────
Figure_Placebo_Lung_final_forest <- (.p_forest / .p_arrows) + plot_layout(heights = c(20, 1))
Figure_Placebo_Lung_final_forest

ggsave(
  filename = file.path(Figure_path, "Figure_Placebo_Lung_final_forest_IQR_with_outcomes.pdf"),
  plot     = Figure_Placebo_Lung_final_forest,
  width    = 10, height = 6, units = "in"
)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# PART C : Models for ∆ ACQ
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

# 4.	Models for ∆ ACQ

# 4.1 Simple model with only ACQ_score_0W as covariable

# Predictor → random-slope variable mapping
# (scaled/log vars use their raw counterpart; binary/categorical use themselves)
.predictors_simple <- list(
  list(pred = "Age_per_10",                      slope = "Age"),
  list(pred = "Sex",                             slope = "Sex"),
  list(pred = "BMI_per_5",                       slope = "BMI"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no"),
  list(pred = "Airborne_allergen_sensitisation", slope = "Airborne_allergen_sensitisation"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis"),
  list(pred = "Eczema",                          slope = "Eczema"),
  list(pred = "CRSwNP",                          slope = "CRSwNP"),
  list(pred = "CRSsNP",                          slope = "CRSsNP"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num"),  # no covariable for this one
  list(pred = "FEV1_preBD_L_0W_per_10",          slope = "FEV1_preBD_L_0W"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_log10",                     slope = "BEC_log10"),
  list(pred = "FeNO_log10",                    slope = "FeNO_log10"),
  list(pred = "IgE_log10",                       slope = "IgE_log10"),
  list(pred = "Placebo_response_Attack",                    slope = "Placebo_response_Attack"),
  list(pred = "Placebo_response_Lung_per_100mL",                       slope = "Placebo_response_Lung_per_100mL")
)

.tables_simple <- vector("list", length(.predictors_simple))

for (i in seq_along(.predictors_simple)) {
  pred  <- .predictors_simple[[i]]$pred
  slope <- .predictors_simple[[i]]$slope
  
  # When the predictor IS ACQ_score_0W, drop it from covariables
  cov       <- if (pred == "ACQ_score_0W") NULL else "ACQ_score_0W"
  #cov_slope <- if (pred == "ACQ_score_0W") NULL else "ACQ_score_0W"
  
  .model <- IPD_one_stage(
    data                        = ORACLE_Placebo_ACQ_imputed,
    outcome_var                 = "Placebo_response_ACQ",
    predictor_vars              = pred,
    covariables                 = "ACQ_score_0W",
    imp_col                     = ".imp",
    #followup_offset = "Yes",
    #followup_col    = "Follow_up_duration_days",
    random_intercept_var        = "Enrolled_Trial_name",
    stratified_intercept_var    = "Enrolled_Trial_name",
    predictor_vars_random_slope = slope,
    covariables_random_slope    = "ACQ_score_0W",
    model_type                  = "lm",
    model_performance           = FALSE,
    weighted_intercept          = FALSE,
    parallel                    = TRUE,
    n_cores                     = 2
  )
  
  # Keep only the row for the predictor of interest
  # startsWith() handles dichotomic variables whose term gains a level suffix (e.g. Sex -> Sex1)
  .tables_simple[[i]] <- .model$table[startsWith(.model$table$term, pred), ]
  rm(.model); gc()
}

Coefficients_Placebo_ACQ_simple_model <- dplyr::bind_rows(.tables_simple)

#######################################################################################################################################
#######################################################################################################################################
# FIGURE 2.C - In the manuscript
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 4.3 Predictors ∆ACQ in the final model
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Map every term back to its predictor / slope pair ────────────────────────
.predictor_map_final <- list(
  list(pred = "Age_per_10",                      slope = "Age",                             term = "Age_per_10"),
  list(pred = "Sex",                             slope = "Sex",                             term = "SexMale"),
  list(pred = "BMI_per_5",                       slope = "BMI",                             term = "BMI_per_5"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no",          term = "Smoking_history_yes_noYes"),
  list(pred = "Airborne_allergen_sensitisation",  slope = "Airborne_allergen_sensitisation", term = "Airborne_allergen_sensitisationYes"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis",               term = "Allergic_RhinitisYes"),
  list(pred = "Eczema",                          slope = "Eczema",                          term = "EczemaYes"),
  list(pred = "CRSwNP",                          slope = "CRSwNP",                          term = "CRSwNPYes"),
  list(pred = "CRSsNP",                          slope = "CRSsNP",                          term = "CRSsNPYes"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease",             term = "Psychiatric_diseaseYes"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric",               term = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W",                    term = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num",              term = "Attack_history_num"),
  list(pred = "FEV1_preBD_L_0W_per_10",        slope = "FEV1_preBD_L_0W",               term = "FEV1_preBD_L_0W_per_10"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10",    term = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_log10",                       slope = "BEC_log10",                       term = "BEC_log10"),
  list(pred = "FeNO_log10",                      slope = "FeNO_log10",                      term = "FeNO_log10"),
  list(pred = "IgE_log10",                       slope = "IgE_log10",                       term = "IgE_log10"),
  list(pred = "Placebo_response_Attack",                      slope = "Placebo_response_Attack",                      term = "Placebo_response_Attack"),
  list(pred = "Placebo_response_Lung_per_100mL",                       slope = "Placebo_response_Lung_per_100mL",                       term = "Placebo_response_Lung_per_100mL")
)

# ── Build covariate pool from predictors significant in the simple models ─────
# Each predictor is tested individually, adjusted for all significant predictors
# except itself (to avoid self-adjustment).
.sig_terms <- Coefficients_Placebo_ACQ_simple_model |>
  filter(p.value < 0.05) |>
  pull(term)

.sig_entries    <- purrr::keep(.predictor_map_final, ~ .x$term %in% .sig_terms)
.sig_cov_preds  <- purrr::map_chr(.sig_entries, "pred")
.sig_cov_slopes <- purrr::map_chr(.sig_entries, "slope")

# ── Loop: test each predictor adjusted for significant covariables ─────────────
.tables_final <- vector("list", length(.predictor_map_final))
.cov_used     <- vector("list", length(.predictor_map_final))

for (i in seq_along(.predictor_map_final)) {
  pred  <- .predictor_map_final[[i]]$pred
  slope <- .predictor_map_final[[i]]$slope
  
  # Remove current predictor from covariate pool (no self-adjustment)
  .excl     <- which(.sig_cov_preds == pred)
  cov       <- if (length(.excl) > 0) .sig_cov_preds[-.excl]  else .sig_cov_preds
  #cov_slope <- if (length(.excl) > 0) .sig_cov_slopes[-.excl] else .sig_cov_slopes
  
  cov       <- if (length(cov) == 0) NULL else cov
  #cov_slope <- if (length(cov_slope) == 0) NULL else cov_slope
  .cov_used[[i]] <- list(predictor = pred, covariates = cov)
  
  .model <- IPD_one_stage(
    data                        = ORACLE_Placebo_ACQ_imputed,,
    outcome_var                 = "Placebo_response_ACQ",
    predictor_vars              = pred,
    covariables                 = cov,
    imp_col                     = ".imp",
    #followup_offset             = "Yes",
    #followup_col                = "Follow_up_duration_days",
    random_intercept_var        = "Enrolled_Trial_name",
    stratified_intercept_var    = "Enrolled_Trial_name",
    predictor_vars_random_slope = slope,
    covariables_random_slope    = "ACQ_score_0W",
    model_type                  = "lm",
    model_performance           = FALSE,
    weighted_intercept          = FALSE,
    parallel                    = TRUE,
    n_cores                     = 2
  )
  
  .tables_final[[i]] <- .model$table[startsWith(.model$table$term, pred), ]
  rm(.model); gc()
}

Coefficients_Placebo_ACQ_final_model <- dplyr::bind_rows(.tables_final)
rm(.tables_final, .sig_entries, .predictor_map_final)
Coefficients_Placebo_ACQ_final_model

cov_table <- purrr::map_dfr(.cov_used, function(x) {
  tibble::tibble(
    predictor  = x$predictor,
    covariates = if (is.null(x$covariates)) "none" else paste(x$covariates, collapse = ", ")
  )
})

#######################################################################################################################################
# 4.4.  Forest plot – Final model (∆ ACQ)
#--------------------------------------------------------------------------------------------------------------------------------------

# ── Label and category mapping ────────────────────────────────────────────────
df_plot <- Coefficients_Placebo_ACQ_final_model |>
  mutate(
    label = case_match(
      term,
      "Age_per_10"                         ~ "Age (per 10 years)",
      "SexMale"                            ~ "Sex (Male vs. Female)",
      "BMI_per_5"                          ~ "BMI (per 5 kg/m\u00b2)",
      "Smoking_history_yes_noYes"          ~ "Smoking history",
      "Airborne_allergen_sensitisationYes" ~ "Airborne allergen sensitisation",
      "Allergic_RhinitisYes"               ~ "Allergic rhinitis",
      "EczemaYes"                          ~ "Eczema",
      "CRSwNPYes"                          ~ "CRS with nasal polyps",
      "CRSsNPYes"                          ~ "CRS without nasal polyps",
      "Psychiatric_diseaseYes"             ~ "Psychiatric disease",
      "GINA_step_numeric"                  ~ "GINA step",
      "ACQ_score_0W"                       ~ "ACQ score (baseline)",
      "Attack_history_num"                 ~ "Exacerbation history (count)",
      "FEV1_preBD_L_0W_per_10"             ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"       ~ "FEV1 reversibility (per 10% increase)",
      "BEC_log10"                        ~ "Blood eosinophils (per 10 fold increase)",
      "FeNO_log10"                       ~ "FeNO (per 10 fold increase)",
      "IgE_log10"                          ~ "Total IgE (log10)",
      "Placebo_response_Attack"                       ~ "Delta number of attack",
      "Placebo_response_Lung_per_100mL"                          ~ "Delta FEV1 (per 100 mL)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5","Smoking_history_yes_noYes")
      ~ "Demographic",
      c("Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline ACQ function",
      c("BEC_log10", "FeNO_log10", "IgE_log10","Placebo_response_Attack"   ,"Placebo_response_Lung_per_100mL" )
      ~ "Inflammatory biomarkers"
    ),
    category    = factor(category, levels = c(
      "Demographic", "Comorbidities", "Asthma history",
      "Baseline ACQ function", "Inflammatory biomarkers"
    )),
    # Drop unused category levels (categories absent from final model)
    category    = droplevels(category),
    label       = factor(label, levels = rev(label)),
    p_label     = case_when(p.value < 0.001 ~ "<0.001", TRUE ~ sprintf("%.3f", p.value)),
    rc_label    = sprintf("%.2f [%.2f; %.2f]", estimate, `2.5 %`, `97.5 %`),
    significant = p.value < 0.05
  )

# ── Build flat data frame: spacer + header rows inserted per category ─────────
.cats     <- levels(df_plot$category)
.all_rows <- list()

for (.cat in .cats) {
  if (.cat != .cats[1]) {
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = paste0(".spacer.", gsub(" ", "_", .cat)), type = "spacer",
      row_label = "", estimate = NA_real_, lo = NA_real_, hi = NA_real_,
      p_label = "", rc_label = "", significant = FALSE
    )
  }
  .all_rows[[length(.all_rows) + 1]] <- tibble(
    y_id = paste0(".header.", gsub(" ", "_", .cat)), type = "header",
    row_label = .cat, estimate = NA_real_, lo = NA_real_, hi = NA_real_,
    p_label = "", rc_label = "", significant = FALSE
  )
  .cat_vars <- df_plot |> filter(category == .cat)
  for (.k in seq_len(nrow(.cat_vars))) {
    .v <- .cat_vars[.k, ]
    .all_rows[[length(.all_rows) + 1]] <- tibble(
      y_id = as.character(.v$label), type = "variable",
      row_label = as.character(.v$label),
      estimate = .v$estimate, lo = .v$`2.5 %`, hi = .v$`97.5 %`,
      p_label = .v$p_label, rc_label = .v$rc_label, significant = .v$significant
    )
  }
}

.df_flat <- bind_rows(.all_rows) |>
  mutate(y_id = factor(y_id, levels = rev(y_id)))

.df_vars <- .df_flat |> filter(type == "variable")

# Bold headers via ggtext markdown, plain variables, blank spacers
.axis_labels <- setNames(
  case_when(
    .df_flat$type == "header"   ~ paste0("**", .df_flat$row_label, "**"),
    .df_flat$type == "variable" ~ .df_flat$row_label,
    TRUE                         ~ ""
  ),
  as.character(.df_flat$y_id)
)

# ── Forest plot ───────────────────────────────────────────────────────────────
.x_pval <- 1.35   # data-unit x for p-value column (outside xlim of -1 to 1)
.n_rows  <- nrow(.df_flat)

.p_forest <- ggplot(.df_flat, aes(y = y_id)) +
  geom_blank(aes(x = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_errorbarh(
    data = .df_vars,
    aes(xmin = lo, xmax = hi, x = estimate),
    height = 0.3, linewidth = 0.55, color = "gray30"
  ) +
  geom_point(
    data = .df_vars,
    aes(x = estimate, fill = significant),
    shape = 22, size = 3, color = "gray20"
  ) +
  geom_text(
    data = .df_vars,
    aes(x = .x_pval, label = rc_label),
    size = 3, hjust = 0.5
  ) +
  annotate("text",
           x = .x_pval, y = .n_rows + 1,
           label = "aRC [CI95]", fontface = "bold", size = 3.2, hjust = 0.5
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "steelblue4"), guide = "none") +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.25)) +
  scale_y_discrete(labels = .axis_labels) +
  coord_cartesian(xlim = c(-1, 1), clip = "off") +
  labs(
    x = expression(Delta ~ "ACQ-5 Score"),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y   = element_markdown(size = 9, hjust = 1, color = "black"),
    plot.margin   = margin(20, 115, 0, 5),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ── Directional arrow strip (below x-axis title) ──────────────────────────────
.p_arrows <- ggplot() +
  annotate("segment",
           x = -0.02, xend = -1, y = 0.85, yend = 0.85,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = -0.5, y = 0.40, label = "Reduced placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  annotate("segment",
           x = 0.02, xend = 1, y = 0.85, yend = 0.85,
           arrow = arrow(length = unit(0.18, "cm"), ends = "last", type = "closed"),
           color = "gray30", linewidth = 0.6
  ) +
  annotate("text", x = 0.5, y = 0.40, label = "Greater placebo response",
           size = 2.9, hjust = 0.5, color = "gray30"
  ) +
  scale_x_continuous(limits = c(-1.1, 1.1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(0, 85, 5, 5))

# ── Combine and save ──────────────────────────────────────────────────────────
Figure_Placebo_ACQ_final_forest <- (.p_forest / .p_arrows) + plot_layout(heights = c(20, 1))
Figure_Placebo_ACQ_final_forest

ggsave(
  filename = file.path(Figure_path, "Figure_Placebo_ACQ_final_forest_with_outcomes.pdf"),
  plot     = Figure_Placebo_ACQ_final_forest,
  width    = 10, height = 6, units = "in"
)
