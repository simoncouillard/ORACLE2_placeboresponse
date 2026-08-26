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
library(survivalAnalysis)
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

#######################################################################################################################################
#--------------------------------------------------------------------------------------------------------------------------------------
# 1.2.	Install and load the homemade package MIAnalysis
install.packages("devtools")
library(devtools)
required_packages <- c("dplyr", "MASS", "Hmisc", "mice", "rlang", "survival")
devtools::install_github("SamLar27/MIAnalysis", force = TRUE)

install.packages("doParallel")
devtools::install_github("SamLar27/MIAnalysis", force = TRUE)

library(MIAnalysis)

# VALIDATE : "MIAnalysis" appears in the Packages (section in R right menu) ?

#--------------------------------------------------------------------------------------------------------------------------------------

### DISCLOSURE FOR MIAnalysis (R Package created by SML)
# MIAnalysis is a publicly available R package developed by SML.
# All functions, source code, and documentation are openly accessible,ensuring full transparency and reproducibility of the analyses performed.
# For the package every analytical step can be inspected,audited, and reproduced independently by any user with access to R.
#
# Key functions used in this project:
#   - IPD_one_stage()  : one-stage individual patient data meta-analysis
#                        across trials and multiply imputed datasets
#   - MI_estimates()   : pooling of estimates across imputed datasets
#                        using Rubin's rules
#
# Package availability: [GitHub / CRAN URL]


# Check package info
packageDescription("MIAnalysis")

# Check version
packageVersion("MIAnalysis")

# Full package help
help(package = "MIAnalysis")

# Specific function help
?MI_estimates
?IPD_one_stage

# List all exported functions
ls("package:MIAnalysis")

# View source code of key functions
getAnywhere(MI_estimates)
getAnywhere(IPD_one_stage)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 2.	Defining the path
#--------------------------------------------------------------------------------------------------------------------------------------
# 2.1. Loading path
Load_path <- file.path(
"/Users/joel/ORACLE - Placebo - Exacerbations/ORACLE_Placebo/Data/Raw_data/Main_FCS_single_pmm_logreg_polyreg_pmm"
)

#--------------------------------------------------------------------------------------------------------------------------------------
# 2.2. Saving data path
Saving_path <- file.path(
 "/Users/joel/ORACLE - Placebo - Exacerbations/Raw_data/New_Dataset"
)

#--------------------------------------------------------------------------------------------------------------------------------------
# 2.3. Saving figure path
Figure_path <- file.path(
  "/Users/joel/ORACLE - Placebo - Exacerbations/Project/Figures_Final"
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
unique(ORACLE_after_imputation$Enrolled_Trial_name)
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
  # FEV1 reversibility 10% decrease
  dplyr::mutate(FEV1_reversibility_0W_per_10 = FEV1_reversibility_0W / 10) %>%
  # FVC pre-BD per 10% decrease
  dplyr::mutate(FVC_preBD_PCT_0W_per_10 = -FVC_preBD_PCT_0W / 10) %>%
  # Tif (%) pre-BD per 10% decrease
  dplyr::mutate(Tif_preBD_PCT_0W_per_10 = -Tif_preBD_PCT_0W / 10) %>%
  # Tif (raw) pre-BD per 0.1 decrease in absolute units
  dplyr::mutate(Tif_preBD_0W_per_10 = -Tif_preBD_0W / 0.1) %>%
  # FEV1 pre-BD L per 0.1 decrease in absolute units
  dplyr::mutate(FEV1_preBD_L_0W_per_10 = -FEV1_preBD_L_0W / 0.1) %>%
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

# Flowchart – attack dataset
flowchart_attack <- DiagrammeR::grViz(sprintf("
digraph consort {
  graph [layout = dot, rankdir = TB, fontname = 'Helvetica', splines = ortho, nodesep = 0.6]
  node [shape = rectangle, fontname = 'Helvetica', fontsize = 13,
        margin = '0.3,0.2', width = 5.0, height = 0.8]
  A  [label = 'ORACLE Dataset\\n%d subjects from %d RCTs',
      fillcolor = '#DDEEFF', style = filled]
  B  [label = 'Blinded placebo arm\\n%d subjects from %d RCTs',
      fillcolor = '#DDEEFF', style = filled]
  C  [label = 'Dataset for placebo response\\non asthma exacerbations\\n%d subjects from %d RCTs',
      fillcolor = '#DDEEFF', style = filled]
  E1 [label = 'Excluded (n = %d)\\nNot a blinded placebo arm',
      fillcolor = '#FFDDDD', style = filled, width = 4.5]
  E2 [label = 'Excluded (n = %d)\\nMissing exacerbation history\\n or count during follow-up',
      fillcolor = '#FFDDDD', style = filled, width = 4.5]
  node [shape = point, width = 0.01, style = invis, label = '']
  mid1; mid2
  A    -> mid1 [arrowhead = none]
  mid1 -> B
  mid1 -> E1
  B    -> mid2 [arrowhead = none]
  mid2 -> C
  mid2 -> E2
  { rank = same; mid1; E1 }
  { rank = same; mid2; E2 }
}
",
total_n, total_rcts,
n_placebo, rcts_placebo,
n_attack, rcts_attack,
n_excl_placebo,
n_excl_attack
))
flowchart_attack

svg_tmp <- tempfile(fileext = ".svg")
writeLines(DiagrammeRsvg::export_svg(flowchart_attack), svg_tmp)
rsvg::rsvg_pdf(svg_tmp, file = file.path(Figure_path, "Flowchart_attack.pdf"))

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

# Flowchart – lung function dataset
flowchart_lung <- DiagrammeR::grViz(sprintf("
digraph consort {
  graph [layout = dot, rankdir = TB, fontname = 'Helvetica', splines = ortho, nodesep = 0.6]
  node [shape = rectangle, fontname = 'Helvetica', fontsize = 13,
        margin = '0.3,0.2', width = 5.0, height = 0.8]
  A  [label = 'ORACLE Dataset\\n%d subjects from %d RCTs',
      fillcolor = '#DDEEFF', style = filled]
  B  [label = 'Blinded placebo arm\\n%d subjects from %d RCTs',
      fillcolor = '#DDEEFF', style = filled]
  C  [label = 'Dataset for placebo response\\non Lung Function\\n%d subjects from %d RCTs',
      fillcolor = '#DDEEFF', style = filled]
  E1 [label = 'Excluded (n = %d)\\nNot a blinded placebo arm',
      fillcolor = '#FFDDDD', style = filled, width = 4.5]
  E2 [label = 'Excluded (n = %d)\\nMissing FEV1 at baseline\\nor at 52 wk follow-up',
      fillcolor = '#FFDDDD', style = filled, width = 4.5]
  node [shape = point, width = 0.01, style = invis, label = '']
  mid1; mid2
  A    -> mid1 [arrowhead = none]
  mid1 -> B
  mid1 -> E1
  B    -> mid2 [arrowhead = none]
  mid2 -> C
  mid2 -> E2
  { rank = same; mid1; E1 }
  { rank = same; mid2; E2 }
}
",
total_n, total_rcts,
n_placebo, rcts_placebo,
n_lung, rcts_lung,
n_excl_placebo,
n_excl_lung
))
flowchart_lung

svg_tmp <- tempfile(fileext = ".svg")
writeLines(DiagrammeRsvg::export_svg(flowchart_lung), svg_tmp)
rsvg::rsvg_pdf(svg_tmp, file = file.path(Figure_path, "Flowchart_lung.pdf"))

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

# Flowchart – symptom dataset
flowchart_symptom <- DiagrammeR::grViz(sprintf("
digraph consort {
  graph [layout = dot, rankdir = TB, fontname = 'Helvetica', splines = ortho, nodesep = 0.6]
  node [shape = rectangle, fontname = 'Helvetica', fontsize = 13,
        margin = '0.3,0.2', width = 5.0, height = 0.8]
  A  [label = 'ORACLE Dataset\\n%d subjects from %d RCTs',
      fillcolor = '#DDEEFF', style = filled]
  B  [label = 'Blinded placebo arm\\n%d subjects from %d RCTs',
      fillcolor = '#DDEEFF', style = filled]
  C  [label = 'Dataset for placebo response\\non symptoms control\\n%d subjects from %d RCTs',
      fillcolor = '#DDEEFF', style = filled]
  E1 [label = 'Excluded (n = %d)\\nNot a blinded placebo arm',
      fillcolor = '#FFDDDD', style = filled, width = 4.5]
  E2 [label = 'Excluded (n = %d)\\nMissing ACQ at baseline\\nor at 24-26 wk follow-up',
      fillcolor = '#FFDDDD', style = filled, width = 4.5]
  node [shape = point, width = 0.01, style = invis, label = '']
  mid1; mid2
  A    -> mid1 [arrowhead = none]
  mid1 -> B
  mid1 -> E1
  B    -> mid2 [arrowhead = none]
  mid2 -> C
  mid2 -> E2
  { rank = same; mid1; E1 }
  { rank = same; mid2; E2 }
}
",
total_n, total_rcts,
n_placebo, rcts_placebo,
n_symptom, rcts_symptom,
n_excl_placebo,
n_excl_symptom
))
flowchart_symptom

svg_tmp <- tempfile(fileext = ".svg")
writeLines(DiagrammeRsvg::export_svg(flowchart_symptom), svg_tmp)
rsvg::rsvg_pdf(svg_tmp, file = file.path(Figure_path, "Flowchart_symptom.pdf"))

# Subject counts summary
cat(sprintf(
  "ORACLE_placebo_attack  : %d subjects from %d RCTs\nORACLE_placebo_lung    : %d subjects from %d RCTs\nORACLE_placebo_symptom : %d subjects from %d RCTs\n",
  n_attack,  rcts_attack,
  n_lung,    rcts_lung,
  n_symptom, rcts_symptom
))

#--------------------------------------------------------------------------------------------------------------------------------------
# 6.4  Combined 3-panel flowchart (A / B / C)

# Export each flowchart to a temporary PNG via SVG
.fc_to_png <- function(fc, width = 900) {
  svg_f <- tempfile(fileext = ".svg")
  png_f <- tempfile(fileext = ".png")
  writeLines(DiagrammeRsvg::export_svg(fc), svg_f)
  rsvg::rsvg_png(svg_f, png_f, width = width)
  png_f
}

png_attack  <- .fc_to_png(flowchart_attack)
png_lung    <- .fc_to_png(flowchart_lung)
png_symptom <- .fc_to_png(flowchart_symptom)


# Build one ggplot panel per PNG
.png_panel <- function(path) {
  img <- png::readPNG(path)
  ggplot2::ggplot() +
    ggplot2::annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::theme_void()
}

panel_attack  <- .png_panel(png_attack)
panel_lung    <- .png_panel(png_lung)
panel_symptom <- .png_panel(png_symptom)

# Assemble into a 3-panel figure with A / B / C labels
combined_flowchart <- cowplot::plot_grid(
  panel_attack, panel_lung, panel_symptom,
  labels     = c("A: Asthma Attacks Placebo Flowchart", "B: Lung function Placebo Flowchart", "C: Symptoms Score Placebo Flowchart"),
  label_size = 16,
  label_x    = 0.02,
  label_y    = 0.92,
  hjust      = 0,
  nrow       = 1,
  rel_widths = c(1, 1, 1)
)
combined_flowchart

# Save combined figure as PDF
cowplot::save_plot(
  filename  = file.path(Figure_path, "Flowchart_combined.pdf"),
  plot      = combined_flowchart,
  base_width  = 18,
  base_height = 8
)
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 7.	Generating the final datasets
#     Each dataset is split into a raw (non-imputed, .imp == 0) and
#     an imputed (.imp 1–10) version, based on the subjects selected in section 4.

# 7.1.	Attack datasets
ORACLE_Placebo_Attack_Raw <- ORACLE_placebo_attack |>
  dplyr::filter(.imp == 0)

ORACLE_Placebo_Attack_imputed <- ORACLE_placebo_attack |>
  dplyr::filter(.imp %in% 1:10)
nrow(ORACLE_Placebo_Attack_Raw)

# 7.2.	Lung function datasets
ORACLE_Placebo_Lung_Raw <- ORACLE_placebo_lung |>
  dplyr::filter(.imp == 0)

ORACLE_Placebo_Lung_imputed <- ORACLE_placebo_lung |>
  dplyr::filter(.imp %in% 1:10)
nrow(ORACLE_Placebo_Lung_Raw)

# 7.3.	Symptom (ACQ) datasets
ORACLE_Placebo_ACQ_Raw <- ORACLE_placebo_symptom |>
  dplyr::filter(.imp == 0)

ORACLE_Placebo_ACQ_imputed <- ORACLE_placebo_symptom |>
  dplyr::filter(.imp %in% 1:10)
nrow(ORACLE_Placebo_ACQ_Raw)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

# TABLE 1

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 8.	Creation of Tables

#--------------------------------------------------------------------------------------------------------------------------------------
# 8.1.	Table 1 – Descriptive characteristics across the 3 analysis datasets

# Helper: add ICS indicator columns and build tbl_summary
# Note: var_type is dropped by tbl_merge(), so missing-row removal must happen

summary(ORACLE_Placebo_ACQ_Raw$GINA_step)
# on each individual tbl_summary before merging.
.make_tbl1 <- function(data) {
  gtsummary::tbl_summary(
    data |>
      dplyr::mutate(GINA_step = droplevels(GINA_step)),
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
    type  = list(Sex ~ "dichotomous"),
    value = list(Sex ~ "Female"),
    statistic    = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    missing_text = "(Missing)",
    label = list(
      Age                            ~ "Age (y)",
      Sex                            ~ "Female sex",
      BMI                            ~ "BMI (kg/m²)",
      Smoking_history_yes_no         ~ "Former smokers",
      Psychiatric_disease            ~ "Psychiatric disease",
      Eczema                         ~ "Eczema",
      Allergic_Rhinitis              ~ "Allergic Rhinitis",
      Airborne_allergen_sensitisation ~ "Airborne allergen sensitisation",
      CRSsNP                         ~ "CRSsNP",
      CRSwNP                         ~ "CRSwNP",
      GINA_step                      ~ "GINA treatment step",
      mOCS                           ~ "mOCS",
      Attack_history_cat             ~ "Severe attack in the past 12 mo",
      ICU_or_ETI_history             ~ "History of ICU admission or intubation",
      ACQ_score_0W                   ~ "ACQ-5 score",
      FEV1_preBD_PCT_0W              ~ "FEV1 (% of predicted)",
      FEV1_preBD_L_0W                ~ "FEV1 (L.)",
      Tif_preBD_0W                   ~ "FEV1/FVC",
      FEV1_reversibility_0W          ~ "FEV1 reversibility (%)",
      BEC                            ~ "Blood eosinophil count (×10⁹ cells per L)",
      FeNO                           ~ "FeNO (ppb)",
      Follow_up_duration_days        ~ "Follow-up duration (days)",
      Attack_number_during_followup_0_1_2 ~ "Number of severe attacks during follow-up",
      Placebo_response_Attack        ~ "Change in asthma attack rate by year",
      FEV1_preBD_L_52W               ~ "FEV1 (L) pre-BD at 52 weeks of follow-up",
      Placebo_response_Lung          ~ "FEV1 (ml) change at 52 weeks",
      ACQ_score_26_25_24             ~ "ACQ-5 at 24-26 weeks of follow-up",
      Placebo_response_ACQ           ~ "ACQ-5 score change at 24-26W"
    )
  ) |>
    # Remove missing-count rows for categorical/dichotomous variables here,
    # while var_type is still available (it is dropped by tbl_merge)
    gtsummary::modify_table_body(function(tbl) {
      tbl |>
        dplyr::filter(!(row_type == "missing" & var_type %in% c("categorical", "dichotomous")))
    })
}

tbl1_attack <- .make_tbl1(ORACLE_Placebo_Attack_Raw)
tbl1_lung   <- .make_tbl1(ORACLE_Placebo_Lung_Raw)
tbl1_acq    <- .make_tbl1(ORACLE_Placebo_ACQ_Raw)

# Merge the 3 tables side by side
Table1_Placebo <- gtsummary::tbl_merge(
  tbls        = list(tbl1_attack, tbl1_lung, tbl1_acq),
  tab_spanner = c(
    "**Asthma Exacerbations**",
    "**Lung Function**",
    "**Symptom Control**"
  )
) |>
  gtsummary::modify_table_body(function(tbl) {
    
    # --- Outcome variables: dataset-specific masking ---
    # Variables to replace with "--" per dataset column
    # stat_1_1 = Asthma Exacerbations, stat_1_2 = Lung Function, stat_1_3 = Symptom Control
    all_outcome_vars <- c(
      "Follow_up_duration_days", "Attack_number_during_followup_0_1_2",
      "Placebo_response_Attack", "FEV1_preBD_L_52W",
      "Placebo_response_Lung",   "ACQ_score_26_25_24", "Placebo_response_ACQ"
    )
    vars_hide_attack <- c(
      "FEV1_preBD_L_52W", "Placebo_response_Lung",
      "ACQ_score_26_25_24", "Placebo_response_ACQ"
    )
    vars_hide_lung <- c(
      "Follow_up_duration_days", "Attack_number_during_followup_0_1_2",
      "Placebo_response_Attack", "ACQ_score_26_25_24", "Placebo_response_ACQ"
    )
    vars_hide_acq <- c(
      "Follow_up_duration_days", "Attack_number_during_followup_0_1_2",
      "Placebo_response_Attack", "FEV1_preBD_L_52W", "Placebo_response_Lung"
    )
    
    tbl <- tbl |>
      # Drop (Missing) rows for all outcome variables across all datasets
      dplyr::filter(!(variable %in% all_outcome_vars & row_type == "missing")) |>
      # Replace irrelevant outcome cells with "--"
      dplyr::mutate(
        stat_0_1 = dplyr::if_else(variable %in% vars_hide_attack, "--", stat_0_1),
        stat_0_2 = dplyr::if_else(variable %in% vars_hide_lung,   "--", stat_0_2),
        stat_0_3 = dplyr::if_else(variable %in% vars_hide_acq,    "--", stat_0_3)
      )
    
    # --- Section headers to insert ---
    headers <- list(
      list(before = "Age",                     label = "Demographic"),
      list(before = "Psychiatric_disease",     label = "Comorbidities"),
      list(before = "GINA_step",               label = "Medication"),
      list(before = "Attack_history_cat",      label = "Exacerbation history"),
      list(before = "ACQ_score_0W",            label = "Symptoms"),
      list(before = "FEV1_preBD_L_0W",       label = "Lung function pre-BD"),
      list(before = "BEC",                     label = "Inflammatory biomarkers"),
      list(before = "Follow_up_duration_days", label = "Outcomes/follow-up data")
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
      tbl <- dplyr::bind_rows(tbl[seq_len(idx - 1), ], new_row, tbl[idx:nrow(tbl), ])
    }
    
    tbl |>
      dplyr::mutate(label = dplyr::case_when(
        startsWith(variable, "header_")     ~ label,
        row_type %in% c("level", "missing") ~ paste0("\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0", label),
        TRUE                                ~ paste0("\u00a0\u00a0\u00a0\u00a0", label)
      ))
  }) |>
  gtsummary::modify_table_styling(
    columns    = label,
    rows       = startsWith(variable, "header_"),
    text_format = "bold"
  ) |>
  gtsummary::modify_header(label ~ "") |>
  gtsummary::modify_footnote(all_stat_cols() ~ NA)

Table1_Placebo

# Save as Word
footer_text <- paste0(
  "Data are represented as Median (Q1, Q3); n / N (%).\n",
  "Abbreviations: ",
  "ACQ, Asthma Control Questionnaire; ",
  "BMI, Body Mass Index; ",
  "CRSsNP, Chronic Rhinosinusitis without Nasal Polyps; ",
  "CRSwNP, Chronic Rhinosinusitis with Nasal Polyps; ",
  "ETI, Endotracheal Intubation; ",
  "FeNO, Fractional exhaled Nitric Oxide; ",
  "FEV1, Forced Expiratory Volume in 1 second; ",
  "FVC, Forced Vital Capacity; ",
  "ICU, Intensive Care Unit; ",
  "mOCS, Maintenance Oral Corticosteroids; ",
  "ppb, parts per billion; ",
  "preBD, pre-Bronchodilator."
)

Table1_Placebo |>
  gtsummary::as_flex_table() |>
  flextable::fontsize(size = 9, part = "all") |>
  flextable::padding(padding = 2, part = "all") |>
  flextable::add_footer_lines(values = footer_text) |>
  flextable::fontsize(size = 8, part = "footer") |>
  flextable::set_table_properties(layout = "autofit", width = 1) |>
  flextable::save_as_docx(path = file.path(Figure_path, "Table1_Placebo.docx"))

#--------------------------------------------------------------------------------------------------------------------------------------
# 8.2.	Table 2 – Trial name, ethnicity and region across the 3 analysis datasets

.make_tbl2 <- function(data) {
  gtsummary::tbl_summary(
    data |>
      dplyr::mutate(
        Ethnicity = dplyr::recode(Ethnicity,
                                  "American_Indian_or_Alaska_Native"          = "American Indian or Alaska Native",
                                  "Black_or_African_American"                 = "Black or African American",
                                  "Native_Hawaiian_or_other_Pacific_Islander" = "Native Hawaiian or other Pacific Islander",
                                  "Multiple"                                  = "Multiple ethnicities"
        ),
        Region = dplyr::recode(Region,
                               "North_America" = "North America",
                               "South_America" = "South America",
                               "South_Africa"  = "South Africa"
        )
      ),
    include = c(Enrolled_Trial_name, Ethnicity, Region),
    statistic    = list(all_categorical() ~ "{n} / {N} ({p}%)"),
    missing_text = "(Missing)",
    label = list(
      Enrolled_Trial_name ~ "Trial name",
      Ethnicity           ~ "Ethnicity",
      Region              ~ "Region"
    )
  ) |>
    gtsummary::modify_table_body(function(tbl) {
      tbl |>
        dplyr::mutate(label = dplyr::case_when(
          row_type %in% c("level", "missing") ~
            paste0("\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0", label),
          TRUE ~ label
        ))
    })
}

tbl2_attack <- .make_tbl2(ORACLE_Placebo_Attack_Raw)
tbl2_lung   <- .make_tbl2(ORACLE_Placebo_Lung_Raw)
tbl2_acq    <- .make_tbl2(ORACLE_Placebo_ACQ_Raw)

Table2_Placebo <- gtsummary::tbl_merge(
  tbls        = list(tbl2_attack, tbl2_lung, tbl2_acq),
  tab_spanner = c(
    "**Asthma Exacerbations**",
    "**Lung Function**",
    "**Symptom Control**"
  )
) |>
  gtsummary::modify_header(label ~ "")

Table2_Placebo

Table2_Placebo |>
  gtsummary::as_flex_table() |>
  flextable::fontsize(size = 9, part = "all") |>
  flextable::padding(padding = 2, part = "all") |>
  flextable::save_as_docx(path = file.path(Figure_path, "Table2_Placebo_trials.docx"))

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 9.	Saving the datasets

save(ORACLE_Placebo_Attack_Raw,     file = file.path(Saving_path, "ORACLE_Placebo_Attack_Raw.RData"))
save(ORACLE_Placebo_Attack_imputed, file = file.path(Saving_path, "ORACLE_Placebo_Attack_imputed.RData"))
save(ORACLE_Placebo_Lung_Raw,       file = file.path(Saving_path, "ORACLE_Placebo_Lung_Raw.RData"))
save(ORACLE_Placebo_Lung_imputed,   file = file.path(Saving_path, "ORACLE_Placebo_Lung_imputed.RData"))
save(ORACLE_Placebo_ACQ_Raw,        file = file.path(Saving_path, "ORACLE_Placebo_ACQ_Raw.RData"))
save(ORACLE_Placebo_ACQ_imputed,    file = file.path(Saving_path, "ORACLE_Placebo_ACQ_imputed.RData"))

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

# FIGURE 1 - BOXPLOT

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 10.	Distribution of placebo responses across trials

# 10.1	Boxplot (SML's version)

# Trial order: grouped by severity (2, 1, 0), alphabetical within each group
.trial_severity <- dplyr::bind_rows(
  dplyr::select(ORACLE_Placebo_Attack_Raw, Enrolled_Trial_name, Trial_severe_attack_min),
  dplyr::select(ORACLE_Placebo_Lung_Raw,   Enrolled_Trial_name, Trial_severe_attack_min),
  dplyr::select(ORACLE_Placebo_ACQ_Raw,    Enrolled_Trial_name, Trial_severe_attack_min)
) |>
  dplyr::distinct(Enrolled_Trial_name, .keep_all = TRUE) |>
  dplyr::mutate(Enrolled_Trial_name = as.character(Enrolled_Trial_name)) |>
  dplyr::arrange(dplyr::desc(Trial_severe_attack_min), Enrolled_Trial_name)  # ← alphabetical ASC within group
.all_trials <- rev(.trial_severity$Enrolled_Trial_name)

# Shared fill scale: Trial_severe_attack_min (0 = none, 1 = >=1, 2 = >=2)
.severity_levels <- c("0", "1", "2")
.severity_fill <- ggplot2::scale_fill_manual(
  name   = "Asthma attack selection criteria",
  values = c("0" = "grey70", "1" = "goldenrod2", "2" = "darkred"),
  labels = c("0" = "No criteria",
             "1" = "\u2265 1 attack in the past 12 mo",
             "2" = "\u2265 2 attack in the past 12 mo"),
  drop   = FALSE
)

# Helper to build each panel
.make_panel <- function(data, x_var, x_label, xlim_range,
                        show_y = FALSE, show_legend = FALSE, label_fmt = "%.1f",
                        panel_title = NULL, pathway_black = FALSE) {  # ← new flag
  
  .n_trials       <- dplyr::n_distinct(data$Enrolled_Trial_name)
  .n_participants <- nrow(data)
  
  .missing_df <- data.frame(
    Enrolled_Trial_name = factor(
      setdiff(.all_trials, unique(data$Enrolled_Trial_name)),
      levels = .all_trials
    )
  )
  
  p <- data |>
    dplyr::mutate(
      Enrolled_Trial_name     = factor(Enrolled_Trial_name, levels = .all_trials),
      Trial_severe_attack_min = factor(Trial_severe_attack_min, levels = .severity_levels)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = {{ x_var }}, y = Enrolled_Trial_name,
                                 fill = Trial_severe_attack_min)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::stat_summary(fun = mean, geom = "point",
                          shape = 21, size = 2.5, fill = "white", color = "gray20") +
    
    # --- Bold mean labels: white for all trials (or all non-PATHWAY if flag is on) ---
    ggplot2::stat_summary(
      data = if (pathway_black) ~ dplyr::filter(., Enrolled_Trial_name != "PATHWAY") else ~ .,
      fun  = mean, geom = "text",
      ggplot2::aes(label = ggplot2::after_stat(sprintf(label_fmt, x))),
      size = 3.5, vjust = -0.7, color = "white", fontface = "bold"
    ) +
    
    ggplot2::geom_text(
      data = .missing_df,
      ggplot2::aes(x = 0, y = Enrolled_Trial_name, label = "MISSING"),
      inherit.aes = FALSE, color = "gray50", fontface = "italic", size = 3
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::coord_cartesian(xlim = xlim_range) +
    .severity_fill +
    ggplot2::labs(
      title    = panel_title,
      subtitle = paste0("N (trials) = ", .n_trials, ",  n (participants) = ", .n_participants),
      x        = x_label,
      y        = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y      = if (show_y) ggplot2::element_text(size = 12) else ggplot2::element_blank(),
      axis.ticks.y     = if (show_y) ggplot2::element_line()          else ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(size = 12),
      axis.title.x     = ggplot2::element_text(size = 13, face = "bold"),
      plot.title       = ggplot2::element_text(size = 13, face = "bold",  hjust = 0.5),  # ← centered
      plot.subtitle    = ggplot2::element_text(size = 11, color = "gray30", hjust = 0.5),  # ← centered
      legend.text      = ggplot2::element_text(size = 11),
      legend.title     = ggplot2::element_text(size = 12, face = "bold"),
      legend.position  = if (show_legend) "bottom" else "none"
    )
  
  # --- Conditionally add black bold label for PATHWAY only ---
  if (pathway_black) {
    p <- p + ggplot2::stat_summary(
      data = ~ dplyr::filter(., Enrolled_Trial_name == "PATHWAY"),
      fun  = mean, geom = "text",
      ggplot2::aes(label = ggplot2::after_stat(sprintf(label_fmt, x))),
      size = 3.5, vjust = -0.7, color = "black", fontface = "bold"
    )
  }
  
  p
}

# --- Panel A: pathway_black = TRUE ---
p_attack <- .make_panel(
  dplyr::filter(ORACLE_Placebo_Attack_Raw, !is.na(Placebo_response_Attack)),
  Placebo_response_Attack, "\u2206Asthma Attack Number", c(-6, 6),
  show_y = TRUE, show_legend = FALSE,
  panel_title = "A) Asthma Attacks Placebo Change",
  pathway_black = TRUE  # ← only here
)

# --- Panel B: pathway_black = FALSE (default) ---
p_lung <- .make_panel(
  dplyr::filter(ORACLE_Placebo_Lung_Raw, !is.na(Placebo_response_Lung)),
  Placebo_response_Lung, "\u2206FEV1 (mL)", c(-1000, 1000),
  show_y = FALSE, show_legend = TRUE, label_fmt = "%.0f",
  panel_title = "B) Lung Placebo Change"
)

# --- Panel C: pathway_black = FALSE (default) ---
p_acq <- .make_panel(
  dplyr::filter(ORACLE_Placebo_ACQ_Raw, !is.na(Placebo_response_ACQ)),
  Placebo_response_ACQ, "\u2206ACQ-5", c(-5, 5),
  show_y = FALSE, show_legend = FALSE,
  panel_title = "C) ACQ-5 Placebo Change"
)

# --- Combine and save ---
distribution_plot <- (p_attack | p_lung | p_acq)

distribution_plot

showtext::showtext_auto()
ggplot2::ggsave(
  filename = file.path(Figure_path, "Distribution_placebo_response.pdf"),
  plot     = distribution_plot,
  width    = 15,
  height   = 10
)
showtext::showtext_auto(enable = FALSE)


# 10.2	Violin plot

.make_panel_violin <- function(data, x_var, x_label, xlim_range,
                               show_y = FALSE, show_legend = FALSE, label_fmt = "%.1f") {
  .missing_df <- data.frame(
    Enrolled_Trial_name = factor(
      setdiff(.all_trials, unique(data$Enrolled_Trial_name)),
      levels = .all_trials
    )
  )
  
  data |>
    dplyr::mutate(
      Enrolled_Trial_name     = factor(Enrolled_Trial_name, levels = .all_trials),
      Trial_severe_attack_min = factor(Trial_severe_attack_min, levels = .severity_levels)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = {{ x_var }}, y = Enrolled_Trial_name,
                                 fill = Trial_severe_attack_min)) +
    ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
    ggplot2::stat_summary(fun = mean, geom = "point",
                          shape = 21, size = 2.5, fill = "white", color = "gray20") +
    ggplot2::stat_summary(fun = mean, geom = "text",
                          ggplot2::aes(label = ggplot2::after_stat(sprintf(label_fmt, x))),
                          size = 3.5, vjust = -0.7, color = "white") +
    ggplot2::geom_text(
      data = .missing_df,
      ggplot2::aes(x = 0, y = Enrolled_Trial_name, label = "MISSING"),
      inherit.aes = FALSE, color = "gray50", fontface = "italic", size = 3
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::coord_cartesian(xlim = xlim_range) +
    .severity_fill +
    ggplot2::labs(x = x_label, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y     = if (show_y) ggplot2::element_text() else ggplot2::element_blank(),
      axis.ticks.y    = if (show_y) ggplot2::element_line() else ggplot2::element_blank(),
      legend.position = if (show_legend) "bottom" else "none"
    )
}

# --- Panel A: Asthma attack rate ---
p_attack_v <- .make_panel_violin(
  dplyr::filter(ORACLE_Placebo_Attack_Raw, !is.na(Placebo_response_Attack)),
  Placebo_response_Attack, "\u2206Asthma Attack Number", c(-6, 6),
  show_y = TRUE, show_legend = FALSE
)

# --- Panel B: FEV1 ---
p_lung_v <- .make_panel_violin(
  dplyr::filter(ORACLE_Placebo_Lung_Raw, !is.na(Placebo_response_Lung)),
  Placebo_response_Lung, "\u2206FEV1 (mL)", c(-1000, 1000),
  show_y = FALSE, show_legend = TRUE, label_fmt = "%.0f"
)

# --- Panel C: ACQ-5 ---
p_acq_v <- .make_panel_violin(
  dplyr::filter(ORACLE_Placebo_ACQ_Raw, !is.na(Placebo_response_ACQ)),
  Placebo_response_ACQ, "\u2206ACQ-5", c(-5, 5),
  show_y = FALSE, show_legend = FALSE
)

# --- Combine and save ---
distribution_violin_plot <- (p_attack_v | p_lung_v | p_acq_v)

distribution_violin_plot

showtext::showtext_auto()
ggplot2::ggsave(
  filename = file.path(Figure_path, "Distribution_placebo_response_violin.pdf"),
  plot     = distribution_violin_plot,
  width    = 15,
  height   = 10
)
showtext::showtext_auto(enable = FALSE)


# 10.3	Density ridge plot

.make_panel_density <- function(data, x_var, x_label, xlim_range,
                                show_y = FALSE, show_legend = FALSE, label_fmt = "%.1f") {
  .missing_df <- data.frame(
    Enrolled_Trial_name = factor(
      setdiff(.all_trials, unique(data$Enrolled_Trial_name)),
      levels = .all_trials
    )
  )
  
  data |>
    dplyr::mutate(
      Enrolled_Trial_name     = factor(Enrolled_Trial_name, levels = .all_trials),
      Trial_severe_attack_min = factor(Trial_severe_attack_min, levels = .severity_levels)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = {{ x_var }}, y = Enrolled_Trial_name,
                                 fill = Trial_severe_attack_min)) +
    ggridges::geom_density_ridges(
      scale = 0.9, alpha = 0.8,
      rel_min_height = 0.01,
      quantile_lines = TRUE, quantiles = 2   # median line
    ) +
    ggplot2::stat_summary(fun = mean, geom = "point",
                          shape = 21, size = 2.5, fill = "white", color = "gray20") +
    ggplot2::stat_summary(fun = mean, geom = "text",
                          ggplot2::aes(label = ggplot2::after_stat(sprintf(label_fmt, x))),
                          size = 3.5, vjust = -0.7, color = "gray20") +
    ggplot2::geom_text(
      data = .missing_df,
      ggplot2::aes(x = 0, y = Enrolled_Trial_name, label = "MISSING"),
      inherit.aes = FALSE, color = "gray50", fontface = "italic", size = 3
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::coord_cartesian(xlim = xlim_range) +
    .severity_fill +
    ggplot2::labs(x = x_label, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y     = if (show_y) ggplot2::element_text() else ggplot2::element_blank(),
      axis.ticks.y    = if (show_y) ggplot2::element_line() else ggplot2::element_blank(),
      legend.position = if (show_legend) "bottom" else "none"
    )
}

# --- Panel A: Asthma attack rate ---
p_attack_d <- .make_panel_density(
  dplyr::filter(ORACLE_Placebo_Attack_Raw, !is.na(Placebo_response_Attack)),
  Placebo_response_Attack, "\u2206Asthma Attack Number", c(-6, 6),
  show_y = TRUE, show_legend = FALSE
)

# --- Panel B: FEV1 ---
p_lung_d <- .make_panel_density(
  dplyr::filter(ORACLE_Placebo_Lung_Raw, !is.na(Placebo_response_Lung)),
  Placebo_response_Lung, "\u2206FEV1 (mL)", c(-1000, 1000),
  show_y = FALSE, show_legend = TRUE, label_fmt = "%.0f"
)

# --- Panel C: ACQ-5 ---
p_acq_d <- .make_panel_density(
  dplyr::filter(ORACLE_Placebo_ACQ_Raw, !is.na(Placebo_response_ACQ)),
  Placebo_response_ACQ, "\u2206ACQ-5", c(-5, 5),
  show_y = FALSE, show_legend = FALSE
)

# --- Combine and save ---
distribution_density_plot <- (p_attack_d | p_lung_d | p_acq_d)

distribution_density_plot

showtext::showtext_auto()
ggplot2::ggsave(
  filename = file.path(Figure_path, "Distribution_placebo_response_density.pdf"),
  plot     = distribution_density_plot,
  width    = 15,
  height   = 10
)
showtext::showtext_auto(enable = FALSE)

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
## BOXPLOT (JSP's version)

Data_Oracle_attack_clean <- ORACLE_Placebo_Attack_Raw

sum(is.na(Data_Oracle_attack_clean$Placebo_response_Attack))
      
## DELTA ASTHMA ATTACK - BOXPLOT
library(dplyr)

# Define the list of trials with missing data
missing_trials_attack <- c("LUTE", "MILLY", "VERSE")

# Define the desired order
ordered_trials_attack <- c(
  "VERSE","MILLY","LUTE","LAVOLTA_2", "LAVOLTA_1",
  "QUEST", "EXTRA", "DRI12544", "STRATOS_2", "STRATOS_1",
  "PATHWAY", "NAVIGATOR", "LUSTER_2", "LUSTER_1",
  "DREAM", "BENRAP2B", "AZISAST"
)

# Define colors based on criteria
trial_colors_attack <- c(
  "MISSING"                = "black",
  "No criteria"            = "#A0A0A0",
  "With ≥1 asthma attack"  = "#FFA000",
  "With ≥2 asthma attacks" = "#D84315"
)

# Create a mapping of trials to criteria
trial_criteria_map_attack <- data.frame(
  Enrolled_Trial_name = c(
    "VERSE", "MILLY","LUTE","LAVOLTA_2", "LAVOLTA_1",
    "QUEST", "EXTRA", "DRI12544", "STRATOS_2", "STRATOS_1",
    "PATHWAY", "NAVIGATOR", "LUSTER_2", "LUSTER_1",
    "DREAM", "BENRAP2B", "AZISAST"
  ),
  Attack_History_Criteria = c(
    "MISSING", "MISSING", "MISSING","No criteria", "No criteria",
    "With ≥1 asthma attack", "With ≥1 asthma attack", "With ≥1 asthma attack",
    "With ≥2 asthma attacks", "With ≥2 asthma attacks", "With ≥2 asthma attacks",
    "With ≥2 asthma attacks", "With ≥2 asthma attacks", "With ≥2 asthma attacks",
    "With ≥2 asthma attacks", "With ≥2 asthma attacks", "With ≥2 asthma attacks"
  ),
  stringsAsFactors = FALSE
)

# Join criteria to plot data
Data_Oracle_attack_plot <- Data_Oracle_attack_clean %>%
  dplyr::left_join(trial_criteria_map_attack, by = "Enrolled_Trial_name") %>%
  dplyr::mutate(
    Enrolled_Trial_name = factor(
      Enrolled_Trial_name,
      levels = ordered_trials_attack
    ),
    Attack_History_Criteria = factor(
      Attack_History_Criteria,
      levels = c(
        "No criteria",
        "With ≥1 asthma attack",
        "With ≥2 asthma attacks",
        "MISSING"
      )
    )
  )

# Separate missing trials data for label
missing_label_data <- trial_criteria_map_attack %>%
  dplyr::filter(Enrolled_Trial_name %in% missing_trials_attack) %>%
  dplyr::mutate(
    Enrolled_Trial_name = factor(Enrolled_Trial_name, levels = ordered_trials_attack)
  )

# Non-missing data shortcut
Data_Oracle_attack_plot_nonmissing <- Data_Oracle_attack_plot %>%
  dplyr::filter(!Enrolled_Trial_name %in% missing_trials_attack)

# Plot
Boxplot_Placebo_attack_per_trial <- ggplot() +
  
  # 1. Min-to-max line (bottom layer)
  stat_summary(
    data = Data_Oracle_attack_plot_nonmissing,
    aes(x = Enrolled_Trial_name, y = Placebo_response_Attack),
    fun.min  = min,
    fun.max  = max,
    geom     = "linerange",
    linewidth = 0.6,
    color    = "black"
  ) +
  
  # 2. Median crossbar for zero-IQR trials
  stat_summary(
    data = Data_Oracle_attack_plot_nonmissing,
    aes(x = Enrolled_Trial_name, y = Placebo_response_Attack),
    fun  = median,
    geom = "crossbar",
    width     = 0.8,
    linewidth = 0.4,
    color    = "black",
    fill     = NA
  ) +
  
  # 3. Boxplot on top of lines
  geom_boxplot(
    data = Data_Oracle_attack_plot_nonmissing,
    aes(
      x    = Enrolled_Trial_name,
      y    = Placebo_response_Attack,
      fill = Attack_History_Criteria
    ),
    width          = 0.8,
    linewidth      = 0.6,
    outlier.colour = NA,
    coef           = 100,
    alpha          = 0.8,
    color          = "black"
  ) +
  
  # 4. Mean point (on top of box)
  stat_summary(
    data = Data_Oracle_attack_plot_nonmissing,
    aes(x = Enrolled_Trial_name, y = Placebo_response_Attack),
    fun    = mean,
    geom   = "point",
    shape  = 23,
    size   = 1,
    stroke = 0.4,
    fill   = "white",
    colour = "black"
  ) +
  
  # 5. Mean label (topmost layer)
  stat_summary(
    data = Data_Oracle_attack_plot_nonmissing,
    aes(x = Enrolled_Trial_name, y = Placebo_response_Attack),
    fun.data = function(x) data.frame(
      y     = mean(x, na.rm = TRUE),
      label = round(mean(x, na.rm = TRUE), 1)
    ),
    geom     = "text",
    size     = 4,
    color    = "black",
    fontface = "bold",
    vjust    = -0.5
  ) +
  
  geom_text(
    data = missing_label_data,
    aes(x = Enrolled_Trial_name, y = 0, label = "MISSING"),
    size     = 4,
    color    = "black",
    fontface = "italic"
  ) +
  
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
  
  scale_fill_manual(
    values = trial_colors_attack,
    name   = "Attack History Criteria",
    drop   = FALSE
  ) +
  
  scale_x_discrete(limits = ordered_trials_attack, drop = FALSE) +
  
  scale_y_continuous(
    limits = c(-5, 5),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  
  labs(
    title    = element_blank(),
    subtitle = "N = 14   ;   n = 4381",
    x        = NULL,
    y        = "Delta Asthma Attacks"
  ) +
  
  theme_minimal() +
  theme(
    legend.position      = "bottom",
    legend.justification = "left",
    legend.box           = "horizontal",
    legend.title         = element_text(size = 12, face = "bold"),
    legend.text          = element_text(size = 12),
    axis.line.x          = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x         = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length.x  = unit(3, "pt"),
    axis.title.x         = element_text(size = 16, face = "bold"),
    axis.title.y         = element_text(size = 14, face = "bold"),
    axis.text.x          = element_text(face = "bold", size = 12, margin = margin(t = 12)),
    axis.text.y          = element_text(size = 12, face = "bold"),
    panel.grid.major.y   = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.major.x   = element_blank(),
    panel.background     = element_rect(fill = "white", color = "white"),
    plot.background      = element_rect(fill = "white", color = "white")
  ) +
  
  coord_flip()

Boxplot_Placebo_attack_per_trial

pdf("/Users/joel/ORACLE - Placebo - Exacerbations/Project/Boxplot_Placebo_Attacks_per_trial.pdf",
    width = 15, height = 12)

print(Boxplot_Placebo_attack_per_trial)

dev.off()

# Save the plot
ggsave("/Users/joel/ORACLE - Placebo - Exacerbations/Project/Plot_attack_placebo_per_trial.png",
       plot = Boxplot_Placebo_attack_per_trial, width=10, height =10, dpi = 300)

#======================================================================================================================================
Data_Oracle_lung_clean<- ORACLE_Placebo_Lung_Raw

sum(is.na(Data_Oracle_lung_clean$Placebo_response_Lung))

##DELTA FEV1 - BOXPLOT
Studies_included_lung
# Define the list of trials with missing data
missing_trials_lung <- c("VERSE","MILLY", "LUTE", "EXTRA", "DRI12544",
                         "LUSTER_2", "LUSTER_1", "DREAM", "AZISAST")

# Define the desired order (for use only in the plot)
ordered_trials <- c("VERSE","MILLY", "LUTE", "LAVOLTA_2", "LAVOLTA_1", 
                    "QUEST", "EXTRA", "DRI12544", "STRATOS_2", "STRATOS_1", 
                    "PATHWAY", "NAVIGATOR", "LUSTER_2", "LUSTER_1", 
                    "DREAM", "BENRAP2B", "AZISAST")

# Define the color gradient
trial_colors <- c(
  "AZISAST"    = "#FFB6C1",  # LightPink
  "BENRAP2B"   = "#FF69B4",  # HotPink
  "DREAM"      = "#FF1493",  # DeepPink
  "LUSTER_1"   = "#BA55D3",  # MediumOrchid (lighter than #9932CC)
  "LUSTER_2"   = "#9370DB",  # MediumPurple (more readable than #8B008B)
  "NAVIGATOR"  = "#3399FF",  # Lighter blue (as requested)
  "PATHWAY"    = "#87CEEB",  # SkyBlue
  "STRATOS_1"  = "#40E0D0",  # Turquoise
  "STRATOS_2"  = "#AFEEEE",  # PaleTurquoise
  "DRI12544"   = "#90EE90",  # LightGreen
  "QUEST"      = "#66CDAA",  # MediumAquamarine (soft green)
  "EXTRA"      = "#F0E68C",  # Khaki
  "LAVOLTA_1"  = "#ADFF2F",  # GreenYellow
  "LAVOLTA_2"  = "#FFD700",  # Gold (brighter than DarkGoldenrod)
  "LUTE"       = "#FFA500",  # Orange (better for text contrast)
  "VERSE"      = "#FF6347"   # Tomato (lighter than #FF4500)
)

##FOR ERS
trial_colors <- c(
  "VERSE"      = "#A0A0A0",   # No criteria
  "MILLY"       = "#A0A0A0",   # No criteria
  "LUTE"       = "#A0A0A0",   # No criteria
  
  "LAVOLTA_2"  = "#A0A0A0",   # No criteria
  "LAVOLTA_1"  = "#A0A0A0",   # No criteria
  
  "QUEST"      = "#FFA000",   # ≥1 attack
  "EXTRA"      = "#FFA000",   # ≥1 attack
  "DRI12544"   = "#FFA000",   # ≥1 attack
  
  "STRATOS_2"  = "#D84315",   # ≥2 attacks
  "STRATOS_1"  = "#D84315",   # ≥2 attacks
  "PATHWAY"    = "#D84315",   # ≥2 attacks
  "NAVIGATOR"  = "#D84315",   # ≥2 attacks
  "LUSTER_2"   = "#D84315",   # ≥2 attacks
  "LUSTER_1"   = "#D84315",   # ≥2 attacks
  "DREAM"      = "#D84315",   # ≥2 attacks
  "BENRAP2B"   = "#D84315",   # ≥2 attacks
  "AZISAST"    = "#D84315"    # ≥2 attacks
)

# Plot
Boxplot_Placebo_lung_per_trial <- ggplot() +
  
  # Boxplots for non-missing trials
  geom_boxplot(
    data = Data_Oracle_lung_clean_imp0 %>% filter(!Enrolled_Trial_name %in% missing_trials_lung),
    aes(x = Enrolled_Trial_name, y = Placebo_response_Lung, fill = Enrolled_Trial_name),
    width = 0.8,
    linewidth = 0.4,
    outlier.colour = NA,
    alpha = 0.8,
    color = "black"
  ) +
  
  # Mean diamond
  stat_summary(
    data = Data_Oracle_lung %>% filter(!Enrolled_Trial_name %in% missing_trials_lung),
    aes(x = Enrolled_Trial_name, y = Placebo_response_Lung),
    fun = mean,
    geom = "point",
    shape = 23,
    size = 1,
    stroke = 0.4,
    fill = "white",
    colour = "black"
  ) +
  
  # Mean value labels
  stat_summary(
    data = Data_Oracle_lung %>% filter(!Enrolled_Trial_name %in% missing_trials_lung),
    aes(x = Enrolled_Trial_name, y = Placebo_response_Lung),
    fun.data = \(x) data.frame(
      y = mean(x, na.rm = TRUE),
      label = round(mean(x, na.rm = TRUE), 1)
    ),
    geom = "text",
    size = 4,
    color = "black",
    fontface = "bold",
    vjust = -0.5 # moves label above point
  ) +
  
  # Add "MISSING" labels
  geom_text(
    data = data.frame(Enrolled_Trial_name = missing_trials_lung),
    aes(x = Enrolled_Trial_name, y = 0, label = "MISSING"),
    size = 4,
    color = "black",
    fontface = "italic"
  ) +
  
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
  
  scale_fill_manual(values = trial_colors, guide = "none") +
  
  # Force the x-axis order here only (no mutate before)
  scale_x_discrete(limits = ordered_trials, drop = FALSE) +
  scale_y_continuous(limits = c(-1000, 1000), expand = expansion(mult = c(0.05, 0.1))) +
  
  labs(
    title = "Lung Placebo Change",
    subtitle = "N = 8   ;   n = 2 590",
    x = NULL,
    y = "∆ FEV1 (ml)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 20, face = "bold", hjust = .5),
    axis.line.x          = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x         = element_line(color = "black", linewidth = 0.6),
    plot.subtitle = element_text(size = 11, hjust = .5),
    axis.title.x  = element_text(size = 16, face = "bold"),
    axis.title.y  = element_text(size = 14, face = "bold"),
    axis.text.x   = element_text(size = 12),
    axis.text.y   = element_blank(),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey90", linewidth = .5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "white", color = "white"),
    plot.background    = element_rect(fill = "white", color = "white")
  ) +
  coord_flip()


Boxplot_Placebo_lung_per_trial

pdf("/Users/joel/ORACLE - Placebo - Exacerbations/Project/Boxplot_Placebo_lung_per_trial.pdf",
    width = 15, height = 12)

print(Boxplot_Placebo_lung_per_trial)

dev.off()

#======================================================================================================================================
Data_Oracle_ACQ_clean<- ORACLE_Placebo_ACQ_Raw

##DELTA ACQ - BOXPLOT
Studies_included_ACQ
# Define the list of trials with missing data
missing_trials_acq <- c("QUEST", 
                        "EXTRA", "DRI12544", "PATHWAY",
                        "NAVIGATOR", "LUSTER_2", "LUSTER_1", "DREAM", "AZISAST")

# Define the desired order (for use only in the plot)
ordered_trials <- c("VERSE","MILLY", "LUTE", "LAVOLTA_2", "LAVOLTA_1", "QUEST", 
                    "EXTRA", "DRI12544", "STRATOS_2", "STRATOS_1", "PATHWAY",
                    "NAVIGATOR", "LUSTER_2", "LUSTER_1", "DREAM", "BENRAP2B", "AZISAST")

# Define the color gradient
trial_colors <- c(
  "AZISAST"    = "#FFB6C1",  # LightPink
  "BENRAP2B"   = "#FF69B4",  # HotPink
  "DREAM"      = "#FF1493",  # DeepPink
  "LUSTER_1"   = "#BA55D3",  # MediumOrchid (lighter than #9932CC)
  "LUSTER_2"   = "#9370DB",  # MediumPurple (more readable than #8B008B)
  "NAVIGATOR"  = "#3399FF",  # Lighter blue (as requested)
  "PATHWAY"    = "#87CEEB",  # SkyBlue
  "STRATOS_1"  = "#40E0D0",  # Turquoise
  "STRATOS_2"  = "#AFEEEE",  # PaleTurquoise
  "DRI12544"   = "#90EE90",  # LightGreen
  "QUEST"      = "#66CDAA",  # MediumAquamarine (soft green)
  "EXTRA"      = "#F0E68C",  # Khaki
  "LAVOLTA_1"  = "#ADFF2F",  # GreenYellow
  "LAVOLTA_2"  = "#FFD700",  # Gold (brighter than DarkGoldenrod)
  "LUTE"       = "#FFA500",  # Orange (better for text contrast)
  "VERSE"      = "#FF6347"   # Tomato (lighter than #FF4500)
)

trial_colors <- c(
  "VERSE"      = "#A0A0A0",   # No criteria
  "MILLY"       = "#A0A0A0",   # No criteria
  "LUTE"       = "#A0A0A0",   # No criteria
  
  "LAVOLTA_2"  = "#A0A0A0",   # No criteria
  "LAVOLTA_1"  = "#A0A0A0",   # No criteria
  
  "QUEST"      = "#FFA000",   # ≥1 attack
  "EXTRA"      = "#FFA000",   # ≥1 attack
  "DRI12544"   = "#FFA000",   # ≥1 attack
  
  "STRATOS_2"  = "#D84315",   # ≥2 attacks
  "STRATOS_1"  = "#D84315",   # ≥2 attacks
  "PATHWAY"    = "#D84315",   # ≥2 attacks
  "NAVIGATOR"  = "#D84315",   # ≥2 attacks
  "LUSTER_2"   = "#D84315",   # ≥2 attacks
  "LUSTER_1"   = "#D84315",   # ≥2 attacks
  "DREAM"      = "#D84315",   # ≥2 attacks
  "BENRAP2B"   = "#D84315",   # ≥2 attacks
  "AZISAST"    = "#D84315"    # ≥2 attacks
)


# Plot
Boxplot_Placebo_ACQ_per_trial <- ggplot() +
  
  # Boxplots for non-missing trials
  geom_boxplot(
    data = Data_Oracle_ACQ_clean_imp0 %>% filter(!Enrolled_Trial_name %in% missing_trials_acq),
    aes(x = Enrolled_Trial_name, y = Placebo_response_ACQ, fill = Enrolled_Trial_name),
    width = 0.8,
    linewidth = 0.4,
    outlier.colour = NA,
    alpha = 0.8,
    color = "black"
  ) +
  
  # Mean diamond
  stat_summary(
    data = Data_Oracle_ACQ %>% filter(!Enrolled_Trial_name %in% missing_trials_acq),
    aes(x = Enrolled_Trial_name, y = Placebo_response_ACQ),
    fun = mean,
    geom = "point",
    shape = 23,
    size = 1,
    stroke = 0.4,
    fill = "white",
    colour = "black"
  ) +
  
  # Mean value labels
  stat_summary(
    data = Data_Oracle_ACQ %>% filter(!Enrolled_Trial_name %in% missing_trials_acq),
    aes(x = Enrolled_Trial_name, y = Placebo_response_ACQ),
    fun.data = \(x) data.frame(
      y = mean(x, na.rm = TRUE),
      label = round(mean(x, na.rm = TRUE), 1)
    ),
    geom = "text",
    size = 4,
    color = "black",
    fontface = "bold",
    vjust = -0.5 # moves label above point
  ) +
  
  # Add "MISSING" labels
  geom_text(
    data = data.frame(Enrolled_Trial_name = missing_trials_acq),
    aes(x = Enrolled_Trial_name, y = 0, label = "MISSING"),
    size = 4,
    color = "black",
    fontface = "italic"
  ) +
  
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
  
  scale_fill_manual(values = trial_colors, guide = "none") +
  
  # Force the x-axis order here only (no mutate before)
  scale_x_discrete(limits = ordered_trials, drop = FALSE) +
  scale_y_continuous(limits = c(-5, 5), expand = expansion(mult = c(0.05, 0.1))) +
  
  labs(
    title = "ACQ-5 Placebo Change",
    subtitle = "N = 8   ;   n = 1736",
    x = NULL,
    y = "∆ ACQ-5"
  ) +
  
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 20, face = "bold", hjust = .5),
    axis.line.x          = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x         = element_line(color = "black", linewidth = 0.6),
    plot.subtitle = element_text(size = 11, hjust = .5),
    axis.title.x  = element_text(size = 16, face = "bold"),
    axis.title.y  = element_text(size = 14, face = "bold"),
    axis.text.x   = element_text(size = 12),
    axis.text.y   = element_blank() ,
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey90", linewidth = .5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "white", color = "white"),
    plot.background    = element_rect(fill = "white", color = "white")
  ) +
  coord_flip()

Boxplot_Placebo_ACQ_per_trial

pdf("/Users/joel/ORACLE - Placebo - Exacerbations/Project/Boxplot_Placebo_ACQ_per_trial.pdf",
    width = 15, height = 12)

print(Boxplot_Placebo_ACQ_per_trial)

dev.off()

#======================================================================================================================================

##MERGE - 3 boxplots

# Combine the three plots
library(patchwork)

Boxplot_Placebo_3Outcomes_per_trial <- 
  (Boxplot_Placebo_attack_per_trial | 
     Boxplot_Placebo_lung_per_trial | 
     Boxplot_Placebo_ACQ_per_trial) +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Placebo Change Across Trials in Three Outcomes",
    theme = theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5))
  )

Boxplot_Placebo_3Outcomes_per_trial

pdf("/Users/joel/ORACLE - Placebo - Exacerbations/Project/Boxplot_Placebo_3Outcomes_per_trial.pdf",
    width = 15, height = 12)

print(Boxplot_Placebo_3Outcomes_per_trial)

dev.off()

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

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

# Usually, n of detectCores minus 2
plan(multisession, workers = 6)

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
  list(pred = "FEV1_preBD_L_0W_per_10",        slope = "FEV1_preBD_L_0W"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_log10",                       slope = "BEC_log10"),
  list(pred = "FeNO_log10",                      slope = "FeNO_log10"),
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
    followup_offset = "Yes",
    followup_col    = "Follow_up_duration_days",
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
      "FEV1_preBD_L_0W_per_10"           ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"             ~ "FEV1 reversibility (per 10% increase)",
      "BEC_log10"                          ~ "Blood eosinophils (log10)",
      "FeNO_log10"                         ~ "FeNO (log10)",
      "IgE_log10"                          ~ "Total IgE (log10)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5")
      ~ "Demographic",
      c("Smoking_history_yes_noYes", "Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline lung function",
      c("BEC_log10", "FeNO_log10", "IgE_log10")
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
Figure_Placebo_Attack_simple_forest

ggsave(
  filename = file.path(Figure_path, "Figure_Placebo_Attack_simple_forest.pdf"),
  plot     = Figure_Placebo_Attack_simple_forest,
  
  width    = 8, height = 8, units = "in"
)
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
print(cov_table, width = Inf)
View(cov_table)

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
      "FEV1_preBD_L_0W_per_10"           ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"       ~ "FEV1 reversibility (per 10% increase)",
      "BEC_log10"                          ~ "Blood eosinophils (log10)",
      "FeNO_log10"                         ~ "FeNO (log10)",
      "IgE_log10"                          ~ "Total IgE (log10)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5")
      ~ "Demographic",
      c("Smoking_history_yes_noYes", "Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline lung function",
      c("BEC_log10", "FeNO_log10", "IgE_log10")
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
  filename = file.path(Figure_path, "Figure_Placebo_Attack_final_forest.pdf"),
  plot     = Figure_Placebo_Attack_final_forest,
  width    = 10, height = 6, units = "in"
)

#######################################################################################################################################
# 4.5 saving the results from multrivariable model
#--------------------------------------------------------------------------------------------------------------------------------------
save(Coefficients_Placebo_Attack_final_model,
     file = file.path(Figure_path, "Coefficients_Placebo_Attack_final_model.RData"))


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
  list(pred = "Age_per_10",                     slope = "Age"),
  list(pred = "Sex",                             slope = "Sex"),
  list(pred = "BMI_per_5",                       slope = "BMI"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no"),
  list(pred = "Airborne_allergen_sensitisation",  slope = "Airborne_allergen_sensitisation"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis"),
  list(pred = "Eczema",                          slope = "Eczema"),
  list(pred = "CRSwNP",                          slope = "CRSwNP"),
  list(pred = "CRSsNP",                          slope = "CRSsNP"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num"),  # no covariable for this one
  list(pred = "FEV1_preBD_L_0W_per_10",        slope = "FEV1_preBD_L_0W"),
  list(pred = "FEV1_reversibility_0W_per_10",           slope = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_log10",                       slope = "BEC_log10"),
  list(pred = "FeNO_log10",                      slope = "FeNO_log10"),
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
      "FEV1_preBD_L_0W_per_10"           ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"             ~ "FEV1 reversibility (per 10% increase)",
      "BEC_log10"                          ~ "Blood eosinophils (log10)",
      "FeNO_log10"                         ~ "FeNO (log10)",
      "IgE_log10"                          ~ "Total IgE (log10)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5")
      ~ "Demographic",
      c("Smoking_history_yes_noYes", "Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline lung function",
      c("BEC_log10", "FeNO_log10", "IgE_log10")
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
.x_pval <- 540  # data-unit x for p-value column (outside xlim of -160 to 160)
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
  scale_x_continuous(breaks = seq(-400, 400, by = 40)) +
  scale_y_discrete(labels = .axis_labels) +
  coord_cartesian(xlim = c(-400, 400), clip = "off") +
  labs(
    x        = expression(Delta ~ "FEV (mL)"),
    y        = NULL,
    title    = "Predictors of placebo response \u2014 Lung function",
    subtitle = "Simple models (each predictor adjusted for FEV1 pre-BD (L))"
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
Figure_Placebo_Lung_simple_forest <- (.p_forest / .p_arrows) + plot_layout(heights = c(20, 1))
Figure_Placebo_Lung_simple_forest

ggsave(
  filename = file.path(Figure_path, "Figure_Placebo_Lung_simple_forest.pdf"),
  plot     = Figure_Placebo_Lung_simple_forest,
  
  width    = 10, height = 6, units = "in"
)

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
  list(pred = "FEV1_preBD_L_0W_per_10",        slope = "FEV1_preBD_L_0W",               term = "FEV1_preBD_L_0W_per_10"),
  list(pred = "FEV1_reversibility_0W_per_10",    slope = "FEV1_reversibility_0W_per_10",    term = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_log10",                       slope = "BEC_log10",                       term = "BEC_log10"),
  list(pred = "FeNO_log10",                      slope = "FeNO_log10",                      term = "FeNO_log10"),
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

print(cov_table, width = Inf)

View(cov_table)

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
      "FEV1_preBD_L_0W_per_10"           ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"       ~ "FEV1 reversibility (per 10% increase)",
      "BEC_log10"                          ~ "Blood eosinophils (log10)",
      "FeNO_log10"                         ~ "FeNO (log10)",
      "IgE_log10"                          ~ "Total IgE (log10)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5")
      ~ "Demographic",
      c("Smoking_history_yes_noYes", "Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline lung function",
      c("BEC_log10", "FeNO_log10", "IgE_log10")
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
  filename = file.path(Figure_path, "Figure_Placebo_Lung_final_forest.pdf"),
  plot     = Figure_Placebo_Lung_final_forest,
  width    = 10, height = 6, units = "in"
)

#######################################################################################################################################
# 4.5 saving the results from multrivariable model
#--------------------------------------------------------------------------------------------------------------------------------------
save(Coefficients_Placebo_Lung_final_model,
     file = file.path(Figure_path, "Coefficients_Placebo_Lung_final_model.RData"))


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
  list(pred = "Age_per_10",                     slope = "Age"),
  list(pred = "Sex",                             slope = "Sex"),
  list(pred = "BMI_per_5",                       slope = "BMI"),
  list(pred = "Smoking_history_yes_no",          slope = "Smoking_history_yes_no"),
  list(pred = "Airborne_allergen_sensitisation",  slope = "Airborne_allergen_sensitisation"),
  list(pred = "Allergic_Rhinitis",               slope = "Allergic_Rhinitis"),
  list(pred = "Eczema",                          slope = "Eczema"),
  list(pred = "CRSwNP",                          slope = "CRSwNP"),
  list(pred = "CRSsNP",                          slope = "CRSsNP"),
  list(pred = "Psychiatric_disease",             slope = "Psychiatric_disease"),
  list(pred = "GINA_step_numeric",               slope = "GINA_step_numeric"),
  list(pred = "ACQ_score_0W",                    slope = "ACQ_score_0W"),
  list(pred = "Attack_history_num",              slope = "Attack_history_num"),  # no covariable for this one
  list(pred = "FEV1_preBD_L_0W_per_10",        slope = "FEV1_preBD_L_0W"),
  list(pred = "FEV1_reversibility_0W_per_10",           slope = "FEV1_reversibility_0W_per_10"),
  list(pred = "BEC_log10",                       slope = "BEC_log10"),
  list(pred = "FeNO_log10",                      slope = "FeNO_log10"),
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
      "Allergic_RhinitisYes"              ~ "Allergic rhinitis",
      "EczemaYes"                          ~ "Eczema",
      "CRSwNPYes"                          ~ "CRS with nasal polyps",
      "CRSsNPYes"                          ~ "CRS without nasal polyps",
      "Psychiatric_diseaseYes"             ~ "Psychiatric disease",
      "GINA_step_numeric"                  ~ "GINA step",
      "ACQ_score_0W"                       ~ "ACQ score (baseline)",
      "Attack_history_num"                 ~ "Exacerbation history (count)",
      "FEV1_preBD_L_0W_per_10"           ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"             ~ "FEV1 reversibility (per 10% increase)",
      "BEC_log10"                          ~ "Blood eosinophils (log10)",
      "FeNO_log10"                         ~ "FeNO (log10)",
      "IgE_log10"                          ~ "Total IgE (log10)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5")
      ~ "Demographic",
      c("Smoking_history_yes_noYes", "Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline ACQ function",
      c("BEC_log10", "FeNO_log10", "IgE_log10")
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
.x_pval <- 0.950   # data-unit x for p-value column (outside xlim of -0.5 to 0.5)
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
  scale_x_continuous(breaks = seq(-0.75, 0.75, by = 0.25)) +
  scale_y_discrete(labels = .axis_labels) +
  coord_cartesian(xlim = c(-0.75, 0.75), clip = "off") +
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
Figure_Placebo_ACQ_simple_forest

ggsave(
  filename = file.path(Figure_path, "Figure_Placebo_ACQ_simple_forest.pdf"),
  plot     = Figure_Placebo_ACQ_simple_forest,
  
  width    = 8, height = 8, units = "in"
)

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

print(cov_table, width = Inf)

View(cov_table)

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
      "Allergic_RhinitisYes"              ~ "Allergic rhinitis",
      "EczemaYes"                          ~ "Eczema",
      "CRSwNPYes"                          ~ "CRS with nasal polyps",
      "CRSsNPYes"                          ~ "CRS without nasal polyps",
      "Psychiatric_diseaseYes"             ~ "Psychiatric disease",
      "GINA_step_numeric"                  ~ "GINA step",
      "ACQ_score_0W"                       ~ "ACQ score (baseline)",
      "Attack_history_num"                 ~ "Exacerbation history (count)",
      "FEV1_preBD_L_0W_per_10"           ~ "FEV1 pre-BD (L., per 0.1 decrease)",
      "FEV1_reversibility_0W_per_10"       ~ "FEV1 reversibility (per 10% increase)",
      "BEC_log10"                          ~ "Blood eosinophils (log10)",
      "FeNO_log10"                         ~ "FeNO (log10)",
      "IgE_log10"                          ~ "Total IgE (log10)"
    ),
    category = case_match(
      term,
      c("Age_per_10", "SexMale", "BMI_per_5")
      ~ "Demographic",
      c("Smoking_history_yes_noYes", "Airborne_allergen_sensitisationYes",
        "Allergic_RhinitisYes", "EczemaYes", "CRSwNPYes", "CRSsNPYes",
        "Psychiatric_diseaseYes")
      ~ "Comorbidities",
      c("GINA_step_numeric", "ACQ_score_0W", "Attack_history_num")
      ~ "Asthma history",
      c("FEV1_preBD_L_0W_per_10", "FEV1_reversibility_0W_per_10")
      ~ "Baseline ACQ function",
      c("BEC_log10", "FeNO_log10", "IgE_log10")
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
  filename = file.path(Figure_path, "Figure_Placebo_ACQ_final_forest.pdf"),
  plot     = Figure_Placebo_ACQ_final_forest,
  width    = 10, height = 6, units = "in"
)
#######################################################################################################################################
# 4.5 saving the results from multrivariable model
#--------------------------------------------------------------------------------------------------------------------------------------
save(Coefficients_Placebo_ACQ_final_model,
     file = file.path(Saving_path, "Coefficients_Placebo_ACQ_final_model.RData"))

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

Coefficients_Placebo_Attack_final_model
Coefficients_Placebo_Lung_final_model
Coefficients_Placebo_ACQ_final_model

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 3.  Multi-panel forest plot: predictors of placebo response across 3 models
#--------------------------------------------------------------------------------------------------------------------------------------

# 3.1  Y-axis structure with category headers and spacers (top → bottom reading order)
#      Spacers use different numbers of spaces to keep factor levels unique.
#      Category headers use **bold** markdown syntax (rendered via ggtext).
.y_structure <- tibble::tribble(
  ~y_display,                                                       ~term,
  "**Demographic**",                                               NA_character_,
  "Age (per 10 years)",                                             "Age_per_10",
  "Sex (Male vs Female)",                                           "SexMale",
  "Body Mass Index (per 5 kg/m\u00b2)",                             "BMI_per_5",
  "History of smoking (current/former vs never)",                                               "Smoking_history_yes_noYes",
  "",                                                               NA_character_,
  "**Comorbidities**",                                              NA_character_,
  "Airborne allergen sensitisation*",                               "Airborne_allergen_sensitisationYes",
  "Allergic Rhinitis*",                                             "Allergic_RhinitisYes",
  "Eczema*",                                                        "EczemaYes",
  "CRS without nasal polyps*",                                      "CRSsNPYes",
  "CRS with nasal polyps*",                                         "CRSwNPYes",
  "Psychiatric disease*",                                           "Psychiatric_diseaseYes",
  " ",                                                              NA_character_,
  "**Asthma history**",                                             NA_character_,
  "Treatment step (per step increase)",                             "GINA_step_numeric",
  "Baseline ACQ-5 score",                                           "ACQ_score_0W",
  "Number of asthma attack in the past 12 mo",                      "Attack_history_num",
  "  ",                                                             NA_character_,
  "**Baseline lung function**",                                     NA_character_,
  "FEV1 (L., per 0.1 decrease)",                                       "FEV1_preBD_L_0W_per_10",
  "FEV1% Reversibility to BD (per 10% increase)",                   "FEV1_reversibility_0W_per_10",
  "   ",                                                            NA_character_,
  "**Inflammatory biomarkers**",                                    NA_character_,
  "BEC (per 10 fold increase)",                 "BEC_log10",
  "FeNO (per 10 fold increase)",         "FeNO_log10",
  "Total IgE (per 10 fold increase)",                  "IgE_log10"
)

# rev() so the first row of the tibble (Demographics) appears at the top of the plot
.y_levels <- rev(.y_structure$y_display)

# 3.2  Prep helper: join structure onto model coefficients, flag significance
.prep_model_v2 <- function(df) {
  df |>
    dplyr::rename(ci_low = `2.5 %`, ci_high = `97.5 %`) |>
    dplyr::mutate(sig = p.value < 0.05) |>
    dplyr::right_join(.y_structure, by = "term") |>
    dplyr::mutate(y_display = factor(y_display, levels = .y_levels))
}

.d_attack <- .prep_model_v2(Coefficients_Placebo_Attack_final_model)
.d_lung   <- .prep_model_v2(Coefficients_Placebo_Lung_final_model)
.d_acq    <- .prep_model_v2(Coefficients_Placebo_ACQ_final_model)

# 3.3  Panel builder
.make_fp_v2 <- function(data, title, show_y = FALSE, show_legend = FALSE,
                        sig_colour = "darkred") {
  data_pts <- dplyr::filter(data, !is.na(estimate))
  
  ggplot2::ggplot(data, ggplot2::aes(x = estimate, y = y_display)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbarh(
      data = data_pts,
      ggplot2::aes(xmin = ci_low, xmax = ci_high),
      height = 0.3, linewidth = 0.5, colour = "grey30"
    ) +
    ggplot2::geom_point(
      data = data_pts,
      ggplot2::aes(colour = sig, fill = sig),
      shape = 21, size = 2.0, stroke = 0.6
    ) +
    ggplot2::scale_colour_manual(
      name   = NULL,
      values = c("TRUE" = sig_colour, "FALSE" = "grey40"),
      labels = c("TRUE" = "p < 0.05", "FALSE" = "p \u2265 0.05")
    ) +
    ggplot2::scale_fill_manual(
      name   = NULL,
      values = c("TRUE" = sig_colour, "FALSE" = "white"),
      labels = c("TRUE" = "p < 0.05", "FALSE" = "p \u2265 0.05")
    ) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(x = title, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y        = if (show_y) ggtext::element_markdown(size = 9, colour = "black") else ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      legend.position    = "none",
      panel.grid.major.y = ggplot2::element_line(colour = "grey92"),
      panel.grid.minor   = ggplot2::element_blank()
    )
}

# 3.4  Explicit x-axis limits (shared between main panels and arrow panels)
.xlim_attack <- c(-0.8,  0.8)
.xlim_lung   <- c(-400, 400)
.xlim_acq    <- c(-0.8,  0.8)

# l = 5 on p1 to keep space for y-axis labels; l = 1 on p2/p3 to close inter-column gap
.p1 <- .make_fp_v2(.d_attack, expression(Delta~"Asthma Attack Rate"), show_y = TRUE,
                   sig_colour = "darkred") +
  ggplot2::coord_cartesian(xlim = .xlim_attack) +
  ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 0, 5))
.p2 <- .make_fp_v2(.d_lung,   expression(Delta~"FEV1 (mL)"),          show_y = FALSE,
                   sig_colour = "darkblue") +
  ggplot2::coord_cartesian(xlim = .xlim_lung) +
  ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 0, 1))
.p3 <- .make_fp_v2(.d_acq,    expression(Delta~"ACQ-5"),              show_y = FALSE,
                   sig_colour = "darkgreen") +
  ggplot2::coord_cartesian(xlim = .xlim_acq) +
  ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 0, 1))

# 3.5  Arrow panel builder — x limits matched to the main panel above it
#      Arrows at top of panel (nearest to x-axis title); labels below the arrow.
.make_arrow_panel <- function(xlim_range) {
  .gap <- xlim_range[2] * 0.0625   # proportional center gap (e.g. 0.05 for ±0.8)
  
  ggplot2::ggplot() +
    # Arrows at top of panel, nearest to x-axis title
    ggplot2::annotate("segment",
                      x = -.gap, xend = xlim_range[1] * 0.87, y = 0.85, yend = 0.85,
                      arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                      linewidth = 0.4
    ) +
    ggplot2::annotate("segment",
                      x =  .gap, xend = xlim_range[2] * 0.87, y = 0.85, yend = 0.85,
                      arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                      linewidth = 0.4
    ) +
    # Labels below the arrow (vjust = 1 anchors text bottom just under the arrow)
    ggplot2::annotate("text",
                      x = xlim_range[1] * 0.5, y = 0.78,
                      label = "Reduced placebo response",
                      vjust = 1, hjust = 0.5, size = 2.6
    ) +
    ggplot2::annotate("text",
                      x = xlim_range[2] * 0.5, y = 0.78,
                      label = "Greater placebo response",
                      vjust = 1, hjust = 0.5, size = 2.6
    ) +
    ggplot2::scale_x_continuous(
      limits = xlim_range,
      expand = ggplot2::expansion(mult = 0.02)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 0, b = 5))
}

.arr1 <- .make_arrow_panel(.xlim_attack)
.arr2 <- .make_arrow_panel(.xlim_lung)
.arr3 <- .make_arrow_panel(.xlim_acq)

# 3.6  aRC (95% CI) text panel builder
.make_text_panel <- function(data, fmt = "%.2f") {
  data_text <- data |>
    dplyr::mutate(
      rc_label = dplyr::if_else(
        !is.na(estimate),
        sprintf(paste0(fmt, " (", fmt, ", ", fmt, ")"), estimate, ci_low, ci_high),
        ""
      )
    )
  
  ggplot2::ggplot(data_text, ggplot2::aes(x = 0.05, y = y_display)) +
    # Column header placed at the topmost row for perfect alignment with the panel
    ggplot2::geom_text(
      data  = dplyr::filter(data_text,
                            as.integer(y_display) == max(as.integer(y_display), na.rm = TRUE)),
      ggplot2::aes(label = "aRC (95% CI)"),
      hjust = 0, size = 2.8, fontface = "bold", colour = "black"
    ) +
    # Values for predictor rows
    ggplot2::geom_text(
      data  = dplyr::filter(data_text, !is.na(estimate)),
      ggplot2::aes(label = rc_label),
      hjust = 0, size = 2.6, colour = "black"
    ) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, 1),
                                expand = ggplot2::expansion(mult = 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(5, 1, 0, 2))
}

.txt1 <- .make_text_panel(.d_attack, fmt = "%.2f")
.txt2 <- .make_text_panel(.d_lung,   fmt = "%.0f")
.txt3 <- .make_text_panel(.d_acq,    fmt = "%.2f")

# Flat layout — avoids nested patchwork padding that creates large inter-column gaps.
# Row 1: p | txt | p | txt | p | txt
# Row 2: arr | spacer | arr | spacer | arr | spacer
forest_plot <- patchwork::wrap_plots(
  .p1,  .txt1, .p2,  .txt2, .p3,  .txt3,
  .arr1, patchwork::plot_spacer(),
  .arr2, patchwork::plot_spacer(),
  .arr3, patchwork::plot_spacer(),
  nrow    = 2,
  byrow   = TRUE,
  widths  = c(5, 1.5, 5, 1.5, 5, 1.5),
  heights = c(10, 1.5)
)

forest_plot

ggplot2::ggsave(
  filename = file.path(Figure_path, "Multivariate_forrestplot_3models.pdf"),
  plot     = forest_plot,
  width    = 16,
  height   = 9
)

###############################################################################
###############################################################################
###############################################################################
# 4.  Multi-panel forest plot: predictors of placebo response — vertical layout
#     Panels stacked top (attack) → middle (lung) → bottom (ACQ)
#     Each panel shows its own y-axis labels; text column to the right.

.vp1 <- .make_fp_v2(.d_attack, expression(Delta~"Asthma Attack Rate"), show_y = TRUE,
                    sig_colour = "darkred") +
  ggplot2::coord_cartesian(xlim = .xlim_attack) +
  ggplot2::labs(tag = "A") +
  ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 0, 5),
                 plot.tag    = ggplot2::element_text(face = "bold", size = 14))

.vp2 <- .make_fp_v2(.d_lung,   expression(Delta~"FEV1 (mL)"),          show_y = TRUE,
                    sig_colour = "darkblue") +
  ggplot2::coord_cartesian(xlim = .xlim_lung) +
  ggplot2::labs(tag = "B") +
  ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 0, 5),
                 plot.tag    = ggplot2::element_text(face = "bold", size = 14))

.vp3 <- .make_fp_v2(.d_acq,    expression(Delta~"ACQ-5"),               show_y = TRUE,
                    sig_colour = "darkgreen") +
  ggplot2::coord_cartesian(xlim = .xlim_acq) +
  ggplot2::labs(tag = "C") +
  ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 0, 5),
                 plot.tag    = ggplot2::element_text(face = "bold", size = 14))

# Arrow panels (same builder as section 3)
.varr1 <- .make_arrow_panel(.xlim_attack)
.varr2 <- .make_arrow_panel(.xlim_lung)
.varr3 <- .make_arrow_panel(.xlim_acq)

# Text panels (same builder as section 3)
.vtxt1 <- .make_text_panel(.d_attack, fmt = "%.2f")
.vtxt2 <- .make_text_panel(.d_lung,   fmt = "%.0f")
.vtxt3 <- .make_text_panel(.d_acq,    fmt = "%.2f")

# Vertical layout:
# Col 1 (forest panel)  Col 2 (text)
# Row 1:  vp1            vtxt1
# Row 2:  varr1          spacer
# Row 3:  vp2            vtxt2
# Row 4:  varr2          spacer
# Row 5:  vp3            vtxt3
# Row 6:  varr3          spacer
forest_plot_vertical <- patchwork::wrap_plots(
  .vp1,   .vtxt1,
  .varr1, patchwork::plot_spacer(),
  .vp2,   .vtxt2,
  .varr2, patchwork::plot_spacer(),
  .vp3,   .vtxt3,
  .varr3, patchwork::plot_spacer(),
  ncol    = 2,
  byrow   = TRUE,
  widths  = c(5, 1.5),
  heights = c(10, 1.5, 10, 1.5, 10, 1.5)
)

forest_plot_vertical

ggplot2::ggsave(
  filename = file.path(Figure_path, "Multivariate_forrestplot_3models_vertical.pdf"),
  plot     = forest_plot_vertical,
  #device   = cairo_pdf,
  width    = 9,
  height   = 12
)

ggplot2::ggsave(
  filename = file.path(Figure_path, "Multivariate_forrestplot_3models_vertical.png"),
  plot     = forest_plot_vertical,
  width    = 9,
  height   = 20
)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

# FIGURE 3 : SPLINE CURVES

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

#######################################################################################################################################
#######################################################################################################################################
# Figure 3 - A.1
#######################################################################################################################################

# 4.	Spline curve of BEC for asthma attack risk

summary(ORACLE_Placebo_Attack_imputed$Placebo_response_ACQ)

covariables_list<- c(
  "Sex",
  "Smoking_history_yes_no",
  "GINA_step_numeric",
  "ACQ_score_0W",
  "Attack_history_num",
  "FEV1_preBD_L_0W_per_10",
  "FEV1_reversibility_0W_per_10",
  "FeNO_log10"
)

#--------------------------------------------------------------------------------------------------------------------------------------
# 
# ── Step 1: Fit model ─────────────────────────────────────────────────────────
Model_BEC_rcs5k <- IPD_one_stage(
  data           = ORACLE_Placebo_Attack_imputed,
  outcome_var    = "Placebo_response_Attack",
  predictor_vars = "BEC",
  covariables    = covariables_list,
  
  imp_col         = ".imp",
  model_type      = "lm",
  followup_offset = "Yes",
  followup_col    = "Follow_up_duration_days",
  
  stratified_intercept_var    = "Enrolled_Trial_name",
  random_intercept_var        = "Enrolled_Trial_name",
  predictor_vars_random_slope = "BEC_log10",
  covariables_random_slope    = "Attack_history_num",
  
  spline_terms            = c("BEC"),
  spline_knots_percentile = c(5, 27.5, 50, 72.5, 95),
  
  model_performance  = TRUE,
  weighted_intercept = TRUE,
  parallel           = TRUE,
  n_cores            = 2
)

# ── Step 2: Extract model components ─────────────────────────────────────────
spline_info      <- attr(Model_BEC_rcs5k, "spline_info")
models_list      <- attr(Model_BEC_rcs5k, "models")
knots_bec       <- spline_info[["BEC"]]$knots
basis_names_bec <- spline_info[["BEC"]]$basis_names

dat_ref <- ORACLE_Placebo_Attack_imputed[ORACLE_Placebo_Attack_imputed$.imp == 1L, ]
dat_ref

# ── Step 3: Build prediction frame ────────────────────────────────────────────
# x_seq must be on the BEC scale (10^9/L) — the same scale used to fit the model
# and to compute knots_bec. BEC_log10 was only used for modelling internally.
x_seq <- seq(
  quantile(dat_ref$BEC, 0.01, na.rm = TRUE),
  quantile(dat_ref$BEC, 0.99, na.rm = TRUE),
  length.out = 100L
)

pred_BEC_rcs5k <- data.frame(BEC = x_seq)
pred_BEC_rcs5k

for (v in covariables_list) {
  x <- dat_ref[[v]]
  pred_BEC_rcs5k[[v]] <- if (is.factor(x)) {
    factor(names(which.max(table(x))), levels = levels(x))
  } else {
    median(x, na.rm = TRUE)
  }
}

pred_BEC_rcs5k$Follow_up_duration_days <- 365
# Must use the reference level (level 1) of Enrolled_Trial_name so that
# predict() is anchored at beta_0. link_shift = beta_w - beta_0 then correctly
# shifts to the population-weighted intercept. Using any other trial level
# would double-count its fixed effect on top of the shift.
pred_BEC_rcs5k$Enrolled_Trial_name <- factor(
  levels(dat_ref$Enrolled_Trial_name)[1],
  levels = levels(dat_ref$Enrolled_Trial_name)
)

pred_BEC_rcs5k[, basis_names_bec] <- as.matrix(rms::rcs(x_seq, parms = knots_bec))
pred_BEC_rcs5k

# ── Step 4: Rubin-pooled predictions (response scale, lm / identity link) ─────
# lmerMod does not support se.fit in predict(); SEs are derived from
# the fixed-effects variance-covariance matrix: SE = sqrt(diag(X V X')).
# Trial dummy columns in mm are replaced with population weights so that
# SE reflects Var(β_weighted) — matching the link_shift applied to the mean —
# rather than Var(β_AZISAST_reference_trial), which inflates CIs by ~4×.
n_imp    <- length(models_list)
fits_mat <- matrix(NA_real_, 100L, n_imp)
ses_mat  <- matrix(NA_real_, 100L, n_imp)

trial_counts_bec  <- table(dat_ref$Enrolled_Trial_name)
trial_weights_bec <- trial_counts_bec / sum(trial_counts_bec)

for (j in seq_len(n_imp)) {
  m <- models_list[[j]]
  if (is.null(m)) next
  tryCatch({
    fit_val <- predict(m, newdata = pred_BEC_rcs5k, re.form = NA)
    mm_raw   <- model.matrix(delete.response(terms(reformulas::nobars(formula(m)))), data = pred_BEC_rcs5k)
    vcov_nm  <- rownames(as.matrix(vcov(m)))
    mm       <- matrix(0, nrow(mm_raw), length(vcov_nm), dimnames = list(NULL, vcov_nm))
    mm[, intersect(vcov_nm, colnames(mm_raw))] <- mm_raw[, intersect(vcov_nm, colnames(mm_raw))]
    # Replace trial dummy columns with population weights
    trial_cols <- grep("Enrolled_Trial_name", vcov_nm, value = TRUE)
    for (col in trial_cols) {
      trial_name <- sub(".*\\)", "", col)
      if (trial_name %in% names(trial_weights_bec))
        mm[, col] <- as.numeric(trial_weights_bec[trial_name])
    }
    se_val   <- sqrt(diag(mm %*% as.matrix(vcov(m)) %*% t(mm)))
    fits_mat[, j] <- fit_val
    ses_mat[, j]  <- se_val
  }, error = function(e) NULL)
}

n_ok     <- rowSums(!is.na(fits_mat))
mean_fit <- rowMeans(fits_mat, na.rm = TRUE)
U_bar    <- rowMeans(ses_mat^2, na.rm = TRUE)
B        <- apply(fits_mat, 1L, var, na.rm = TRUE)
T_se     <- sqrt(U_bar + (1 + 1 / n_ok) * B)

# ── Step 4b: Shift to weighted intercept ──────────────────────────────────────
# predict(..., re.form = NA) anchors at beta_0 (reference trial fixed intercept).
# Shift to the study-prevalence-weighted intercept for population-level prediction.
beta_0     <- Model_BEC_rcs5k$Intercept$estimate[Model_BEC_rcs5k$Intercept$term == "(Intercept)"]
beta_w     <- Model_BEC_rcs5k$Weighted_intercept$estimate
link_shift <- beta_w - beta_0
mean_fit_w <- mean_fit + link_shift



# lm uses identity link: predictions are already on the response scale
pred_BEC_rcs5k$prediction <- mean_fit_w
pred_BEC_rcs5k$lower_ci   <- mean_fit_w - 1.96 * T_se
pred_BEC_rcs5k$upper_ci   <- mean_fit_w + 1.96 * T_se
pred_BEC_rcs5k

# ── Step 5: Plot ──────────────────────────────────────────────────────────────
# plotmath expression() labels are used for axis titles so that Delta (∆) and
# 10^9 render correctly in every device (RStudio viewer, pdf, cairo_pdf)
# without Unicode conversion warnings.
plt_BEC_rcs5k <- ggplot(pred_BEC_rcs5k, aes(x = BEC, y = prediction)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              fill = "darkred", alpha = 0.3) +
  geom_line(color = "darkred", linewidth = 1, linetype = "solid") +
  scale_x_log10(breaks = c(0.1, 0.15, 0.3, 0.6, 1.0, 1.5),
                limits = c(0.1, 1.50)) +
  scale_y_continuous(limits = c(0.5, 1.5),
                     breaks = seq(0.5, 1.5, by = 0.1)) +
  labs(
    x = expression(bold("Blood eosinophil count (" %*% 10^9 * "/L)")),
    y = "Placebo Response on\nAsthma Attack annualized rate"
  ) +
  theme_minimal() +
  theme(
    axis.title.x     = element_text(face = "bold", size = 12),
    axis.title.y     = element_text(face = "bold", size = 12),
    axis.text.x      = element_text(face = "bold", size = 10),
    axis.text.y      = element_text(face = "bold", size = 10),
    axis.line.x      = element_line(colour = "black", linewidth = 1),
    axis.line.y      = element_line(colour = "black", linewidth = 1),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.2, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(t = 10, r = 10, b = 10, l = 10)
  )

# ── Step 6: Assemble result list ──────────────────────────────────────────────
Spline_BEC_rcs5k <- list(
  model       = Model_BEC_rcs5k,
  predictions = pred_BEC_rcs5k,
  knots       = knots_bec,
  plot        = plt_BEC_rcs5k
)
Spline_BEC_rcs5k
Spline_BEC_rcs5k$plot
Spline_BEC_rcs5k$knots

Spline_BEC_rcs5k_plot <- Spline_BEC_rcs5k$plot +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
Spline_BEC_rcs5k_plot


ggsave(filename = file.path(Figure_path, "Spline_BEC_rcs5k_Attack.png"),
       plot     = Spline_BEC_rcs5k_plot,
       width    = 5, height = 4)


ggsave(filename = file.path(Figure_path, "Spline_BEC_rcs5k_Attack.pdf"),
       plot     = Spline_BEC_rcs5k_plot,
       width    = 5, height = 4)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# Figure 3 - A.2
#######################################################################################################################################

# 5.	Spline curve of FeNO for asthma attack risk

covariables_list_feno <- c(
  "Sex",
  "Smoking_history_yes_no",
  "GINA_step_numeric",
  "ACQ_score_0W",
  "Attack_history_num",
  "FEV1_preBD_L_0W_per_10",
  "FEV1_reversibility_0W_per_10",
  "BEC_log10"          # BEC_log10 replaces FeNO_log10 (FeNO is now the predictor)
)

# ── Step 1: Fit model ─────────────────────────────────────────────────────────
Model_FeNO_rcs5k <- IPD_one_stage(
  data           = ORACLE_Placebo_Attack_imputed,
  outcome_var    = "Placebo_response_Attack",
  predictor_vars = "FeNO",
  covariables    = covariables_list_feno,
  
  imp_col         = ".imp",
  model_type      = "lm",
  followup_offset = "Yes",
  followup_col    = "Follow_up_duration_days",
  
  stratified_intercept_var    = "Enrolled_Trial_name",
  random_intercept_var        = "Enrolled_Trial_name",
  predictor_vars_random_slope = "FeNO_log10",
  covariables_random_slope    = "Attack_history_num",
  
  spline_terms            = c("FeNO"),
  spline_knots_percentile = c(5, 27.5, 50, 72.5, 95),
  
  model_performance  = TRUE,
  weighted_intercept = TRUE,
  parallel           = TRUE,
  n_cores            = 2
)

# ── Step 2: Extract model components ─────────────────────────────────────────
spline_info_feno      <- attr(Model_FeNO_rcs5k, "spline_info")
models_list_feno      <- attr(Model_FeNO_rcs5k, "models")
knots_feno            <- spline_info_feno[["FeNO"]]$knots
basis_names_feno      <- spline_info_feno[["FeNO"]]$basis_names

dat_ref_feno <- ORACLE_Placebo_Attack_imputed[ORACLE_Placebo_Attack_imputed$.imp == 1L, ]

# ── Step 3: Build prediction frame ────────────────────────────────────────────
# x_seq on the FeNO scale (ppb) — the same scale used to fit the model.
x_seq_feno <- seq(
  quantile(dat_ref_feno$FeNO, 0.01, na.rm = TRUE),
  quantile(dat_ref_feno$FeNO, 0.99, na.rm = TRUE),
  length.out = 100L
)

pred_FeNO_rcs5k <- data.frame(FeNO = x_seq_feno)

for (v in covariables_list_feno) {
  x <- dat_ref_feno[[v]]
  pred_FeNO_rcs5k[[v]] <- if (is.factor(x)) {
    factor(names(which.max(table(x))), levels = levels(x))
  } else {
    median(x, na.rm = TRUE)
  }
}

pred_FeNO_rcs5k$Follow_up_duration_days <- 365
pred_FeNO_rcs5k$Enrolled_Trial_name <- factor(
  levels(dat_ref_feno$Enrolled_Trial_name)[1],
  levels = levels(dat_ref_feno$Enrolled_Trial_name)
)

pred_FeNO_rcs5k[, basis_names_feno] <- as.matrix(rms::rcs(x_seq_feno, parms = knots_feno))

# ── Step 4: Rubin-pooled predictions ─────────────────────────────────────────
n_imp_feno    <- length(models_list_feno)
fits_mat_feno <- matrix(NA_real_, 100L, n_imp_feno)
ses_mat_feno  <- matrix(NA_real_, 100L, n_imp_feno)

trial_counts_feno  <- table(dat_ref_feno$Enrolled_Trial_name)
trial_weights_feno <- trial_counts_feno / sum(trial_counts_feno)

for (j in seq_len(n_imp_feno)) {
  mj <- models_list_feno[[j]]
  if (is.null(mj)) next
  tryCatch({
    fit_val  <- predict(mj, newdata = pred_FeNO_rcs5k, re.form = NA)
    mm_raw   <- model.matrix(delete.response(terms(reformulas::nobars(formula(mj)))),
                             data = pred_FeNO_rcs5k)
    vcov_nm  <- rownames(as.matrix(vcov(mj)))
    mm       <- matrix(0, nrow(mm_raw), length(vcov_nm), dimnames = list(NULL, vcov_nm))
    mm[, intersect(vcov_nm, colnames(mm_raw))] <- mm_raw[, intersect(vcov_nm, colnames(mm_raw))]
    trial_cols <- grep("Enrolled_Trial_name", vcov_nm, value = TRUE)
    for (col in trial_cols) {
      trial_name <- sub(".*\\)", "", col)
      if (trial_name %in% names(trial_weights_feno))
        mm[, col] <- as.numeric(trial_weights_feno[trial_name])
    }
    ses_mat_feno[, j]  <- sqrt(diag(mm %*% as.matrix(vcov(mj)) %*% t(mm)))
    fits_mat_feno[, j] <- fit_val
  }, error = function(e) NULL)
}

n_ok_feno     <- rowSums(!is.na(fits_mat_feno))
mean_fit_feno <- rowMeans(fits_mat_feno, na.rm = TRUE)
U_bar_feno    <- rowMeans(ses_mat_feno^2, na.rm = TRUE)
B_feno        <- apply(fits_mat_feno, 1L, var, na.rm = TRUE)
T_se_feno     <- sqrt(U_bar_feno + (1 + 1 / n_ok_feno) * B_feno)

# ── Step 4b: Shift to weighted intercept ─────────────────────────────────────
beta_0_feno     <- Model_FeNO_rcs5k$Intercept$estimate[Model_FeNO_rcs5k$Intercept$term == "(Intercept)"]
beta_w_feno     <- Model_FeNO_rcs5k$Weighted_intercept$estimate
link_shift_feno <- beta_w_feno - beta_0_feno
mean_fit_w_feno <- mean_fit_feno + link_shift_feno

pred_FeNO_rcs5k$prediction <- mean_fit_w_feno
pred_FeNO_rcs5k$lower_ci   <- mean_fit_w_feno - 1.96 * T_se_feno
pred_FeNO_rcs5k$upper_ci   <- mean_fit_w_feno + 1.96 * T_se_feno
pred_FeNO_rcs5k

# ── Step 5: Plot ──────────────────────────────────────────────────────────────
plt_FeNO_rcs5k <- ggplot(pred_FeNO_rcs5k, aes(x = FeNO, y = prediction)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              fill = "darkred", alpha = 0.3) +
  geom_line(color = "darkred", linewidth = 1, linetype = "solid") +
  scale_x_log10(breaks = c(10, 25, 25, 50,75, 100, 200),
                limits = c(
                  10,
                  200
                )) +
  scale_y_continuous(limits = c(0.5, 1.5),
                     breaks = seq(0.5, 1.5, by = 0.1)) +
  labs(
    x = expression(bold("Fractional exhaled nitric oxide (ppb)")),
    y = "Placebo Response on\nAsthma Attack annualized rate"
  ) +
  theme_minimal() +
  theme(
    axis.title.x     = element_text(face = "bold", size = 12),
    axis.title.y     = element_text(face = "bold", size = 12),
    axis.text.x      = element_text(face = "bold", size = 10),
    axis.text.y      = element_text(face = "bold", size = 10),
    axis.line.x      = element_line(colour = "black", linewidth = 1),
    axis.line.y      = element_line(colour = "black", linewidth = 1),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.2, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(t = 10, r = 10, b = 10, l = 10)
  )

# ── Step 6: Assemble result list ──────────────────────────────────────────────
Spline_FeNO_rcs5k <- list(
  model       = Model_FeNO_rcs5k,
  predictions = pred_FeNO_rcs5k,
  knots       = knots_feno,
  plot        = plt_FeNO_rcs5k
)
Spline_FeNO_rcs5k
Spline_FeNO_rcs5k$plot
Spline_FeNO_rcs5k$knots

Spline_FeNO_rcs5k_plot <- Spline_FeNO_rcs5k$plot +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
Spline_FeNO_rcs5k_plot

ggsave(filename = file.path(Figure_path, "Spline_FeNO_rcs5k_Attack.png"),
       plot     = Spline_FeNO_rcs5k_plot,
       width    = 5, height = 4)


ggsave(filename = file.path(Figure_path, "Spline_FeNO_rcs5k_Attack.pdf"),
       plot     = Spline_FeNO_rcs5k_plot,
       width    = 5, height = 4)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# Figure 3 - B
#######################################################################################################################################

# 6.	Spline curve of FeNO for lung function

covariables_list_lung_feno <- c(
  "Age_per_10",
  "BMI_per_5",
  "GINA_step_numeric",
  "FEV1_preBD_L_0W_per_10",
  "FEV1_reversibility_0W_per_10"
)

# ── Step 1: Fit model ─────────────────────────────────────────────────────────
Model_Lung_FeNO_rcs5k <- IPD_one_stage(
  data           = ORACLE_Placebo_Lung_imputed,
  outcome_var    = "Placebo_response_Lung",
  predictor_vars = "FeNO",
  covariables    = covariables_list_lung_feno,
  
  imp_col         = ".imp",
  model_type      = "lm",
  followup_offset = "No",
  
  stratified_intercept_var    = "Enrolled_Trial_name",
  random_intercept_var        = "Enrolled_Trial_name",
  predictor_vars_random_slope = "FeNO_log10",
  covariables_random_slope    = "FEV1_preBD_L_0W",
  
  spline_terms            = c("FeNO"),
  spline_knots_percentile = c(5, 27.5, 50, 72.5, 95),
  
  model_performance  = TRUE,
  weighted_intercept = TRUE,
  parallel           = TRUE,
  n_cores            = 2
)

# ── Step 2: Extract model components ─────────────────────────────────────────
spline_info_lung_feno      <- attr(Model_Lung_FeNO_rcs5k, "spline_info")
models_list_lung_feno      <- attr(Model_Lung_FeNO_rcs5k, "models")
knots_lung_feno            <- spline_info_lung_feno[["FeNO"]]$knots
basis_names_lung_feno      <- spline_info_lung_feno[["FeNO"]]$basis_names

dat_ref_lung_feno <- ORACLE_Placebo_Lung_imputed[ORACLE_Placebo_Lung_imputed$.imp == 1L, ]
dat_ref_lung_feno
# ── Step 3: Build prediction frame ────────────────────────────────────────────
x_seq_lung_feno <- seq(
  quantile(dat_ref_lung_feno$FeNO, 0.01, na.rm = TRUE),
  quantile(dat_ref_lung_feno$FeNO, 0.99, na.rm = TRUE),
  length.out = 100L
)

pred_Lung_FeNO_rcs5k <- data.frame(FeNO = x_seq_lung_feno)


for (v in covariables_list_lung_feno) {
  x <- dat_ref_lung_feno[[v]]
  pred_Lung_FeNO_rcs5k[[v]] <- if (is.factor(x)) {
    factor(names(which.max(table(x))), levels = levels(x))
  } else {
    median(x, na.rm = TRUE)
  }
}

# Use droplevels() to get only trials actually present in the lung dataset.
# The inherited factor has 22 levels but only 8 are populated; picking level [1]
# without dropping gives AZISAST (0 patients), which the lung model never saw.
lung_feno_present_levels <- levels(droplevels(dat_ref_lung_feno$Enrolled_Trial_name))
pred_Lung_FeNO_rcs5k$Enrolled_Trial_name <- factor(
  lung_feno_present_levels[1],
  levels = lung_feno_present_levels
)

pred_Lung_FeNO_rcs5k[, basis_names_lung_feno] <- as.matrix(
  rms::rcs(x_seq_lung_feno, parms = knots_lung_feno)
)
pred_Lung_FeNO_rcs5k
# ── Step 4: Rubin-pooled predictions ─────────────────────────────────────────
n_imp_lung_feno    <- length(models_list_lung_feno)
fits_mat_lung_feno <- matrix(NA_real_, 100L, n_imp_lung_feno)
ses_mat_lung_feno  <- matrix(NA_real_, 100L, n_imp_lung_feno)


trial_counts_lung_feno  <- table(droplevels(dat_ref_lung_feno$Enrolled_Trial_name))
trial_weights_lung_feno <- trial_counts_lung_feno / sum(trial_counts_lung_feno)


for (j in seq_len(n_imp_lung_feno)) {
  mj <- models_list_lung_feno[[j]]
  if (is.null(mj)) next
  tryCatch({
    fit_val  <- predict(mj, newdata = pred_Lung_FeNO_rcs5k, re.form = NA)
    mm_raw   <- model.matrix(delete.response(terms(reformulas::nobars(formula(mj)))),
                             data = pred_Lung_FeNO_rcs5k)
    vcov_nm  <- rownames(as.matrix(vcov(mj)))
    mm       <- matrix(0, nrow(mm_raw), length(vcov_nm), dimnames = list(NULL, vcov_nm))
    mm[, intersect(vcov_nm, colnames(mm_raw))] <- mm_raw[, intersect(vcov_nm, colnames(mm_raw))]
    trial_cols <- grep("Enrolled_Trial_name", vcov_nm, value = TRUE)
    for (col in trial_cols) {
      trial_name <- sub(".*\\)", "", col)
      if (trial_name %in% names(trial_weights_lung_feno))
        mm[, col] <- as.numeric(trial_weights_lung_feno[trial_name])
    }
    ses_mat_lung_feno[, j]  <- sqrt(diag(mm %*% as.matrix(vcov(mj)) %*% t(mm)))
    fits_mat_lung_feno[, j] <- fit_val
  }, error = function(e) NULL)
}

n_ok_lung_feno     <- rowSums(!is.na(fits_mat_lung_feno))
mean_fit_lung_feno <- rowMeans(fits_mat_lung_feno, na.rm = TRUE)
U_bar_lung_feno    <- rowMeans(ses_mat_lung_feno^2, na.rm = TRUE)
B_lung_feno        <- apply(fits_mat_lung_feno, 1L, var, na.rm = TRUE)
T_se_lung_feno     <- sqrt(U_bar_lung_feno + (1 + 1 / n_ok_lung_feno) * B_lung_feno)

# ── Step 4b: Shift to weighted intercept ─────────────────────────────────────
beta_0_lung_feno     <- Model_Lung_FeNO_rcs5k$Intercept$estimate[
  Model_Lung_FeNO_rcs5k$Intercept$term == "(Intercept)"
]
beta_w_lung_feno     <- Model_Lung_FeNO_rcs5k$Weighted_intercept$estimate
link_shift_lung_feno <- beta_w_lung_feno - beta_0_lung_feno
mean_fit_w_lung_feno <- mean_fit_lung_feno + link_shift_lung_feno

pred_Lung_FeNO_rcs5k$prediction <- mean_fit_w_lung_feno
pred_Lung_FeNO_rcs5k$lower_ci   <- mean_fit_w_lung_feno - 1.96 * T_se_lung_feno
pred_Lung_FeNO_rcs5k$upper_ci   <- mean_fit_w_lung_feno + 1.96 * T_se_lung_feno
pred_Lung_FeNO_rcs5k

# ── Step 5: Plot ──────────────────────────────────────────────────────────────
plt_Lung_FeNO_rcs5k <- ggplot(pred_Lung_FeNO_rcs5k, aes(x = FeNO, y = prediction)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              fill = "darkblue", alpha = 0.3) +
  geom_line(color = "darkblue", linewidth = 1, linetype = "solid") +
  scale_x_log10(breaks = c(10, 15, 25, 50, 75, 100, 200),
                limits = c(
                  9.88078,
                  200
                )) +
  scale_y_continuous(limits = c(0, 250),
                     breaks = seq(0, 250, by = 25)) +
  labs(
    x = expression(bold("Fractional exhaled nitric oxide (ppb)")),
    y = "Placebo Response on\nFEV1 (mL)"
  ) +
  theme_minimal() +
  theme(
    axis.title.x     = element_text(face = "bold", size = 12),
    axis.title.y     = element_text(face = "bold", size = 12),
    axis.text.x      = element_text(face = "bold", size = 10),
    axis.text.y      = element_text(face = "bold", size = 10),
    axis.line.x      = element_line(colour = "black", linewidth = 1),
    axis.line.y      = element_line(colour = "black", linewidth = 1),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.2, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(t = 10, r = 10, b = 10, l = 10)
  )

# ── Step 6: Assemble result list ──────────────────────────────────────────────
Spline_Lung_FeNO_rcs5k <- list(
  model       = Model_Lung_FeNO_rcs5k,
  predictions = pred_Lung_FeNO_rcs5k,
  knots       = knots_lung_feno,
  plot        = plt_Lung_FeNO_rcs5k
)
Spline_Lung_FeNO_rcs5k
Spline_Lung_FeNO_rcs5k$plot
Spline_Lung_FeNO_rcs5k$knots

Spline_Lung_FeNO_rcs5k_plot <- Spline_Lung_FeNO_rcs5k$plot +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
Spline_Lung_FeNO_rcs5k_plot

ggsave(filename = file.path(Figure_path, "Spline_Lung_FeNO_rcs5k.png"),
       plot     = Spline_Lung_FeNO_rcs5k_plot,
       width    = 5, height = 4)

ggsave(filename = file.path(Figure_path, "Spline_Lung_FeNO_rcs5k.pdf"),
       plot     = Spline_Lung_FeNO_rcs5k_plot,
       width    = 5, height = 4)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# Figure 3 - C
#######################################################################################################################################

# 7.	Spline curve of FeNO for ACQ

covariables_list_ACQ_feno <- c(
  "Age_per_10",
  "GINA_step_numeric",
  "ACQ_score_0W",
  "Attack_history_num",
  "FEV1_preBD_L_0W_per_10"
)

# ── Step 1: Fit model ─────────────────────────────────────────────────────────
Model_ACQ_FeNO_rcs5k <- IPD_one_stage(
  data           = ORACLE_Placebo_ACQ_imputed,
  outcome_var    = "Placebo_response_ACQ",
  predictor_vars = "FeNO",
  covariables    = covariables_list_ACQ_feno,
  
  imp_col         = ".imp",
  model_type      = "lm",
  followup_offset = "No",
  
  stratified_intercept_var    = "Enrolled_Trial_name",
  random_intercept_var        = "Enrolled_Trial_name",
  predictor_vars_random_slope = "FeNO_log10",
  covariables_random_slope    = "ACQ_score_0W",
  
  spline_terms            = c("FeNO"),
  spline_knots_percentile = c(5, 27.5, 50, 72.5, 95),
  
  model_performance  = TRUE,
  weighted_intercept = TRUE,
  parallel           = TRUE,
  n_cores            = 2
)

# ── Step 2: Extract model components ─────────────────────────────────────────
spline_info_ACQ_feno      <- attr(Model_ACQ_FeNO_rcs5k, "spline_info")
models_list_ACQ_feno      <- attr(Model_ACQ_FeNO_rcs5k, "models")
knots_ACQ_feno            <- spline_info_ACQ_feno[["FeNO"]]$knots
basis_names_ACQ_feno      <- spline_info_ACQ_feno[["FeNO"]]$basis_names

dat_ref_ACQ_feno <- ORACLE_Placebo_ACQ_imputed[ORACLE_Placebo_ACQ_imputed$.imp == 1L, ]

# ── Step 3: Build prediction frame ────────────────────────────────────────────
x_seq_ACQ_feno <- seq(
  quantile(dat_ref_ACQ_feno$FeNO, 0.01, na.rm = TRUE),
  quantile(dat_ref_ACQ_feno$FeNO, 0.99, na.rm = TRUE),
  length.out = 100L
)

pred_ACQ_FeNO_rcs5k <- data.frame(FeNO = x_seq_ACQ_feno)

for (v in covariables_list_ACQ_feno) {
  x <- dat_ref_ACQ_feno[[v]]
  pred_ACQ_FeNO_rcs5k[[v]] <- if (is.factor(x)) {
    factor(names(which.max(table(x))), levels = levels(x))
  } else {
    median(x, na.rm = TRUE)
  }
}

# droplevels() to use only the 8 trials present in the ACQ dataset
ACQ_feno_present_levels <- levels(droplevels(dat_ref_ACQ_feno$Enrolled_Trial_name))
pred_ACQ_FeNO_rcs5k$Enrolled_Trial_name <- factor(
  ACQ_feno_present_levels[1],
  levels = ACQ_feno_present_levels
)

pred_ACQ_FeNO_rcs5k[, basis_names_ACQ_feno] <- as.matrix(
  rms::rcs(x_seq_ACQ_feno, parms = knots_ACQ_feno)
)

# ── Step 4: Rubin-pooled predictions ─────────────────────────────────────────
n_imp_ACQ_feno    <- length(models_list_ACQ_feno)
fits_mat_ACQ_feno <- matrix(NA_real_, 100L, n_imp_ACQ_feno)
ses_mat_ACQ_feno  <- matrix(NA_real_, 100L, n_imp_ACQ_feno)

trial_counts_ACQ_feno  <- table(droplevels(dat_ref_ACQ_feno$Enrolled_Trial_name))
trial_weights_ACQ_feno <- trial_counts_ACQ_feno / sum(trial_counts_ACQ_feno)

for (j in seq_len(n_imp_ACQ_feno)) {
  mj <- models_list_ACQ_feno[[j]]
  if (is.null(mj)) next
  tryCatch({
    fit_val  <- predict(mj, newdata = pred_ACQ_FeNO_rcs5k, re.form = NA)
    mm_raw   <- model.matrix(delete.response(terms(reformulas::nobars(formula(mj)))),
                             data = pred_ACQ_FeNO_rcs5k)
    vcov_nm  <- rownames(as.matrix(vcov(mj)))
    mm       <- matrix(0, nrow(mm_raw), length(vcov_nm), dimnames = list(NULL, vcov_nm))
    mm[, intersect(vcov_nm, colnames(mm_raw))] <- mm_raw[, intersect(vcov_nm, colnames(mm_raw))]
    trial_cols <- grep("Enrolled_Trial_name", vcov_nm, value = TRUE)
    for (col in trial_cols) {
      trial_name <- sub(".*\\)", "", col)
      if (trial_name %in% names(trial_weights_ACQ_feno))
        mm[, col] <- as.numeric(trial_weights_ACQ_feno[trial_name])
    }
    ses_mat_ACQ_feno[, j]  <- sqrt(diag(mm %*% as.matrix(vcov(mj)) %*% t(mm)))
    fits_mat_ACQ_feno[, j] <- fit_val
  }, error = function(e) NULL)
}

n_ok_ACQ_feno     <- rowSums(!is.na(fits_mat_ACQ_feno))
mean_fit_ACQ_feno <- rowMeans(fits_mat_ACQ_feno, na.rm = TRUE)
U_bar_ACQ_feno    <- rowMeans(ses_mat_ACQ_feno^2, na.rm = TRUE)
B_ACQ_feno        <- apply(fits_mat_ACQ_feno, 1L, var, na.rm = TRUE)
T_se_ACQ_feno     <- sqrt(U_bar_ACQ_feno + (1 + 1 / n_ok_ACQ_feno) * B_ACQ_feno)

# ── Step 4b: Shift to weighted intercept ─────────────────────────────────────
beta_0_ACQ_feno     <- Model_ACQ_FeNO_rcs5k$Intercept$estimate[
  Model_ACQ_FeNO_rcs5k$Intercept$term == "(Intercept)"
]
beta_w_ACQ_feno     <- Model_ACQ_FeNO_rcs5k$Weighted_intercept$estimate
link_shift_ACQ_feno <- beta_w_ACQ_feno - beta_0_ACQ_feno
mean_fit_w_ACQ_feno <- mean_fit_ACQ_feno + link_shift_ACQ_feno

pred_ACQ_FeNO_rcs5k$prediction <- mean_fit_w_ACQ_feno
pred_ACQ_FeNO_rcs5k$lower_ci   <- mean_fit_w_ACQ_feno - 1.96 * T_se_ACQ_feno
pred_ACQ_FeNO_rcs5k$upper_ci   <- mean_fit_w_ACQ_feno + 1.96 * T_se_ACQ_feno
pred_ACQ_FeNO_rcs5k
# ── Step 5: Plot ──────────────────────────────────────────────────────────────
plt_ACQ_FeNO_rcs5k <- ggplot(pred_ACQ_FeNO_rcs5k, aes(x = FeNO, y = prediction)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              fill = "darkgreen", alpha = 0.3) +
  geom_line(color = "darkgreen", linewidth = 1, linetype = "solid") +
  scale_x_log10(breaks = c(10,15, 25, 50,75, 100, 200),
                limits = c(
                  9.55045,
                  200
                )) +
  scale_y_continuous(limits = c(0.5, 1.5),
                     breaks = seq(0.5, 1.5, by = 0.1)) +
  labs(
    x = expression(bold("Fractional exhaled nitric oxide (ppb)")),
    y = "Placebo Response on\nACQ-5 Score"
  ) +
  theme_minimal() +
  theme(
    axis.title.x     = element_text(face = "bold", size = 12),
    axis.title.y     = element_text(face = "bold", size = 12),
    axis.text.x      = element_text(face = "bold", size = 10),
    axis.text.y      = element_text(face = "bold", size = 10),
    axis.line.x      = element_line(colour = "black", linewidth = 1),
    axis.line.y      = element_line(colour = "black", linewidth = 1),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "gray90", linewidth = 0.2, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(t = 10, r = 10, b = 10, l = 10)
  )

# ── Step 6: Assemble result list ──────────────────────────────────────────────
Spline_ACQ_FeNO_rcs5k <- list(
  model       = Model_ACQ_FeNO_rcs5k,
  predictions = pred_ACQ_FeNO_rcs5k,
  knots       = knots_ACQ_feno,
  plot        = plt_ACQ_FeNO_rcs5k
)
Spline_ACQ_FeNO_rcs5k
Spline_ACQ_FeNO_rcs5k$plot
Spline_ACQ_FeNO_rcs5k$knots

Spline_ACQ_FeNO_rcs5k_plot <- Spline_ACQ_FeNO_rcs5k$plot +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
Spline_ACQ_FeNO_rcs5k_plot

ggsave(filename = file.path(Figure_path, "Spline_ACQ_FeNO_rcs5k.png"),
       plot     = Spline_ACQ_FeNO_rcs5k_plot,
       width    = 5, height = 4)

ggsave(filename = file.path(Figure_path, "Spline_ACQ_FeNO_rcs5k.pdf"),
       plot     = Spline_ACQ_FeNO_rcs5k_plot,
       width    = 5, height = 4)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# 8.  Multipanel plot: all four spline curves

library(patchwork)

Multipanel_Splines <- (Spline_BEC_rcs5k_plot  | Spline_FeNO_rcs5k_plot) /
  (Spline_Lung_FeNO_rcs5k_plot | Spline_ACQ_FeNO_rcs5k_plot) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 14))
  )

Multipanel_Splines

ggsave(filename = file.path(Figure_path, "Multipanel_Splines.png"),
       plot     = Multipanel_Splines,
       width    = 10, height = 8)

ggsave(filename = file.path(Figure_path, "Multipanel_Splines.pdf"),
       plot     = Multipanel_Splines,
       width    = 10, height = 8)



