#######################################################################################################################################

# SCRIPT FOR THE SUPPLEMENTARY MATERIAL (TABLE S_ + FIGURE S_) FOR THE PLACEBO-ASTHMA ANALYSIS
#JSP and SML et al.

## PlaceboAsthma_MAIN_Script_clean = should be run before.

#######################################################################################################################################
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#load the packages
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

#Load the homemade package MIAnalysis and OSI Calculator
install.packages("devtools")
library(devtools)
required_packages <- c("dplyr", "MASS", "Hmisc", "mice", "rlang", "survival")
devtools::install_github("SamLar27/MIAnalysis", force = TRUE)

install.packages("doParallel")
devtools::install_github("SamLar27/MIAnalysis", force = TRUE)

library(MIAnalysis)

#======================================================================================================================================
#======================================================================================================================================
## PREPARATION OF THE DATASET (MULTIPLE IMPUTATION)
#======================================================================================================================================
#======================================================================================================================================
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Import the DATA from the file
setwd("/Users/joel/ORACLE - Placebo - Exacerbations/Raw_data")

#======================================================================================================================================
## NEW IMPUTATED DATA (.imp1-10)
#FCS_single_pmm_logreg_pmm (for variables : continuous_dichotomic_count/integer)
# New imputated dataset (.imp=0, original, so .imp=1-10, imputated (no NA)) 
#======================================================================================================================================
#Import the DATA imputed from the file
##Selecting the working directory
user_name <- Sys.info()[["user"]]

Load_path <- file.path(
  "/Users/joel/ORACLE - Placebo - Exacerbations/ORACLE_Placebo/Data/Raw_data/Main_FCS_single_pmm_logreg_polyreg_pmm"
)

Load_file <- "ORACLE_MI_FULL_long_imp0_to_m.rds"

full_path <- file.path(Load_path, Load_file)

file.exists(full_path)

data_long <- readRDS(full_path)

# Convert to classic data.frame (not tibble)
ORACLE_after_imputation_sup <- as.data.frame(data_long)
nrow(ORACLE_after_imputation)
ncol(ORACLE_after_imputation)

data_imputated <- ORACLE_after_imputation_sup

#======================================================================================================================================
## ORIGINAL DATA - FOR MISSING DATA (PLOT) - SUPPLEMENTARY
#Import the DATA imputed from the file
##Selecting the working directory

user_name <- Sys.info()[["user"]]

Load_path_2 <- file.path(
  "/Users/joel/ORACLE - Placebo - Exacerbations/ORACLE_Placebo/Data/3_After_part_A_Data_preparation"
)

original_file <- "ORACLE_PartA_prepared_2026_05_01.RData"

full_path_2 <- file.path(Load_path_2, original_file)

file.exists(full_path_2)

# Correct way to import .RData
obj_names <- load(full_path_2)

# Check what was loaded
obj_names
ls()

# Extract the object (assuming only one main dataset)
data_original <- get(obj_names[1])

# Convert to data.frame if needed
ORACLE_original <- as.data.frame(data_original)

nrow(ORACLE_original)
ncol(ORACLE_original)

#======================================================================================================================================
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Selecting the current working project for the project
##Selecting the working directory
setwd("/Users/joel/ORACLE - Placebo - Exacerbations/Project")

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#======================================================================================================================================
#TRANSFORMATION of Values
#Putting NA for the missing value
data_imputated[data_imputated== "NA"] <- NA

###Creation of a unique dataframe

data_imputated_all <- data_imputated
identical(colnames(data_imputated_all),
          colnames(data_original)[colnames(data_imputated_all) %in% colnames(data_original)])

unique(data_imputated_all$.imp)
nrow(data_imputated_all)
ncol(data_imputated_all)
colnames(data_imputated_all)

All_data <-data_imputated_all
nrow(All_data)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#EXCLUSION CRITERIA (NO OPEN LABEL)
# CAPTAIN AND PACT (NOT TRUE PLACEBO ARM) KEPT UNTIL ANALYSIS FOR SENSITIVITY (IN SUPPLEMENTARY)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unique(All_data$Enrolled_Trial_name)

# Removing the trials open label
Excluded_trials <- c("Novel_START", "PRACTICAL")

All_data__without_open_label<- All_data %>%
  filter(!Enrolled_Trial_name %in% Excluded_trials )
nrow(All_data__without_open_label)
unique (All_data__without_open_label$Enrolled_Trial_name)

#MAIN DATASET FOR FURTHER ANALYSIS WILL BE CALLED :
Data_Oracle<-All_data__without_open_label

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Identify the categorical data 

colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Sex" )]<- "Gender"
Data_Oracle$Gender<-as.factor(Data_Oracle$Gender)

colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Smoking_history_yes_no" )]<- "Smoking_Statut"
Data_Oracle$Smoking_Statut<-as.factor(Data_Oracle$Smoking_Statut)

colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Airborne_allergen_sensitisation" )]<- "Airborne_allergen_sensibilisation"
Data_Oracle$Airborne_allergen_sensibilisation<-as.factor(Data_Oracle$Airborne_allergen_sensibilisation)

colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Eczema" )]<- "Eczema"
Data_Oracle$Eczema<-as.factor(Data_Oracle$Eczema)

colnames(Data_Oracle)[which(colnames(Data_Oracle)=="Allergic_Rhinitis" )]<- "Allergic_rhinitis"
Data_Oracle$Allergic_rhinitis<-as.factor(Data_Oracle$Allergic_rhinitis)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Add the imputated original biomarkers value (not log10)
Data_Oracle <- Data_Oracle %>%
  mutate(
    FeNO =
      if_else(
        is.na(FeNO) & 
          !is.na(FeNO_log10),
        10^FeNO_log10,
        FeNO
      ),
    
    BEC = 
      if_else(
        is.na(BEC) & !is.na(BEC_log10),
        (10^BEC_log10) / 1e9,
        BEC
      ),
    
    IgE =
      if_else(
        is.na(IgE) & !is.na(IgE_log10),
        10^IgE_log10,
        IgE
      )
  )

#======================================================================================================================================
## Calculation of the parameters by a definite change
Data_Oracle<-Data_Oracle%>%
  #Age per 10 year increase
  mutate(Age_per10=Age/10) %>%
  
  # mutate(Age_per_10_imputated=Age_imputated/10) %>%
  #BMI per 5 increase
  mutate(BMI_per5=BMI/5) %>%
  
  #FEV1 per 10% decrease
  mutate(FEV1_preBD_PCT_0W_per10=-FEV1_preBD_PCT_0W/10) %>%
  
  # FEV1 pre-BD L per 0.1 decrease in absolute units
  dplyr::mutate(FEV1_preBD_L_0W_per_10 = -FEV1_preBD_L_0W / 0.1) %>%
  
  #Reversibility per 10%
  mutate(FEV1_reversibility_0W_per10= FEV1_reversibility_0W/10)

#======================================================================================================================================

##Category by GINA_treatment_step
Data_Oracle<-Data_Oracle %>%
  mutate(Treatment_step= case_when(Treatment_step=="1"~"Step 1",Treatment_step=="2"~"Step 2",Treatment_step=="3"~"Step 3",Treatment_step=="4"~"Step 4",Treatment_step=="5"~"Step 5",TRUE~Treatment_step))
Data_Oracle$Treatment_step
Data_Oracle$Treatment_step <- factor(Data_Oracle$Treatment_step, levels = c("Step 1","Step 2","Step 3","Step 4","Step 5"))

Data_Oracle<-Data_Oracle %>%
  mutate(Treatment_step_combine= case_when(Treatment_step=="Step 1"~"Step 1-2",Treatment_step=="Step 2"~"Step 1-2",TRUE~Treatment_step))
Data_Oracle$Treatment_step_combine <- factor(Data_Oracle$Treatment_step_combine, levels = c("Step 1-2","Step 3","Step 4","Step 5"))

Data_Oracle$Treatment_step_ordinal <- as.numeric(gsub("Step ", "", as.character(Data_Oracle$Treatment_step)))
summary(Data_Oracle$Treatment_step_ordinal)

##Create a treatment step_1-2_vs_3-4-5
Data_Oracle <- Data_Oracle %>%
  mutate(Treatment_step_1_2vs3_5 = factor(case_when(
    Treatment_step %in% c("Step 1", "Step 2") ~ "Step 1-2",
    Treatment_step %in% c("Step 3", "Step 4", "Step 5") ~ "Step 3-5",
    TRUE ~ Treatment_step  # Keep other values as they are
  )))

##Create a treatment step_1-2_vs_3-4_vs_5
Data_Oracle <- Data_Oracle %>%
  mutate(Treatment_step_1_2vs3_4vs5 = factor(case_when(
    Treatment_step %in% c("Step 1", "Step 2") ~ "Step 1-2",
    Treatment_step %in% c("Step 3", "Step 4") ~ "Step 3-4",
    Treatment_step %in% c("Step 5") ~ "Step 5",
    TRUE ~ Treatment_step  # Keep other values as they are
  )))

Data_Oracle <- Data_Oracle %>%
  mutate(Treatment_step_1_2vs3_4vs5 = relevel(Treatment_step_1_2vs3_4vs5, ref = "Step 3-4"))

##Category by ACQ-5
Data_Oracle$ACQ_0W_by_group<-cut(Data_Oracle$ACQ_score_0W,breaks = c(-10,1.5,3,100000),labels=c("<1.5","1.5-3",">3"))

##Category by BMI
Data_Oracle$BMI_by_group<-cut(Data_Oracle$BMI,breaks = c(-10,25,30,35,100000),labels=c("<25","25-30","30-35",">35"))

##Category by Age group
Data_Oracle$Age_by_group<-cut(Data_Oracle$Age,breaks = c(-10,40,50,60,100000),labels=c("<40","40-50","50-60",">60"))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Data preparation for country
Data_Oracle$Country <- factor(Data_Oracle$Country)
levels(Data_Oracle$Country)

Data_Oracle <- Data_Oracle %>%
  mutate(Country = case_when(
    Country %in% c("United_States", "United States") ~ "United States",
    Country %in% c("Hungary", "Hungria") ~ "Hungary",
    Country %in% c("Russia", "Russian Federation") ~ "Russia",
    Country %in% c("South_Korea", "Korea","Korea, Republic of") ~ "South Korea",
    Country %in% c("South_Africa", "South Africa") ~ "South Africa",
    Country %in% c("UK - CMD", "United Kingdom") ~ "United Kingdom",
    TRUE ~ Country
  ))
Data_Oracle$Country <- factor(Data_Oracle$Country)
#Data_Oracle$Country <- relevel(Data_Oracle$Country, ref = "United States")

Data_Oracle <- Data_Oracle %>%
  mutate(Country_per_region_1 = case_when(
    Country %in% c("United States", "Canada", "Mexico", "North America") ~ "North America",
    Country %in% c("Argentina", "Chile", "Peru", "Brazil", "South America") ~ "South America",
    Country %in% c("Australia", "New Zealand", "Oceania") ~ "Oceania",
    Country == "South Africa" ~ "Africa South",
    Country %in% c("Italy", "United Kingdom", "Belgium", "Poland", "Romania", "Hungary",
                   "Czech Republic", "Belarus", "Ukraine", "Spain", "Slovakia", "Russia",
                   "France", "Serbia", "Bulgaria", "Germany", "Netherlands", "Latvia",
                   "Lithuania", "Europe") ~ "Europe",
    Country %in% c("Japan", "Korea", "Israel", "Turkey", "Vietnam", "South Korea", "Asia") ~ "Asia",
    TRUE ~ NA
  ))

Data_Oracle <- Data_Oracle %>%
  mutate(Country_per_region_2 = case_when(
    Country %in% c("United States", "Canada", "Mexico") ~ "North America",
    Country %in% c("Argentina", "Chile", "Peru", "Brazil") ~ "South America",
    Country %in% c("Australia", "New_Zealand") ~ "Oceania",
    Country == "South_Africa" ~ "Africa_South",
    Country %in% c("Israel", "Turkey") ~ "Middle East",
    Country %in% c("Japan", "South Korea", "Vietnam") ~ "Asia",
    Country %in% c("United Kingdom", "France", "Germany", "Netherlands", "Belgium") ~ "Europe_Western",
    Country %in% c("Spain", "Italy") ~ "Europe_Southern",
    Country %in% c("Poland", "Hungary", "Hungria", "Czech_Republic", "Slovakia") ~ "Europe_Central Eastern",
    Country %in% c("Belarus", "Ukraine", "Russia", "Lithuania", "Latvia") ~ "Europe_Eastern",
    Country %in% c("Romania", "Bulgaria", "Serbia") ~ "Europe_South Eastern",
    TRUE ~ NA
  ))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#PREPARATION OF DATA FOR THE DIFFERENT ANALYSIS
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Placebo response on asthma attack
#Selection of the studies with attack history in numeric
Studies_included_attack<- c("AZISAST","BENRAP2B","DREAM","DRI12544","EXTRA", "LAVOLTA_1","LAVOLTA_2","LUSTER_1","LUSTER_2","NAVIGATOR","PATHWAY","QUEST","STRATOS_1","STRATOS_2")

Data_Oracle_attack<-Data_Oracle %>%
  filter(Enrolled_Trial_name%in%Studies_included_attack)
unique(Data_Oracle_attack$Enrolled_Trial_name)

nrow(Data_Oracle_attack %>%
       filter(.imp==0))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Preparation for delta asthma attack analysis
sum(is.na(Data_Oracle_attack$Attack_number_during_followup))
sum(is.na(Data_Oracle_attack$Attack_history_num))

#Delta attack
Data_Oracle_attack<- Data_Oracle_attack %>%
  mutate(Delta_attack= Attack_history_num-Attack_number_during_followup)

##Calculate Rate Attack in trial
Data_Oracle_attack <- Data_Oracle_attack %>%
  mutate(
    Attack_rate_before_trial = Attack_history_num / 1,
    Attack_rate_during_trial = Attack_number_during_followup / Follow_up_duration_days,
    Delta_attack_rate = Attack_rate_before_trial - Attack_rate_during_trial,
    Attack_rate_before_trial_log = log(Attack_rate_before_trial + 1)
  )

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Dataset for the main analysis
Data_Oracle_attack_imputated <- Data_Oracle_attack %>%
  # Identify participants where Delta_attack is missing in .imp == 0
  group_by(Sequential_number) %>%
  filter(!any(.imp == 0 & is.na(Delta_attack))) %>%
  ungroup() %>%
  # Keep only .imp = 1 to 10
  filter(.imp != 0)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## For lung function change
#Selection of the studies with the total ∆FEV1
Studies_included_lung<- c("BENRAP2B","LAVOLTA_1","LAVOLTA_2","NAVIGATOR","PATHWAY","QUEST","STRATOS_1","STRATOS_2")

Data_Oracle_lung<-Data_Oracle %>%
  filter(Enrolled_Trial_name%in%Studies_included_lung)

nrow(Data_Oracle_lung %>%
       filter(.imp==0))

Data_Oracle_lung_imp0 <-Data_Oracle_lung %>%
  filter(.imp==0)
sum(is.na(Data_Oracle_lung_imp0$FEV1_preBD_L_0W))
sum(is.na(Data_Oracle_lung_imp0$FEV1_preBD_L_52W))
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##CALCULATION OF DELTA FEV1 (calculated without imputated data)
Data_Oracle_lung <- Data_Oracle_lung %>%
  mutate(delta_FEV1_preBD_L= FEV1_preBD_L_52W-FEV1_preBD_L_0W)%>%
  mutate(delta_FEV1_preBD_mL= (FEV1_preBD_L_52W-FEV1_preBD_L_0W)*1000) %>%
  # Calculated with imputated data
  mutate(delta_FEV1_preBD_PCT= FEV1_preBD_PCT_52W-FEV1_preBD_PCT_0W)

colnames(Data_Oracle_lung)

#Put Treatment step in ordinal
Data_Oracle_lung$Treatment_step_ordinal <- as.numeric(gsub("Step ", "", as.character(Data_Oracle_lung$Treatment_step)))
summary(Data_Oracle_lung$Treatment_step_ordinal)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Dataset for the main analysis
Data_Oracle_lung_clean<- Data_Oracle_lung %>%
  filter(!is.na(delta_FEV1_preBD_mL))
Data_Oracle_lung_imputated <- Data_Oracle_lung_clean %>%
  filter(.imp != 0)
unique(Data_Oracle_lung_imputated$.imp)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##Preparation of different dataset - ACQ
#ACQ : Study selection

#unique(Raw_ACQ_FEV1$Enrolled_Trial_name)
#Studies_included_ACQ<- unique(Raw_ACQ_FEV1$Enrolled_Trial_name)
Studies_included_ACQ <- c("BENRAP2B","LAVOLTA_1","LAVOLTA_2", "LUTE", "VERSE", "MILLY", "STRATOS_1", "STRATOS_2")

Data_Oracle_ACQ<-Data_Oracle %>%
  filter(Enrolled_Trial_name%in%Studies_included_ACQ)

Studies_included_ACQ
unique(Data_Oracle_ACQ$Enrolled_Trial_name)
nrow(Data_Oracle_ACQ%>%
       filter(.imp==0))
Data_Oracle_ACQ_imp0 <-Data_Oracle_ACQ %>%
  filter(.imp==0)
sum(is.na(Data_Oracle_ACQ_imp0$ACQ_score_0W))
sum(is.na(Data_Oracle_ACQ_imp0$ACQ_score_24W))
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Preparation for analysis of ACQ placebo change
sum(is.na(Data_Oracle_imp1to10$ACQ_score_0W))

# Calculate ACQ change using mutate with case_when
Data_Oracle_ACQ <- Data_Oracle_ACQ %>%
  mutate(
    # Determine which follow-up measurement to use with case_when
    ACQ5_followup = case_when(
      !is.na(ACQ_score_24W) ~ ACQ_score_24W,
      !is.na(ACQ_score_25W) ~ ACQ_score_25W,
      !is.na(ACQ_score_26W) ~ ACQ_score_26W,
      TRUE ~ NA_real_  # If neither is available
    ),
    
    # Calculate change from baseline
    delta_ACQ = case_when(
      !is.na(ACQ5_followup) & !is.na(ACQ_score_0W) ~  ACQ_score_0W-ACQ5_followup,
      TRUE ~ NA_real_
    )
  )

#Calculate NA - Delta_ACQ
cols <- c(
  "ACQ_score_24W", 
  "ACQ_score_25W", 
  "ACQ_score_26W"
)

# Count NAs per column
sum(rowSums(is.na(Data_Oracle_ACQ[ , cols])) == 3)
nrow(Data_Oracle_ACQ)
sum(is.na(Data_Oracle_ACQ$delta_ACQ))



#Selection of patients with delta ACQ is not NA
#Data_Oracle_ACQ<-Data_Oracle_ACQ %>%
# filter(!is.na(delta_ACQ))

#Put Treatment step in ordinal
Data_Oracle_ACQ$Treatment_step_ordinal <- as.numeric(gsub("Step ", "", as.character(Data_Oracle_ACQ$Treatment_step)))
summary(Data_Oracle_ACQ$Treatment_step_ordinal)


# Count which follow-up week was actually used
# respecting the priority:
# 24W > 25W > 26W

week_used_count <- Data_Oracle_ACQ %>%
  filter(.imp==0) %>%
  mutate(
    week_used = case_when(
      !is.na(ACQ_score_24W) ~ "24W",
      is.na(ACQ_score_24W) & !is.na(ACQ_score_25W) ~ "25W",
      is.na(ACQ_score_24W) & 
        is.na(ACQ_score_25W) & 
        !is.na(ACQ_score_26W) ~ "26W",
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(
    !is.na(week_used),
    !is.na(ACQ_score_0W)
  ) %>%
  
  count(week_used) %>%
  
  mutate(
    pct = round(n / sum(n) * 100, 1)
  )

print(week_used_count)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Dataset for the main analysis
Data_Oracle_ACQ_clean<- Data_Oracle_ACQ %>%
  filter(!is.na(delta_ACQ))
Data_Oracle_ACQ_imputated <- Data_Oracle_ACQ %>%
  # Identify participants where Delta_attack is missing in .imp == 0
  group_by(Sequential_number) %>%
  filter(!any(.imp == 0 & is.na(delta_ACQ))) %>%
  ungroup() %>%
  # Keep only .imp = 1 to 10
  filter(.imp != 0)
unique(Data_Oracle_ACQ_imputated$.imp)

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================

# SUPPLEMENTARY FIGURES

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# TABLES S_
# Other tables can be downloaded through the following analyses (Figures S_) - below each analyses
# If not present, create the following Word Table displaying the statistics in a table.
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
##TABLE S1.A - TABLE 1 (Characteristics) - Per trial (Asthma attacks Analysis)
#Table S1. Characteristics included in Table 1 with stratification by trial.
#A. Asthma attack (Analysis A)
#B. FEV1 (Analysis B)
#C. ACQ-5 (Analysis C)
#======================================================================================================================================

table_per_trial_attacks <- tbl_summary(
  data_Table_attack,
  by = Enrolled_Trial_name,
  include = c(
    Age, Gender, BMI,
    Smoking_Statut, Atopy_history,
    Airborne_allergen_sensibilisation, Eczema,
    Allergic_rhinitis, CRSsNP, CRSwNP,
    Treatment_step, ACQ_baseline_score_mean,
    Any_severe_attack_previous_12m_0no_1yes,
    Attack_12mo_Nb, Number_severe_asthma_attacks_during_followup,
    maintenance_OCS_prescribed__0no_1yes,
    FEV1_preBD_L_Baseline, FEV1_preBD_PCT_Baseline,
    FEV1_PCT_reversibility_postBD, Blood_Eos_baseline_x10_9_cells_per_L, FeNO_baseline_ppb, Total_IgE,
    Follow_up_duration_days, FEV1PREBD_L_52W,
    FEV1PREBD_PCT_52W, FEV1_FVC_ratio
  ),
  statistic = list(
    # Categorical (multinomial)
    Gender ~ "{n} / {N} ({p}%)",
    Treatment_step ~ "{n} / {N} ({p}%)",
    
    # Dichotomous
    Smoking_Statut ~ "{n} / {N} ({p}%)",
    Atopy_history ~ "{n} / {N} ({p}%)",
    Airborne_allergen_sensibilisation ~ "{n} / {N} ({p}%)",
    Eczema ~ "{n} / {N} ({p}%)",
    Allergic_rhinitis ~ "{n} / {N} ({p}%)",
    CRSsNP ~ "{n} / {N} ({p}%)",
    CRSwNP ~ "{n} / {N} ({p}%)",
    Any_severe_attack_previous_12m_0no_1yes ~ "{n} / {N} ({p}%)",
    maintenance_OCS_prescribed__0no_1yes ~ "{n} / {N} ({p}%)",
    
    # Continuous — mean (IQR)
    Age ~ "{mean} ({p25}-{p75})",
    BMI ~ "{mean} ({p25}-{p75})",
    ACQ_baseline_score_mean ~ "{mean} ({p25}-{p75})",
    FEV1_preBD_L_Baseline ~ "{mean} ({p25}-{p75})",
    FEV1_preBD_PCT_Baseline ~ "{mean} ({p25}-{p75})",
    FEV1_PCT_reversibility_postBD ~ "{mean} ({p25}-{p75})",
    Blood_Eos_baseline_x10_9_cells_per_L ~ "{mean} ({p25}-{p75})",
    FeNO_baseline_ppb ~ "{mean} ({p25}-{p75})",
    Total_IgE ~ "{mean} ({p25}-{p75})",
    Follow_up_duration_days ~ "{mean} ({p25}-{p75})",
    FEV1_FVC_ratio ~ "{mean} ({p25}-{p75})",
    FEV1PREBD_L_52W ~ "{mean} ({p25}-{p75})",
    FEV1PREBD_PCT_52W ~ "{mean} ({p25}-{p75})",
    
    # Count variables (sum only)
    Number_severe_asthma_attacks_during_followup ~ "{sum}",
    Attack_12mo_Nb ~ "{sum}"
  ),
  missing_text = "(Missing)"
)

table_per_trial_attacks
gt_table <- as_gt(table_per_trial_attacks)
gtsave(gt_table, filename = "table_per_trial_attacks.pdf")

#-----------------------------------------------------------------------------------------------------------------------

##TABLE S1.B - TABLE 1 (Characteristics) - Per trial (LUNG ANALYSIS)
table_per_trial_lung <- tbl_summary(
  by = Enrolled_Trial_name,
  data_Table_lung,
  include = c(
    Age, Gender, BMI,
    Smoking_Statut, Atopy_history,
    Airborne_allergen_sensibilisation, Eczema,
    Allergic_rhinitis, CRSsNP, CRSwNP,
    Treatment_step, ACQ_baseline_score_mean,
    Any_severe_attack_previous_12m_0no_1yes,
    Attack_12mo_Nb, Number_severe_asthma_attacks_during_followup,
    maintenance_OCS_prescribed__0no_1yes,
    FEV1_preBD_L_Baseline, FEV1_preBD_PCT_Baseline,
    FEV1_PCT_reversibility_postBD, Blood_Eos_baseline_x10_9_cells_per_L, FeNO_baseline_ppb, Total_IgE,
    Follow_up_duration_days, FEV1PREBD_L_52W,
    FEV1PREBD_PCT_52W, FEV1_FVC_ratio, delta_FEV1_preBD_mL
  ),
  statistic = list(
    # Categorical (multinomial)
    Gender ~ "{n} / {N} ({p}%)",
    Treatment_step ~ "{n} / {N} ({p}%)",
    
    # Dichotomous
    Smoking_Statut ~ "{n} / {N} ({p}%)",
    Atopy_history ~ "{n} / {N} ({p}%)",
    Airborne_allergen_sensibilisation ~ "{n} / {N} ({p}%)",
    Eczema ~ "{n} / {N} ({p}%)",
    Allergic_rhinitis ~ "{n} / {N} ({p}%)",
    CRSsNP ~ "{n} / {N} ({p}%)",
    CRSwNP ~ "{n} / {N} ({p}%)",
    Any_severe_attack_previous_12m_0no_1yes ~ "{n} / {N} ({p}%)",
    maintenance_OCS_prescribed__0no_1yes ~ "{n} / {N} ({p}%)",
    
    # Continuous — mean (IQR)
    Age ~ "{mean} ({p25}-{p75})",
    BMI ~ "{mean} ({p25}-{p75})",
    ACQ_baseline_score_mean ~ "{mean} ({p25}-{p75})",
    FEV1_preBD_L_Baseline ~ "{mean} ({p25}-{p75})",
    FEV1_preBD_PCT_Baseline ~ "{mean} ({p25}-{p75})",
    FEV1_PCT_reversibility_postBD ~ "{mean} ({p25}-{p75})",
    Blood_Eos_baseline_x10_9_cells_per_L ~ "{mean} ({p25}-{p75})",
    FeNO_baseline_ppb ~ "{mean} ({p25}-{p75})",
    Total_IgE ~ "{mean} ({p25}-{p75})",
    Follow_up_duration_days ~ "{mean} ({p25}-{p75})",
    FEV1_FVC_ratio ~ "{mean} ({p25}-{p75})",
    FEV1PREBD_L_52W ~ "{mean} ({p25}-{p75})",
    FEV1PREBD_PCT_52W ~ "{mean} ({p25}-{p75})",
    
    # Count variables (sum only)
    Number_severe_asthma_attacks_during_followup ~ "{sum}",
    Attack_12mo_Nb ~ "{sum}",
    
    # Count variables (sum only)
    delta_FEV1_preBD_mL ~ "{mean} ({p25}-{p75}) ({sd})"
  ),
  missing_text = "(Missing)"
)

table_per_trial_lung

gt_table <- as_gt(table_per_trial_lung)
gtsave(gt_table, filename = "table_per_trial_lung.pdf")

#-----------------------------------------------------------------------------------------------------------------------

##TABLE S1.C - TABLE 1 (Characteristics) - Per trial (ACQ-5 ANALYSIS)

table_per_trial_ACQ <- tbl_summary(
  data_Table_ACQ,
  by = Enrolled_Trial_name,
  include = c(
    Age, Gender, BMI,
    Smoking_Statut, Atopy_history,
    Airborne_allergen_sensibilisation, Eczema,
    Allergic_rhinitis, CRSsNP, CRSwNP,
    Treatment_step, ACQ_baseline_score_mean,
    Any_severe_attack_previous_12m_0no_1yes,
    Attack_12mo_Nb, Number_severe_asthma_attacks_during_followup,
    maintenance_OCS_prescribed__0no_1yes,
    FEV1_preBD_L_Baseline, FEV1_preBD_PCT_Baseline,
    FEV1_PCT_reversibility_postBD, Blood_Eos_baseline_x10_9_cells_per_L, FeNO_baseline_ppb, Total_IgE,
    Follow_up_duration_days, FEV1PREBD_L_52W,
    FEV1PREBD_PCT_52W, FEV1_FVC_ratio, delta_ACQ
  ),
  statistic = list(
    # Categorical (multinomial)
    Gender ~ "{n} / {N} ({p}%)",
    Treatment_step ~ "{n} / {N} ({p}%)",
    
    # Dichotomous
    Smoking_Statut ~ "{n} / {N} ({p}%)",
    Atopy_history ~ "{n} / {N} ({p}%)",
    Airborne_allergen_sensibilisation ~ "{n} / {N} ({p}%)",
    Eczema ~ "{n} / {N} ({p}%)",
    Allergic_rhinitis ~ "{n} / {N} ({p}%)",
    CRSsNP ~ "{n} / {N} ({p}%)",
    CRSwNP ~ "{n} / {N} ({p}%)",
    Any_severe_attack_previous_12m_0no_1yes ~ "{n} / {N} ({p}%)",
    maintenance_OCS_prescribed__0no_1yes ~ "{n} / {N} ({p}%)",
    
    # Continuous — mean (IQR)
    Age ~ "{mean} ({p25}-{p75})",
    BMI ~ "{mean} ({p25}-{p75})",
    ACQ_baseline_score_mean ~ "{mean} ({p25}-{p75})",
    FEV1_preBD_L_Baseline ~ "{mean} ({p25}-{p75})",
    FEV1_preBD_PCT_Baseline ~ "{mean} ({p25}-{p75})",
    FEV1_PCT_reversibility_postBD ~ "{mean} ({p25}-{p75})",
    Blood_Eos_baseline_x10_9_cells_per_L ~ "{mean} ({p25}-{p75})",
    FeNO_baseline_ppb ~ "{mean} ({p25}-{p75})",
    Total_IgE ~ "{mean} ({p25}-{p75})",
    Follow_up_duration_days ~ "{mean} ({p25}-{p75})",
    FEV1_FVC_ratio ~ "{mean} ({p25}-{p75})",
    FEV1PREBD_L_52W ~ "{mean} ({p25}-{p75})",
    FEV1PREBD_PCT_52W ~ "{mean} ({p25}-{p75})",
    
    # Count variables (sum only)
    Number_severe_asthma_attacks_during_followup ~ "{sum}",
    Attack_12mo_Nb ~ "{sum}",
    
    # Count variables (sum only)
    delta_ACQ ~ "{mean} ({p25}-{p75}) ({sd})"
  ),
  missing_text = "(Missing)"
)

table_per_trial_ACQ

gt_table <- as_gt(table_per_trial_ACQ)
gtsave(gt_table, filename = "table_per_trial_ACQ.pdf")

# =============================================================================
# Export Table S1.A, S1.B, S1.C to Word (.docx)
# Requires: gtsummary, flextable, officer
# =============================================================================
# install.packages(c("gtsummary", "flextable", "officer"))

library(gtsummary)
library(flextable)
library(officer)

# =============================================================================
# Helper function: saves a tbl_summary object as a Word document
# =============================================================================

save_gtsummary_docx <- function(tbl, filename, title = NULL) {
  
  # Convert gtsummary → flextable
  ft <- as_flex_table(tbl)
  
  # --- Flextable styling ---
  ft <- ft |>
    # Font throughout
    flextable::font(fontname = "Arial", part = "all") |>
    flextable::fontsize(size = 9, part = "all") |>
    flextable::fontsize(size = 10, part = "header") |>
    
    # Header: bold + blue background + white text
    flextable::bold(part = "header") |>
    flextable::bg(bg = "#4472C4", part = "header") |>
    flextable::color(color = "white", part = "header") |>
    flextable::align(align = "center", part = "header") |>
    
    # Body: alternating row shading
    flextable::bg(i = seq(2, nrow(tbl$table_body), 2),
                  bg = "#F2F2F2", part = "body") |>
    
    # Borders
    flextable::border_outer(part = "all",
                            border = fp_border(color = "#AAAAAA", width = 0.5)) |>
    flextable::border_inner(part = "all",
                            border = fp_border(color = "#CCCCCC", width = 0.5)) |>
    
    # Cell padding
    flextable::padding(padding = 3, part = "all") |>
    
    # Fit to page width (landscape)
    flextable::set_table_properties(layout = "autofit", width = 1)
  
  # --- Word document (landscape) ---
  doc <- read_docx() |>
    body_set_default_section(
      prop_section(
        page_size = page_size(orient = "landscape",
                              width  = 11,   # inches
                              height = 8.5),
        page_margins = page_mar(top    = 0.5,
                                bottom = 0.5,
                                left   = 0.5,
                                right  = 0.5)
      )
    )
  
  # Optional title paragraph
  if (!is.null(title)) {
    doc <- doc |>
      body_add_par(title,
                   style = "heading 1") |>
      body_add_par("", style = "Normal")   # blank line
  }
  
  # Add table
  doc <- doc |>
    body_add_flextable(ft) |>
    body_add_par("", style = "Normal") |>  # blank line after table
    
    # Footnote / abbreviations
    body_add_par(
      paste(
        "Abbreviations: ACQ-5, Asthma Control Questionnaire 5-item;",
        "BD, bronchodilator; BMI, body mass index;",
        "CRSsNP, chronic rhinosinusitis without nasal polyps;",
        "CRSwNP, chronic rhinosinusitis with nasal polyps;",
        "FeNO, fractional exhaled nitric oxide;",
        "FEV1, forced expiratory volume in 1 second;",
        "FVC, forced vital capacity; IgE, immunoglobulin E;",
        "IQR, interquartile range (P25–P75);",
        "OCS, oral corticosteroids; SD, standard deviation."
      ),
      style = "Normal"
    )
  
  print(doc, target = filename)
  message("Saved: ", filename)
}

# =============================================================================
# TABLE S1.A — Asthma Attack Analysis
# =============================================================================

table_per_trial_attacks <- tbl_summary(
  data    = data_Table_attack,
  by      = Enrolled_Trial_name,
  include = c(
    Age, Gender, BMI,
    Smoking_Statut, Atopy_history,
    Airborne_allergen_sensibilisation, Eczema,
    Allergic_rhinitis, CRSsNP, CRSwNP,
    Treatment_step, ACQ_baseline_score_mean,
    Any_severe_attack_previous_12m_0no_1yes,
    Attack_12mo_Nb, Number_severe_asthma_attacks_during_followup,
    maintenance_OCS_prescribed__0no_1yes,
    FEV1_preBD_L_Baseline, FEV1_preBD_PCT_Baseline,
    FEV1_PCT_reversibility_postBD,
    Blood_Eos_baseline_x10_9_cells_per_L, FeNO_baseline_ppb, Total_IgE,
    Follow_up_duration_days, FEV1PREBD_L_52W,
    FEV1PREBD_PCT_52W, FEV1_FVC_ratio
  ),
  statistic = list(
    Gender                                        ~ "{n} / {N} ({p}%)",
    Treatment_step                                ~ "{n} / {N} ({p}%)",
    Smoking_Statut                                ~ "{n} / {N} ({p}%)",
    Atopy_history                                 ~ "{n} / {N} ({p}%)",
    Airborne_allergen_sensibilisation             ~ "{n} / {N} ({p}%)",
    Eczema                                        ~ "{n} / {N} ({p}%)",
    Allergic_rhinitis                             ~ "{n} / {N} ({p}%)",
    CRSsNP                                        ~ "{n} / {N} ({p}%)",
    CRSwNP                                        ~ "{n} / {N} ({p}%)",
    Any_severe_attack_previous_12m_0no_1yes       ~ "{n} / {N} ({p}%)",
    maintenance_OCS_prescribed__0no_1yes          ~ "{n} / {N} ({p}%)",
    Age                                           ~ "{mean} ({p25}-{p75})",
    BMI                                           ~ "{mean} ({p25}-{p75})",
    ACQ_baseline_score_mean                       ~ "{mean} ({p25}-{p75})",
    FEV1_preBD_L_Baseline                         ~ "{mean} ({p25}-{p75})",
    FEV1_preBD_PCT_Baseline                       ~ "{mean} ({p25}-{p75})",
    FEV1_PCT_reversibility_postBD                 ~ "{mean} ({p25}-{p75})",
    Blood_Eos_baseline_x10_9_cells_per_L          ~ "{mean} ({p25}-{p75})",
    FeNO_baseline_ppb                             ~ "{mean} ({p25}-{p75})",
    Total_IgE                                     ~ "{mean} ({p25}-{p75})",
    Follow_up_duration_days                       ~ "{mean} ({p25}-{p75})",
    FEV1_FVC_ratio                                ~ "{mean} ({p25}-{p75})",
    FEV1PREBD_L_52W                               ~ "{mean} ({p25}-{p75})",
    FEV1PREBD_PCT_52W                             ~ "{mean} ({p25}-{p75})",
    Number_severe_asthma_attacks_during_followup  ~ "{sum}",
    Attack_12mo_Nb                                ~ "{sum}"
  ),
  missing_text = "(Missing)"
) |>
  add_overall() |>
  bold_labels()

save_gtsummary_docx(
  tbl      = table_per_trial_attacks,
  filename = "table_per_trial_attacks.docx",
  title    = "Table S1.A — Baseline Characteristics per Trial: Asthma Attack Analysis (Analysis A)"
)

# =============================================================================
# TABLE S1.B — FEV1 / Lung Function Analysis
# =============================================================================

table_per_trial_lung <- tbl_summary(
  data    = data_Table_lung,
  by      = Enrolled_Trial_name,
  include = c(
    Age, Gender, BMI,
    Smoking_Statut, Atopy_history,
    Airborne_allergen_sensibilisation, Eczema,
    Allergic_rhinitis, CRSsNP, CRSwNP,
    Treatment_step, ACQ_baseline_score_mean,
    Any_severe_attack_previous_12m_0no_1yes,
    Attack_12mo_Nb, Number_severe_asthma_attacks_during_followup,
    maintenance_OCS_prescribed__0no_1yes,
    FEV1_preBD_L_Baseline, FEV1_preBD_PCT_Baseline,
    FEV1_PCT_reversibility_postBD,
    Blood_Eos_baseline_x10_9_cells_per_L, FeNO_baseline_ppb, Total_IgE,
    Follow_up_duration_days, FEV1PREBD_L_52W,
    FEV1PREBD_PCT_52W, FEV1_FVC_ratio, delta_FEV1_preBD_mL
  ),
  statistic = list(
    Gender                                        ~ "{n} / {N} ({p}%)",
    Treatment_step                                ~ "{n} / {N} ({p}%)",
    Smoking_Statut                                ~ "{n} / {N} ({p}%)",
    Atopy_history                                 ~ "{n} / {N} ({p}%)",
    Airborne_allergen_sensibilisation             ~ "{n} / {N} ({p}%)",
    Eczema                                        ~ "{n} / {N} ({p}%)",
    Allergic_rhinitis                             ~ "{n} / {N} ({p}%)",
    CRSsNP                                        ~ "{n} / {N} ({p}%)",
    CRSwNP                                        ~ "{n} / {N} ({p}%)",
    Any_severe_attack_previous_12m_0no_1yes       ~ "{n} / {N} ({p}%)",
    maintenance_OCS_prescribed__0no_1yes          ~ "{n} / {N} ({p}%)",
    Age                                           ~ "{mean} ({p25}-{p75})",
    BMI                                           ~ "{mean} ({p25}-{p75})",
    ACQ_baseline_score_mean                       ~ "{mean} ({p25}-{p75})",
    FEV1_preBD_L_Baseline                         ~ "{mean} ({p25}-{p75})",
    FEV1_preBD_PCT_Baseline                       ~ "{mean} ({p25}-{p75})",
    FEV1_PCT_reversibility_postBD                 ~ "{mean} ({p25}-{p75})",
    Blood_Eos_baseline_x10_9_cells_per_L          ~ "{mean} ({p25}-{p75})",
    FeNO_baseline_ppb                             ~ "{mean} ({p25}-{p75})",
    Total_IgE                                     ~ "{mean} ({p25}-{p75})",
    Follow_up_duration_days                       ~ "{mean} ({p25}-{p75})",
    FEV1_FVC_ratio                                ~ "{mean} ({p25}-{p75})",
    FEV1PREBD_L_52W                               ~ "{mean} ({p25}-{p75})",
    FEV1PREBD_PCT_52W                             ~ "{mean} ({p25}-{p75})",
    Number_severe_asthma_attacks_during_followup  ~ "{sum}",
    Attack_12mo_Nb                                ~ "{sum}",
    delta_FEV1_preBD_mL                           ~ "{mean} ({p25}-{p75}) ({sd})"
  ),
  missing_text = "(Missing)"
) |>
  add_overall() |>
  bold_labels()

save_gtsummary_docx(
  tbl      = table_per_trial_lung,
  filename = "table_per_trial_lung.docx",
  title    = "Table S1.B — Baseline Characteristics per Trial: FEV1 Analysis (Analysis B)"
)

# =============================================================================
# TABLE S1.C — ACQ-5 Analysis
# =============================================================================

table_per_trial_ACQ <- tbl_summary(
  data    = data_Table_ACQ,
  by      = Enrolled_Trial_name,
  include = c(
    Age, Gender, BMI,
    Smoking_Statut, Atopy_history,
    Airborne_allergen_sensibilisation, Eczema,
    Allergic_rhinitis, CRSsNP, CRSwNP,
    Treatment_step, ACQ_baseline_score_mean,
    Any_severe_attack_previous_12m_0no_1yes,
    Attack_12mo_Nb, Number_severe_asthma_attacks_during_followup,
    maintenance_OCS_prescribed__0no_1yes,
    FEV1_preBD_L_Baseline, FEV1_preBD_PCT_Baseline,
    FEV1_PCT_reversibility_postBD,
    Blood_Eos_baseline_x10_9_cells_per_L, FeNO_baseline_ppb, Total_IgE,
    Follow_up_duration_days, FEV1PREBD_L_52W,
    FEV1PREBD_PCT_52W, FEV1_FVC_ratio, delta_ACQ
  ),
  statistic = list(
    Gender                                        ~ "{n} / {N} ({p}%)",
    Treatment_step                                ~ "{n} / {N} ({p}%)",
    Smoking_Statut                                ~ "{n} / {N} ({p}%)",
    Atopy_history                                 ~ "{n} / {N} ({p}%)",
    Airborne_allergen_sensibilisation             ~ "{n} / {N} ({p}%)",
    Eczema                                        ~ "{n} / {N} ({p}%)",
    Allergic_rhinitis                             ~ "{n} / {N} ({p}%)",
    CRSsNP                                        ~ "{n} / {N} ({p}%)",
    CRSwNP                                        ~ "{n} / {N} ({p}%)",
    Any_severe_attack_previous_12m_0no_1yes       ~ "{n} / {N} ({p}%)",
    maintenance_OCS_prescribed__0no_1yes          ~ "{n} / {N} ({p}%)",
    Age                                           ~ "{mean} ({p25}-{p75})",
    BMI                                           ~ "{mean} ({p25}-{p75})",
    ACQ_baseline_score_mean                       ~ "{mean} ({p25}-{p75})",
    FEV1_preBD_L_Baseline                         ~ "{mean} ({p25}-{p75})",
    FEV1_preBD_PCT_Baseline                       ~ "{mean} ({p25}-{p75})",
    FEV1_PCT_reversibility_postBD                 ~ "{mean} ({p25}-{p75})",
    Blood_Eos_baseline_x10_9_cells_per_L          ~ "{mean} ({p25}-{p75})",
    FeNO_baseline_ppb                             ~ "{mean} ({p25}-{p75})",
    Total_IgE                                     ~ "{mean} ({p25}-{p75})",
    Follow_up_duration_days                       ~ "{mean} ({p25}-{p75})",
    FEV1_FVC_ratio                                ~ "{mean} ({p25}-{p75})",
    FEV1PREBD_L_52W                               ~ "{mean} ({p25}-{p75})",
    FEV1PREBD_PCT_52W                             ~ "{mean} ({p25}-{p75})",
    Number_severe_asthma_attacks_during_followup  ~ "{sum}",
    Attack_12mo_Nb                                ~ "{sum}",
    delta_ACQ                                     ~ "{mean} ({p25}-{p75}) ({sd})"
  ),
  missing_text = "(Missing)"
) |>
  add_overall() |>
  bold_labels()

save_gtsummary_docx(
  tbl      = table_per_trial_ACQ,
  filename = "table_per_trial_ACQ.docx",
  title    = "Table S1.C — Baseline Characteristics per Trial: ACQ-5 Analysis (Analysis C)"
)

#======================================================================================================================================
#======================================================================================================================================
# Figure S2 - Descriptive statistics summary of Figure 1
#======================================================================================================================================
#======================================================================================================================================

###TABLE WITH THE BOXPLOT (DESCRIPTIVE) STATISTICS
library(dplyr)

compute_stats <- function(data, outcome_col, trial_col = "Enrolled_Trial_name", analysis_label) {
  
  var <- sym(outcome_col)
  
  # Per trial
  per_trial <- data %>%
    filter(!is.na(!!var)) %>%
    group_by(Trial = !!sym(trial_col)) %>%
    summarise(
      N         = n(),
      Mean      = mean(!!var, na.rm = TRUE),
      Median    = median(!!var, na.rm = TRUE),
      SD        = sd(!!var, na.rm = TRUE),
      Q1        = quantile(!!var, 0.25, na.rm = TRUE),
      Q3        = quantile(!!var, 0.75, na.rm = TRUE),
      IQR       = IQR(!!var, na.rm = TRUE),
      Min       = min(!!var, na.rm = TRUE),
      Max       = max(!!var, na.rm = TRUE),
      CI_lower  = Mean - qt(0.975, df = N - 1) * (SD / sqrt(N)),
      CI_upper  = Mean + qt(0.975, df = N - 1) * (SD / sqrt(N)),
      .groups = "drop"
    ) %>%
    mutate(Analysis = analysis_label)
  
  # Pooled row
  pooled <- data %>%
    filter(!is.na(!!var)) %>%
    summarise(
      Trial     = "POOLED",
      N         = n(),
      Mean      = mean(!!var, na.rm = TRUE),
      Median    = median(!!var, na.rm = TRUE),
      SD        = sd(!!var, na.rm = TRUE),
      Q1        = quantile(!!var, 0.25, na.rm = TRUE),
      Q3        = quantile(!!var, 0.75, na.rm = TRUE),
      IQR       = IQR(!!var, na.rm = TRUE),
      Min       = min(!!var, na.rm = TRUE),
      Max       = max(!!var, na.rm = TRUE),
      CI_lower  = Mean - qt(0.975, df = N - 1) * (SD / sqrt(N)),
      CI_upper  = Mean + qt(0.975, df = N - 1) * (SD / sqrt(N)),
      Analysis  = analysis_label
    )
  
  bind_rows(per_trial, pooled)
}

# Run for each outcome (filter to .imp == 0)
stats_attack <- compute_stats(
  Data_Oracle_attack_clean %>% filter(.imp == 0),
  outcome_col = "Delta_attack",
  analysis_label = "Asthma Attacks"
)

stats_lung <- compute_stats(
  Data_Oracle_lung_clean_imp0,
  outcome_col = "delta_FEV1_preBD_mL",
  analysis_label = "FEV1 (mL)"
)

stats_acq <- compute_stats(
  Data_Oracle_ACQ_clean_imp0,
  outcome_col = "delta_ACQ",
  analysis_label = "ACQ-5"
)

# Combine all
stats_all <- bind_rows(stats_attack, stats_lung, stats_acq) %>%
  select(Analysis, Trial, N, Mean, Median, SD, Q1, Q3, IQR, Min, Max, CI_lower, CI_upper) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

# View
print(stats_all)

# Optional: export
write.csv(stats_all, "Descriptive_Stats_All_Outcomes.csv", row.names = FALSE)

# Optional: nice formatted table
library(gt)
stats_all %>%
  gt(groupname_col = "Analysis") %>%
  tab_header(title = "Descriptive Statistics by Trial and Outcome") %>%
  cols_label(
    Trial = "Trial", N = "n", Mean = "Mean", Median = "Median",
    SD = "SD", Q1 = "Q1", Q3 = "Q3", IQR = "IQR",
    Min = "Min", Max = "Max", CI_lower = "95% CI Lower", CI_upper = "95% CI Upper"
  )
library(flextable)
library(officer)

ft <- flextable(stats_all) %>%
  theme_vanilla() %>%
  autofit()

doc <- read_docx() %>%
  body_add_par("Descriptive Statistics by Trial and Outcome", style = "heading 1") %>%
  body_add_flextable(ft)

print(doc, 
      target = "/Users/joel/ORACLE - Placebo - Exacerbations/Project/Descriptive_Stats_All_Outcomes.docx")

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================

#FIGURES S_
#======================================================================================================================================
#======================================================================================================================================

#======================================================================================================================================
# (Figure S0 - Flowchart can be downloaded from PlaceboAsthma_MAIN_Script_clean.R)
#======================================================================================================================================

#======================================================================================================================================
# Figure S1 - Fraction of NA by trial by variable
#======================================================================================================================================
nrow(ORACLE_original)
colnames(ORACLE_original)
unique(ORACLE_Placebo_Attack_Raw$Enrolled_Trial_name)
unique(ORACLE_Placebo_Lung_Raw$Enrolled_Trial_name)
unique(ORACLE_Placebo_ACQ_Raw$Enrolled_Trial_name)

# =============================================================================
# Missing Data Heatmaps — Attack, Lung Function, and ACQ-5 analyses
# Based on ORACLE_original (placebo arm subsets)
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)
library(patchwork)
library(stringr)

# -----------------------------------------------------------------------------
# 0.  Variable mapping: display label → column name(s) in ORACLE_original
#     For variables with pre/post or multiple columns, we flag any row missing
#     across the relevant columns (any NA = missing for that participant).
# -----------------------------------------------------------------------------

base_vars <- list(
  list(label = "Age (years)",                         cols = "Age_con"),
  list(label = "Sex",                                 cols = "Sex"),
  list(label = "Body mass index (kg/m\u00b2)",        cols = "BMI_con"),
  list(label = "Country",                             cols = "Country"),
  list(label = "Smoking history",                     cols = "Smoking_history_yes_no"),
  list(label = "Psychiatric disease",                 cols = "Psychiatric_disease"),
  list(label = "Airborne allergen sensitisation",     cols = "Airborne_allergen_sensitisation"),
  list(label = "Allergic rhinitis",                   cols = "Allergic_Rhinitis"),
  list(label = "Eczema",                              cols = "Eczema"),
  list(label = "CRS without nasal polyps",            cols = "CRSsNP"),
  list(label = "CRS with nasal polyps",               cols = "CRSwNP"),
  list(label = "Intranasal steroid",                  cols = "Intranasal_steroid"),
  list(label = "mOCS",                                cols = "mOCS"),
  list(label = "Prior severe attacks (count)",        cols = "Attack_history_num"),
  list(label = "Prior hospitalisations (yes/no)",     cols = "Hospitalisation_yes_no"),
  list(label = "ICU admission or ETI history",        cols = "ICU_or_ETI_history"),
  list(label = "ACQ score (baseline)",                cols = "ACQ_score_0W"),
  list(label = "Pre-BD FEV1 (% predicted)",     cols = "FEV1_preBD_PCT_0W"),
  list(label = "Pre-BD FEV1/FVC",               cols = c("FEV1_preBD_L_0W", "FVC_preBD_L_0W")),
  list(label = "Post-BD FEV1 (% predicted)",    cols = "FEV1_postBD_PCT_0W"),
  list(label = "FEV\u2081 reversibility (baseline)",  cols = "FEV1_reversibility_0W"),
  list(label = "Blood eosinophils (log10)", cols = "BEC_log10"),
  list(label = "FeNO (log10)",              cols = "FeNO_log10"),
  list(label = "Total IgE (log10)",         cols = "IgE_log10"),
  list(label = "Follow-up duration (days)",       cols = "Log_Follow_up_duration_days")
)

attack_extra_vars <- list(
  list(label = "Severe attacks during follow-up (count)", cols = "Attack_number_during_followup")
)

lung_extra_vars <- list(
  list(label = "FEV1 pre-BD at 52W",             cols = "FEV1_preBD_PCT_52W"),
  list(label = "FEV1/FVC pre-BD at 52W",         cols = c("FEV1_preBD_L_52W", "FVC_preBD_L_52W")),
  list(label = "FEV1 post-BD at 52W",            cols = "FEV1_postBD_PCT_52W"),
  list(label = "FEV1 reversibility at 52W",      cols = "FEV1_reversibility_0W")   # replace with 52W col if available
)

acq_extra_vars <- list(
  list(label = "ACQ score at 24W", cols = "ACQ_score_24W"),
  list(label = "ACQ score at 25W", cols = "ACQ_score_25W"),
  list(label = "ACQ score at 26W", cols = "ACQ_score_26W")
)

# For Pre-BD FEV1/FVC: missingness = any NA in FEV1_preBD_L OR FVC_preBD_L
# For FEV1/FVC at 52W:  missingness = any NA in FEV1_preBD_L_52W OR FVC_preBD_L_52W

# -----------------------------------------------------------------------------
# 1.  Helper: compute % missing for a list of variable specs in a data frame
# -----------------------------------------------------------------------------

compute_missingness <- function(df, var_list) {
  
  # Replace NA-equivalent values with true NA across the whole sub-frame
  recode_na <- function(x) {
    if (length(na_equiv) > 0) {
      x[x %in% na_equiv] <- NA
    }
    x
  }
  
  results <- lapply(var_list, function(v) {
    cols <- v$cols
    cols_exist <- intersect(cols, colnames(df))
    if (length(cols_exist) == 0) {
      miss_pct <- NA_real_
    } else if (length(cols_exist) == 1) {
      vec <- recode_na(df[[cols_exist]])
      miss_pct <- mean(is.na(vec)) * 100
    } else {
      # multi-column: missing if ANY of the cols is NA
      # (i.e. the derived value, e.g. FEV1/FVC ratio, cannot be computed)
      sub <- as.data.frame(lapply(df[, cols_exist, drop = FALSE], recode_na))
      any_miss <- rowSums(is.na(sub)) > 0
      miss_pct <- mean(any_miss) * 100
    }
    data.frame(label = v$label, miss_pct = miss_pct, stringsAsFactors = FALSE)
  })
  bind_rows(results)
}

# -----------------------------------------------------------------------------
# 2.  Helper: build the full missingness table for one analysis
#     Returns a long data frame: trial (+ "Overall"), variable label, miss_pct
# -----------------------------------------------------------------------------

build_miss_table <- function(df, var_list, trial_col = "Enrolled_Trial_name") {
  
  trials <- sort(unique(as.character(df[[trial_col]])))
  n_overall <- nrow(df)
  
  # Per-trial
  per_trial <- lapply(trials, function(tr) {
    sub <- df[as.character(df[[trial_col]]) == tr, ]
    res <- compute_missingness(sub, var_list)
    res$trial <- tr
    res$n     <- nrow(sub)
    res
  })
  
  # Overall
  overall <- compute_missingness(df, var_list)
  overall$trial <- "Overall"
  overall$n     <- n_overall
  
  long <- bind_rows(c(per_trial, list(overall)))
  
  # Preserve variable order (top-to-bottom in plot = first to last in list)
  label_order <- sapply(var_list, function(v) v$label)
  long$label <- factor(long$label, levels = label_order)
  
  # Trial order: Overall first, then alphabetical trials
  trial_order <- c("Overall", trials)
  long$trial <- factor(long$trial, levels = trial_order)
  
  long
}

# -----------------------------------------------------------------------------
# 3.  Subset ORACLE_original by analysis
#     Adjust Treatment_arm filter to match your coding (e.g. "Placebo", "placebo", 0)
# -----------------------------------------------------------------------------

# --- Attack analysis ---
attack_trials <- c("AZISAST","BENRAP2B","DREAM","DRI12544","EXTRA",
                   "LAVOLTA_1","LAVOLTA_2","LUSTER_1","LUSTER_2",
                   "NAVIGATOR","PATHWAY","QUEST","STRATOS_1","STRATOS_2")

df_attack <- ORACLE_original %>%
  filter(Treatment_arm == "Placebo",              # adjust if needed
         Enrolled_Trial_name %in% attack_trials)

# --- Lung function analysis ---
lung_trials <- c("BENRAP2B","LAVOLTA_1","LAVOLTA_2","NAVIGATOR",
                 "PATHWAY","QUEST","STRATOS_1","STRATOS_2")

df_lung <- ORACLE_original %>%
  filter(Treatment_arm == "Placebo",
         Enrolled_Trial_name %in% lung_trials)

# --- ACQ-5 analysis ---
acq_trials <- c("BENRAP2B","LAVOLTA_1","LAVOLTA_2","LUTE",
                "MILLY","STRATOS_1","STRATOS_2","VERSE")

df_acq <- ORACLE_original %>%
  filter(Treatment_arm == "Placebo",
         Enrolled_Trial_name %in% acq_trials)

# -----------------------------------------------------------------------------
# 4.  Build missingness tables
# -----------------------------------------------------------------------------

miss_attack <- build_miss_table(df_attack,
                                c(base_vars, attack_extra_vars))

miss_lung   <- build_miss_table(df_lung,
                                c(base_vars, lung_extra_vars))

miss_acq    <- build_miss_table(df_acq,
                                c(base_vars, acq_extra_vars))

# -----------------------------------------------------------------------------
# 5.  Plotting function
# -----------------------------------------------------------------------------

make_heatmap <- function(long_df, title_text,
                         base_text_size = 8,
                         cell_text_size = 2.2) {
  
  # Build axis labels with n underneath
  trial_labels <- long_df %>%
    distinct(trial, n) %>%
    arrange(trial) %>%
    mutate(axis_label = ifelse(
      trial == "Overall",
      paste0("Overall\n(n=", formatC(n, big.mark = ",", format = "d"), ")"),
      paste0(trial, "\n(n=", formatC(n, big.mark = ",", format = "d"), ")")
    ))
  
  label_map <- setNames(trial_labels$axis_label, trial_labels$trial)
  long_df$trial_label <- factor(label_map[as.character(long_df$trial)],
                                levels = label_map[levels(long_df$trial)])
  
  # Identify which labels are italic (extra variables)
  all_labels <- levels(long_df$label)
  base_labels_vec <- sapply(base_vars, function(v) v$label)
  is_extra <- !all_labels %in% base_labels_vec
  
  face_vec <- ifelse(is_extra, "italic", "plain")
  
  ggplot(long_df, aes(x = trial_label, y = label, fill = miss_pct)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(miss_pct), "N/A",
                                 paste0(round(miss_pct, 1), "%"))),
              size = cell_text_size, color = ifelse(long_df$miss_pct > 55, "white", "grey20"),
              fontface = "plain") +
    scale_fill_gradient(
      low  = "#ffffff",
      high = "#7B0000",
      na.value = "grey85",
      limits = c(0, 100),
      name  = "% missing",
      guide = guide_colorbar(
        title.position = "top",
        barwidth  = unit(0.5, "cm"),
        barheight = unit(6,   "cm"),
        ticks.colour = "grey40"
      )
    ) +
    # Bold separator after Overall column
    geom_vline(xintercept = 1.5, color = "grey40", linewidth = 0.6) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(levels(long_df$label))) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = base_text_size) +
    theme(
      plot.title        = element_text(size = base_text_size + 2, face = "bold",
                                       margin = margin(b = 6)),
      axis.text.x.top   = element_text(angle = 45, hjust = 0, vjust = 0,
                                       size = base_text_size - 0.5,
                                       margin = margin(b = 2)),
      axis.text.y       = element_text(size = base_text_size - 0.5,
                                       face = face_vec,
                                       hjust = 1),
      panel.grid        = element_blank(),
      legend.title      = element_text(size = base_text_size - 1),
      legend.text       = element_text(size = base_text_size - 1.5),
      legend.position   = "right",
      plot.margin       = margin(t = 4, r = 4, b = 4, l = 4)
    )
}

# -----------------------------------------------------------------------------
# 6.  Generate plots
# -----------------------------------------------------------------------------

p_attack <- make_heatmap(miss_attack, "Attack analysis")
p_lung   <- make_heatmap(miss_lung,   "Lung function analysis")
p_acq    <- make_heatmap(miss_acq,    "ACQ-5 analysis")

# -----------------------------------------------------------------------------
# 7.  Export — individual files + combined figure
# -----------------------------------------------------------------------------

ggsave(file.path(Figure_path, "heatmap_attack.pdf"),   p_attack, width = 14, height = 10, units = "in")
ggsave(file.path(Figure_path, "heatmap_lung.pdf"),     p_lung,   width = 11, height = 10, units = "in")
ggsave(file.path(Figure_path, "heatmap_acq.pdf"),      p_acq,    width = 11, height = 10, units = "in")

# Combined stacked figure (useful for supplement)
combined <- p_attack / p_lung / p_acq +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(Figure_path, "heatmap_combined.pdf"), combined, width = 14, height = 30, units = "in")

message("Done. Four PDF files written to: ", Figure_path)
message("  heatmap_attack.pdf")
message("  heatmap_lung.pdf")
message("  heatmap_acq.pdf")
message("  heatmap_combined.pdf")

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S2 - Distribution of placebo change - Bar graph
#======================================================================================================================================
# For Delta Asthma Attacks

#FIGURE S2.A - Bar graph

subset_1 <- Data_Oracle_attack %>%
  filter(.imp == 0, !is.na(Delta_attack))

pdf("Figure_S2A_Histogram_DeltaAsthmaAttacks.pdf",
    width = 10, height = 8)

ggplot(subset_1, aes(x = Delta_attack_rate)) +
  geom_histogram(
    binwidth = 1,
    fill = "steelblue",
    color = "black",
    alpha = 0.7
  ) +
  scale_x_continuous(breaks = seq(-50, 15, 2)) +
  labs(
    title = "Distribution of Placebo Attack",
    x = "Placebo change (∆Attack)",
    y = "Count"
  )

dev.off()

#-----------------------------------------------------------------------------------------------------------------------
# For Delta FEV1

#FIGURE S2.B - Bar graph

subset_2 <- Data_Oracle_lung %>%
  filter(.imp == 0, !is.na(delta_FEV1_preBD_mL))

pdf("Figure_S2B_Histogram_DeltaFEV1.pdf",
    width = 10, height = 8)

ggplot(subset_2, aes(x = delta_FEV1_preBD_mL)) +
  geom_histogram(
    binwidth = 10,
    fill = "steelblue",
    color = "black",
    alpha = 0.7
  ) +
  scale_x_continuous(breaks = seq(-2000, 2000, 500)) +
  labs(
    title = "Distribution of Lung Placebo Change",
    x = "Placebo change (∆FEV1, mL)",
    y = "Count"
  )

dev.off()

#-----------------------------------------------------------------------------------------------------------------------
# For Delta ACQ-5

#FIGURE S2.C - Bar graph

subset_3 <- Data_Oracle_ACQ %>%
  filter(.imp == 0, !is.na(delta_ACQ))

pdf("Figure_S2C_Histogram_DeltaACQ.pdf",
    width = 10, height = 8)

ggplot(subset_3, aes(x = delta_ACQ)) +
  geom_histogram(
    binwidth = 0.5,
    fill = "steelblue",
    color = "black",
    alpha = 0.7
  ) +
  scale_x_continuous(breaks = seq(-5, 6, 0.5)) +
  labs(
    title = "Distribution of ACQ Placebo Change",
    x = "Placebo change (∆ACQ)",
    y = "Count"
  )

dev.off()

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S3 - Q-Q PLOTS WITH STATISTICAL ANNOTATIONS
#======================================================================================================================================

library(ggplot2)
library(e1071)
library(dplyr)

base_path <- "/Users/joel/ORACLE - Placebo - Exacerbations/Project/"

# === HELPER — compute stats for annotation ===
get_qq_stats <- function(x) {
  x <- x[!is.na(x)]
  list(
    n        = length(x),
    mean     = round(mean(x), 3),
    sd       = round(sd(x), 3),
    skew     = round(skewness(x), 3),
    kurt     = round(kurtosis(x), 3),
    shapiro_p = ifelse(
      length(x) > 5000,
      round(shapiro.test(sample(x, 5000))$p.value, 4),
      round(shapiro.test(x)$p.value, 4)
    )
  )
}

# === HELPER — build annotation label ===
make_label <- function(s) {
  paste0(
    "n = ",        s$n,        "\n",
    "Mean = ",     s$mean,     "\n",
    "SD = ",       s$sd,       "\n",
    "Skewness = ", s$skew,     "\n",
    "Kurtosis = ", s$kurt,     "\n",
    "Shapiro p ",  ifelse(s$shapiro_p < 0.001, "< 0.001",
                          paste0("= ", s$shapiro_p))
  )
}

# === HELPER — build Q-Q plot ===
make_qq_plot <- function(data, var, title, axis_label) {
  
  s     <- get_qq_stats(data[[var]])
  label <- make_label(s)
  
  # Theoretical quantile range for annotation placement
  theo_max <- qnorm(0.995)
  
  ggplot(data, aes(sample = .data[[var]])) +
    
    stat_qq(
      color = "steelblue",
      size  = 0.8,
      alpha = 0.5
    ) +
    
    stat_qq_line(
      color     = "red",
      linetype  = "dashed",
      linewidth = 0.8
    ) +
    
    # Statistical annotation box
    annotate(
      "label",
      x          = -Inf,
      y          = Inf,
      label      = label,
      hjust      = -0.1,
      vjust      = 1.1,
      size       = 3.5,
      fill       = "white",
      color      = "grey30",
      label.size = 0.3,
      fontface   = "plain",
      family     = "mono"
    ) +
    
    labs(
      title    = title,
      subtitle = paste0("Imputation 1 (.imp == 0)"),
      x        = "Theoretical quantiles",
      y        = paste0("Sample quantiles — ", axis_label)
    ) +
    
    theme_minimal(base_size = 13) +
    theme(
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 10, color = "grey50", hjust = 0.5),
      axis.title.x     = element_text(size = 12, face = "bold"),
      axis.title.y     = element_text(size = 12, face = "bold"),
      axis.text        = element_text(size = 11),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background  = element_rect(fill = "white", color = "white")
    )
}

#======================================================================================================================================
# FIGURE S3.A — Delta Asthma Attack Rate
#======================================================================================================================================

subset_attack <- Data_Oracle_attack %>%
  filter(.imp == 0, !is.na(Delta_attack))

Fig_E2A <- make_qq_plot(
  data       = subset_attack,
  var        = "Delta_attack",
  title      = "Q-Q Plot — \u2206 Asthma Attacks (Placebo)",
  axis_label = "\u2206 Asthma attacks"
)

ggsave(
  filename = paste0(base_path, "Figure_S3A_QQplot_DeltaAsthmaAttacks.pdf"),
  plot     = Fig_E2A,
  width    = 7,
  height   = 7,
  dpi      = 300
)

ggsave(
  filename = paste0(base_path, "Figure_S3A_QQplot_DeltaAsthmaAttacks.png"),
  plot     = Fig_E2A,
  width    = 7,
  height   = 7,
  dpi      = 300
)

cat("✓ Figure S3A saved.\n")

#======================================================================================================================================
# FIGURE S3.B — Delta FEV1
#======================================================================================================================================

subset_lung <- Data_Oracle_lung %>%
  filter(.imp == 0, !is.na(delta_FEV1_preBD_mL))

Fig_E2B <- make_qq_plot(
  data       = subset_lung,
  var        = "delta_FEV1_preBD_mL",
  title      = "Q-Q Plot — \u2206 FEV1 pre-BD (Placebo)",
  axis_label = "\u2206 FEV1 (mL)"
)

ggsave(
  filename = paste0(base_path, "Figure_S3B_QQplot_DeltaFEV1.pdf"),
  plot     = Fig_E2B,
  width    = 7,
  height   = 7,
  dpi      = 300
)

ggsave(
  filename = paste0(base_path, "Figure_S3B_QQplot_DeltaFEV1.png"),
  plot     = Fig_E2B,
  width    = 7,
  height   = 7,
  dpi      = 300
)

cat("✓ Figure S3B saved.\n")

#======================================================================================================================================
# FIGURE S3.C — Delta ACQ-5
#======================================================================================================================================

subset_acq <- Data_Oracle_ACQ %>%
  filter(.imp == 0, !is.na(delta_ACQ))

Fig_E2C <- make_qq_plot(
  data       = subset_acq,
  var        = "delta_ACQ",
  title      = "Q-Q Plot — \u2206 ACQ-5 (Placebo)",
  axis_label = "\u2206 ACQ-5"
)

ggsave(
  filename = paste0(base_path, "Figure_S2C_QQplot_DeltaACQ.pdf"),
  plot     = Fig_E2C,
  width    = 7,
  height   = 7,
  dpi      = 300
)

ggsave(
  filename = paste0(base_path, "Figure_S2C_QQplot_DeltaACQ.png"),
  plot     = Fig_E2C,
  width    = 7,
  height   = 7,
  dpi      = 300
)

cat("✓ Figure S3C saved.\n")

#======================================================================================================================================
# COMBINED FIGURE S3 — All 3 panels side by side
#======================================================================================================================================

library(patchwork)

Fig_E2_combined <- (Fig_E2A | Fig_E2B | Fig_E2C) +
  plot_annotation(
    title   = "Figure S3 — Normality Assessment of Placebo Change Distributions",
    caption = "Red dashed line = theoretical normal distribution. 
               Shapiro-Wilk test performed on imputation 0 (.imp == 0).
               For n > 5000, Shapiro-Wilk performed on a random sample of 5000 observations.",
    theme   = theme(
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.caption = element_text(size =  9, color = "grey50", hjust = 0.5)
    )
  )

ggsave(
  filename = paste0(base_path, "Figure_S3_QQplots_Combined.pdf"),
  plot     = Fig_E2_combined,
  width    = 21,
  height   = 8,
  dpi      = 300
)

ggsave(
  filename = paste0(base_path, "Figure_S3_QQplots_Combined.png"),
  plot     = Fig_E2_combined,
  width    = 21,
  height   = 8,
  dpi      = 300
)

cat("✓ Figure S3 combined saved.\n")


#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S4 - SEE "Simple model with only Attack_history_num as covariable" (PART A) in the PlaceboAsthma_MAIN_Script_clean
#======================================================================================================================================

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S5 - Asthma attack - Forest plot’s results per trial + I2 (for FeNO and BEC)
#======================================================================================================================================

#======================================================================================================================================
# Figure S5.A - For FeNO (across trial)
#======================================================================================================================================
# Per-trial forest plot: FeNO_log10 → Delta_attack
# Uses MI_estimates() filtered per trial + full-dataset pooled overall
# Adjusted for: Attack_history_num | Random intercept: No (per-trial)
# =============================================================================

library(dplyr)
library(ggplot2)
library(metafor)

# =============================================================================
# 0. SETTINGS
# =============================================================================

data         <- ORACLE_Placebo_Attack_imputed
outcome      <- "Placebo_response_Attack"
predictor    <- "FeNO_log10"
covariate    <- "Attack_history_num"
trial_col    <- "Enrolled_Trial_name"
imp_col      <- ".imp"
followup_col <- "Follow_up_duration_days"


# =============================================================================
# 1. GET TRIAL NAMES
# =============================================================================

trials <- data[[trial_col]][data[[imp_col]] == min(data[[imp_col]])] |>
  unique() |> sort()
cat("Trials:", paste(trials, collapse = ", "), "\n\n")


# =============================================================================
# 2. RUN MI_estimates PER TRIAL
# =============================================================================

run_trial <- function(trial_name) {
  cat("Fitting:", trial_name, "... ")
  
  df <- data[data[[trial_col]] == trial_name, ]
  n  <- nrow(df[df[[imp_col]] == min(df[[imp_col]]), ])
  
  res <- tryCatch(
    MI_estimates(
      data             = df,
      outcome_var      = outcome,
      predictor_vars   = predictor,
      covariables      = covariate,
      imp_col          = imp_col,
      followup_offset  = "Yes",
      followup_col     = followup_col,
      model_type       = "lm"
    ),
    error = function(e) { message("ERROR: ", e$message); NULL }
  )
  if (is.null(res)) return(NULL)
  
  # MI_estimates returns a data.frame when predictor_vars has 1 element
  # Extract the FeNO row
  row <- res[grepl(predictor, res$term), ]
  if (nrow(row) == 0) { message("  predictor row not found"); return(NULL) }
  
  cat("OK (n =", n, ")\n")
  data.frame(
    trial   = trial_name,
    n       = n,
    est     = row$estimate[1],
    se      = row$std.error[1],
    lo95    = row$`2.5 %`[1],
    hi95    = row$`97.5 %`[1],
    p.value = row$p.value[1],
    stringsAsFactors = FALSE
  )
}

results_list <- lapply(trials, run_trial)
results      <- do.call(rbind, Filter(Negate(is.null), results_list))


# =============================================================================
# 3. POOLED OVERALL — full dataset, random intercept on trial
# =============================================================================

cat("\nFitting overall pooled model...\n")

overall_fit <- MI_estimates(
  data                 = data,
  outcome_var          = outcome,
  predictor_vars       = predictor,
  covariables          = covariate,
  imp_col              = imp_col,
  followup_offset      = "Yes",
  followup_col         = followup_col,
  random_intercept_var = trial_col,
  model_type           = "lm"
)

overall_row <- overall_fit[grepl(predictor, overall_fit$term), ]

overall <- data.frame(
  trial   = "Overall",
  n       = sum(results$n),
  est     = overall_row$estimate[1],
  se      = overall_row$std.error[1],
  lo95    = overall_row$`2.5 %`[1],
  hi95    = overall_row$`97.5 %`[1],
  p.value = overall_row$p.value[1],
  stringsAsFactors = FALSE
)

cat(sprintf("Overall: \u03b2 = %.3f [%.3f, %.3f]  p = %s\n",
            overall$est, overall$lo95, overall$hi95,
            ifelse(overall$p.value < 0.001, "<0.001",
                   sprintf("%.3f", overall$p.value))))


# =============================================================================
# 4. HETEROGENEITY (Q-based I², Borenstein H method)
# =============================================================================

meta   <- rma(yi = est, sei = se, data = results, method = "REML")
k_meta <- nrow(results)
w_meta <- 1 / results$se^2
Q_meta <- sum(w_meta * (results$est - weighted.mean(results$est, w_meta))^2)
df     <- k_meta - 1
i2_est <- max(0, (Q_meta - df) / Q_meta)

H      <- sqrt(max(Q_meta, df) / df)
lnH    <- log(H)
se_lnH <- if (Q_meta > df + 1) (log(Q_meta) - log(df)) / (2*(Q_meta - df)) else sqrt(1/(2*df))
H_lo   <- exp(lnH - 1.96 * se_lnH)
H_hi   <- exp(lnH + 1.96 * se_lnH)
i2_lo  <- max(0,        1 - 1/H_lo^2)
i2_hi  <- max(0, min(1, 1 - 1/H_hi^2))

i2_label <- if (Q_meta <= df) {
  sprintf("Heterogeneity: I\u00b2 = 0 (Q(%d) = %.2f, p = %s) \u2014 no between-trial variance detected",
          df, Q_meta,
          ifelse(meta$QEp < 0.001, "<0.001", sprintf("%.3f", meta$QEp)))
} else {
  sprintf("Heterogeneity: I\u00b2 = %.2f [%.2f, %.2f]  |  Q(%d) = %.2f, p = %s",
          i2_est, i2_lo, i2_hi, df, Q_meta,
          ifelse(meta$QEp < 0.001, "<0.001", sprintf("%.3f", meta$QEp)))
}
cat("\n", i2_label, "\n")


# =============================================================================
# 5. PREPARE PLOT DATA
# =============================================================================

# Alphabetical, A at top → reverse for ggplot y-axis (bottom to top)
results <- results |>
  arrange(trial) |>
  mutate(
    weight_pct = n / sum(n) * 100,
    sig        = p.value < 0.05,
    p_label    = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    label      = sprintf("%.2f [%.2f, %.2f]   p=%s   W=%.1f%%",
                         est, lo95, hi95, p_label, weight_pct),
    trial      = factor(trial, levels = rev(trial))
  )

K         <- nrow(results)
trial_y   <- seq(K, 1)
overall_y <- -1

plot_trials <- results |> mutate(y = trial_y)

# Dot size proportional to weight (manual, no size legend)
min_pt <- 2; max_pt <- 7
w_min  <- min(results$weight_pct); w_max <- max(results$weight_pct)
plot_trials <- plot_trials |>
  mutate(pt_size = min_pt + (weight_pct - w_min) / (w_max - w_min) * (max_pt - min_pt))

overall_label <- sprintf("%.2f [%.2f, %.2f]   p=%s",
                         overall$est, overall$lo95, overall$hi95,
                         ifelse(overall$p.value < 0.001, "<0.001",
                                sprintf("%.3f", overall$p.value)))

x_lo    <- min(c(results$lo95, overall$lo95)) - 0.1
x_hi    <- max(c(results$hi95, overall$hi95)) + 0.1
label_x <- x_hi + (x_hi - x_lo) * 0.02

y_breaks <- c(trial_y, overall_y)
y_labels <- c(as.character(results$trial), "Overall")


# =============================================================================
# 6. FOREST PLOT
# =============================================================================

forest_plot <- ggplot() +
  
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  
  geom_hline(yintercept = (overall_y + 1) / 2,
             colour = "grey60", linewidth = 0.35) +
  
  # CI lines
  geom_errorbarh(
    data = plot_trials,
    aes(y = y, xmin = lo95, xmax = hi95),
    height = 0.18, linewidth = 0.5, colour = "grey30"
  ) +
  
  # Points (manual size, no legend)
  geom_point(
    data  = plot_trials,
    aes(x = est, y = y, colour = sig, fill = sig),
    size  = plot_trials$pt_size,
    shape = 21, stroke = 0.55
  ) +
  
  # Overall diamond
  geom_polygon(
    data = data.frame(
      x = c(overall$lo95, overall$est, overall$hi95, overall$est),
      y = c(overall_y, overall_y + 0.42, overall_y, overall_y - 0.42)
    ),
    aes(x = x, y = y),
    fill = "#185FA5", colour = "#0D3D6E", linewidth = 0.5
  ) +
  
  # Trial labels (right)
  geom_text(
    data = plot_trials,
    aes(x = label_x, y = y, label = label, colour = sig),
    hjust = 0, size = 2.7, show.legend = FALSE
  ) +
  
  # Overall label (right)
  annotate("text", x = label_x, y = overall_y, label = overall_label,
           hjust = 0, size = 2.7, colour = "#185FA5", fontface = "bold") +
  
  # I² annotation
  annotate("text",
           x = (x_lo + x_hi) / 2, y = overall_y - 0.85,
           label = i2_label,
           size = 2.5, colour = "grey35", hjust = 0.5, fontface = "italic") +
  
  # Column header
  annotate("text", x = label_x, y = K + 0.85,
           label = "\u03b2 [95% CI]   p-value   Weight",
           hjust = 0, size = 2.7, colour = "grey30", fontface = "bold") +
  
  scale_colour_manual(
    values = c("TRUE" = "#A32D2D", "FALSE" = "#4A4A47"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "p \u2265 0.05"),
    name   = NULL
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#D96060", "FALSE" = "#AAAAAA"),
    guide  = "none"
  ) +
  scale_x_continuous(
    limits = c(x_lo, x_hi + (x_hi - x_lo) * 0.52),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = y_labels,
    expand = expansion(add = c(1.4, 1.2))
  ) +
  
  labs(
    title    = "Effect of FeNO (log10) on Delta Attack rate - Per trial",
    subtitle = paste0(
      "Univariate linear regression, adjusted for ", covariate
    ),
    x = "Adjusted coefficient for FeNO (log10)  (95% CI)",
    y = NULL
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    panel.border       = element_rect(colour = "grey70", linewidth = 0.5),
    axis.text.y        = element_text(size = 15),
    axis.text.x        = element_text(size = 15),
    plot.title         = element_text(face = "bold", size = 20),
    plot.subtitle      = element_text(size = 8.5, colour = "grey40"),
    legend.position    = "bottom",
    legend.text        = element_text(size = 12),
    plot.margin        = margin(10, 10, 10, 10)
  )

print(forest_plot)


# =============================================================================
# 7. SAVE
# =============================================================================

plot_height <- max(6, K * 0.55 + 3.5)

ggsave(file.path(getwd(), "Forest_plot_FeNO_per_trial_MIestimates.pdf"),
       plot = forest_plot, width = 18, height = plot_height, units = "in")

ggsave(file.path(getwd(), "Forest_plot_FeNO_per_trial_MIestimates.png"),
       plot = forest_plot, width = 18, height = plot_height,
       dpi = 300, units = "in")

cat("\nSaved to:", getwd(), "\n")


# =============================================================================
# 8. RETURN
# =============================================================================

list(results = results, overall = overall,
     I2 = i2_est, I2_ci = c(i2_lo, i2_hi), i2_label = i2_label)


#======================================================================================================================================
#======================================================================================================================================
# Figure S5.B - For BEC (across trial)
#======================================================================================================================================
# =============================================================================
# =============================================================================
# Per-trial forest plot: BEC_log10 → Delta_attack
# Uses MI_estimates() filtered per trial + full-dataset pooled overall
# Adjusted for: Attack_history_num | Random intercept: No (per-trial)
# =============================================================================

library(dplyr)
library(ggplot2)
library(metafor)


# =============================================================================
# 0. SETTINGS
# =============================================================================

data         <- ORACLE_Placebo_Attack_imputed
outcome      <- "Placebo_response_Attack"
predictor    <- "BEC_log10"
covariate    <- "Attack_history_num"
trial_col    <- "Enrolled_Trial_name"
imp_col      <- ".imp"
followup_col <- "Follow_up_duration_days"


# =============================================================================
# 1. GET TRIAL NAMES
# =============================================================================

trials <- data[[trial_col]][data[[imp_col]] == min(data[[imp_col]])] |>
  unique() |> sort()
cat("Trials:", paste(trials, collapse = ", "), "\n\n")


# =============================================================================
# 2. RUN MI_estimates PER TRIAL
# =============================================================================

run_trial <- function(trial_name) {
  cat("Fitting:", trial_name, "... ")
  
  df <- data[data[[trial_col]] == trial_name, ]
  n  <- nrow(df[df[[imp_col]] == min(df[[imp_col]]), ])
  
  res <- tryCatch(
    MI_estimates(
      data             = df,
      outcome_var      = outcome,
      predictor_vars   = predictor,
      covariables      = covariate,
      imp_col          = imp_col,
      followup_offset  = "Yes",
      followup_col     = followup_col,
      model_type       = "lm"
    ),
    error = function(e) { message("ERROR: ", e$message); NULL }
  )
  if (is.null(res)) return(NULL)
  
  # MI_estimates returns a data.frame when predictor_vars has 1 element
  # Extract the FeNO row
  row <- res[grepl(predictor, res$term), ]
  if (nrow(row) == 0) { message("  predictor row not found"); return(NULL) }
  
  cat("OK (n =", n, ")\n")
  data.frame(
    trial   = trial_name,
    n       = n,
    est     = row$estimate[1],
    se      = row$std.error[1],
    lo95    = row$`2.5 %`[1],
    hi95    = row$`97.5 %`[1],
    p.value = row$p.value[1],
    stringsAsFactors = FALSE
  )
}

results_list <- lapply(trials, run_trial)
results      <- do.call(rbind, Filter(Negate(is.null), results_list))


# =============================================================================
# 3. POOLED OVERALL — full dataset, random intercept on trial
# =============================================================================

cat("\nFitting overall pooled model...\n")

overall_fit <- MI_estimates(
  data                 = data,
  outcome_var          = outcome,
  predictor_vars       = predictor,
  covariables          = covariate,
  imp_col              = imp_col,
  followup_offset      = "Yes",
  followup_col         = followup_col,
  random_intercept_var = trial_col,
  model_type           = "lm"
)

overall_row <- overall_fit[grepl(predictor, overall_fit$term), ]

overall <- data.frame(
  trial   = "Overall",
  n       = sum(results$n),
  est     = overall_row$estimate[1],
  se      = overall_row$std.error[1],
  lo95    = overall_row$`2.5 %`[1],
  hi95    = overall_row$`97.5 %`[1],
  p.value = overall_row$p.value[1],
  stringsAsFactors = FALSE
)

cat(sprintf("Overall: \u03b2 = %.3f [%.3f, %.3f]  p = %s\n",
            overall$est, overall$lo95, overall$hi95,
            ifelse(overall$p.value < 0.001, "<0.001",
                   sprintf("%.3f", overall$p.value))))


# =============================================================================
# 4. HETEROGENEITY (Q-based I², Borenstein H method)
# =============================================================================

meta   <- rma(yi = est, sei = se, data = results, method = "REML")
k_meta <- nrow(results)
w_meta <- 1 / results$se^2
Q_meta <- sum(w_meta * (results$est - weighted.mean(results$est, w_meta))^2)
df     <- k_meta - 1
i2_est <- max(0, (Q_meta - df) / Q_meta)

H      <- sqrt(max(Q_meta, df) / df)
lnH    <- log(H)
se_lnH <- if (Q_meta > df + 1) (log(Q_meta) - log(df)) / (2*(Q_meta - df)) else sqrt(1/(2*df))
H_lo   <- exp(lnH - 1.96 * se_lnH)
H_hi   <- exp(lnH + 1.96 * se_lnH)
i2_lo  <- max(0,        1 - 1/H_lo^2)
i2_hi  <- max(0, min(1, 1 - 1/H_hi^2))

i2_label <- if (Q_meta <= df) {
  sprintf("Heterogeneity: I\u00b2 = 0 (Q(%d) = %.2f, p = %s) \u2014 no between-trial variance detected",
          df, Q_meta,
          ifelse(meta$QEp < 0.001, "<0.001", sprintf("%.3f", meta$QEp)))
} else {
  sprintf("Heterogeneity: I\u00b2 = %.2f [%.2f, %.2f]  |  Q(%d) = %.2f, p = %s",
          i2_est, i2_lo, i2_hi, df, Q_meta,
          ifelse(meta$QEp < 0.001, "<0.001", sprintf("%.3f", meta$QEp)))
}
cat("\n", i2_label, "\n")


# =============================================================================
# 5. PREPARE PLOT DATA
# =============================================================================

# Alphabetical, A at top → reverse for ggplot y-axis (bottom to top)
results <- results |>
  arrange(trial) |>
  mutate(
    weight_pct = n / sum(n) * 100,
    sig        = p.value < 0.05,
    p_label    = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    label      = sprintf("%.2f [%.2f, %.2f]   p=%s   W=%.1f%%",
                         est, lo95, hi95, p_label, weight_pct),
    trial      = factor(trial, levels = rev(trial))
  )

K         <- nrow(results)
trial_y   <- seq(K, 1)
overall_y <- -1

plot_trials <- results |> mutate(y = trial_y)

# Dot size proportional to weight (manual, no size legend)
min_pt <- 2; max_pt <- 7
w_min  <- min(results$weight_pct); w_max <- max(results$weight_pct)
plot_trials <- plot_trials |>
  mutate(pt_size = min_pt + (weight_pct - w_min) / (w_max - w_min) * (max_pt - min_pt))

overall_label <- sprintf("%.2f [%.2f, %.2f]   p=%s",
                         overall$est, overall$lo95, overall$hi95,
                         ifelse(overall$p.value < 0.001, "<0.001",
                                sprintf("%.3f", overall$p.value)))

x_lo    <- min(c(results$lo95, overall$lo95)) - 0.1
x_hi    <- max(c(results$hi95, overall$hi95)) + 0.1
label_x <- x_hi + (x_hi - x_lo) * 0.02

y_breaks <- c(trial_y, overall_y)
y_labels <- c(as.character(results$trial), "Overall")


# =============================================================================
# 6. FOREST PLOT
# =============================================================================

forest_plot <- ggplot() +
  
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  
  geom_hline(yintercept = (overall_y + 1) / 2,
             colour = "grey60", linewidth = 0.35) +
  
  # CI lines
  geom_errorbarh(
    data = plot_trials,
    aes(y = y, xmin = lo95, xmax = hi95),
    height = 0.18, linewidth = 0.5, colour = "grey30"
  ) +
  
  # Points (manual size, no legend)
  geom_point(
    data  = plot_trials,
    aes(x = est, y = y, colour = sig, fill = sig),
    size  = plot_trials$pt_size,
    shape = 21, stroke = 0.55
  ) +
  
  # Overall diamond
  geom_polygon(
    data = data.frame(
      x = c(overall$lo95, overall$est, overall$hi95, overall$est),
      y = c(overall_y, overall_y + 0.42, overall_y, overall_y - 0.42)
    ),
    aes(x = x, y = y),
    fill = "#185FA5", colour = "#0D3D6E", linewidth = 0.5
  ) +
  
  # Trial labels (right)
  geom_text(
    data = plot_trials,
    aes(x = label_x, y = y, label = label, colour = sig),
    hjust = 0, size = 2.7, show.legend = FALSE
  ) +
  
  # Overall label (right)
  annotate("text", x = label_x, y = overall_y, label = overall_label,
           hjust = 0, size = 2.7, colour = "#185FA5", fontface = "bold") +
  
  # I² annotation
  annotate("text",
           x = (x_lo + x_hi) / 2, y = overall_y - 0.85,
           label = i2_label,
           size = 2.5, colour = "grey35", hjust = 0.5, fontface = "italic") +
  
  # Column header
  annotate("text", x = label_x, y = K + 0.85,
           label = "\u03b2 [95% CI]   p-value   Weight",
           hjust = 0, size = 2.7, colour = "grey30", fontface = "bold") +
  
  scale_colour_manual(
    values = c("TRUE" = "#A32D2D", "FALSE" = "#4A4A47"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "p \u2265 0.05"),
    name   = NULL
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#D96060", "FALSE" = "#AAAAAA"),
    guide  = "none"
  ) +
  scale_x_continuous(
    limits = c(x_lo, x_hi + (x_hi - x_lo) * 0.52),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = y_labels,
    expand = expansion(add = c(1.4, 1.2))
  ) +
  
  labs(
    title    = "Effect of BEC (log10) on Delta Attack rate - per trial",
    subtitle = paste0(
      "Univariate linear regression, adjusted for ", covariate
    ),
    x = "Adjusted coefficient for BEC (log10)  (95% CI)",
    y = NULL
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    panel.border       = element_rect(colour = "grey70", linewidth = 0.5),
    axis.text.y        = element_text(size = 15),
    axis.text.x        = element_text(size = 15),
    plot.title         = element_text(face = "bold", size = 20),
    plot.subtitle      = element_text(size = 8.5, colour = "grey40"),
    legend.position    = "bottom",
    legend.text        = element_text(size = 12),
    plot.margin        = margin(10, 10, 10, 10)
  )

print(forest_plot)


# =============================================================================
# 7. SAVE
# =============================================================================

plot_height <- max(6, K * 0.55 + 3.5)

ggsave(file.path(getwd(), "Forest_plot_BEC_per_trial_MIestimates.pdf"),
       plot = forest_plot, width = 18, height = plot_height, units = "in")

ggsave(file.path(getwd(), "Forest_plot_BEC_per_trial_MIestimates.png"),
       plot = forest_plot, width = 18, height = plot_height,
       dpi = 300, units = "in")

cat("\nSaved to:", getwd(), "\n")


# =============================================================================
# 8. RETURN
# =============================================================================

list(results = results, overall = overall,
     I2 = i2_est, I2_ci = c(i2_lo, i2_hi), i2_label = i2_label)


#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S7 + Table S3/S4 - Asthma attack - Sensitivity analysis : Negative binomial regression
# Compared to linear regression model for the main model (in the manuscript)
#======================================================================================================================================
# === SHIFT DELTA TO NON-NEGATIVE ===
Data_Oracle_attack_shifted <- ORACLE_Placebo_Attack_imputed %>%
  mutate(
    Delta_attack_shifted = Placebo_response_Attack + abs(min(Placebo_response_Attack, na.rm = TRUE))
  )

# === VERIFY ===
min(Data_Oracle_attack_shifted$Delta_attack_shifted, na.rm = TRUE)  # Should be 0
max(Data_Oracle_attack_shifted$Delta_attack_shifted, na.rm = TRUE)
summary(Data_Oracle_attack_shifted$Delta_attack_shifted)

# === CHECK DISTRIBUTION ===
hist(Data_Oracle_attack_shifted$Delta_attack_shifted,
     main = "Shifted Delta Asthma Attacks",
     xlab = "Delta Attack (shifted)",
     col  = "steelblue")

#======================================================================================================================================
# Figure S6.A - Simple model
#======================================================================================================================================

# === Univariate MODEL (MODEL 1 - Supplementary) - ASTHMA ATTACKS ===
Univariate_deltaattack_model_nb <- MI_estimates(
  data        = Data_Oracle_attack_shifted,
  outcome_var = "Delta_attack_shifted",
  predictor_vars = c(
    "Age_per_10",
    "Sex",
    "BMI_per_5",
    "Smoking_history_yes_no",
    "Airborne_allergen_sensitisation",
    "Allergic_Rhinitis",
    "Eczema",
    "CRSwNP",
    "CRSsNP",
    "Psychiatric_disease",
    "GINA_step_numeric",
    "ACQ_score_0W",
    "Attack_history_num",
    "FEV1_preBD_L_0W_per_10",
    "FEV1_reversibility_0W_per_10",
    "BEC_log10",
    "FeNO_log10",
    "IgE_log10"
  ),
  covariables          = NULL,
  imp_col              = ".imp",
  followup_offset      = "Yes",
  followup_col         = "Follow_up_duration_days",
  random_intercept_var = "Enrolled_Trial_name",
  model_type           = "nb"
)

# === EXTRACT & LABEL RESULTS ===
Univariate_deltaattack_results_nb <- attr(Univariate_deltaattack_model_nb, "combined_results")

row.names(Univariate_deltaattack_results_nb) <- c(
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
  "FEV1 pre-BD (L., per 0.1 decrease)",
  "FEV1 reversibility (per 10% increase)",
  # Inflammatory biomarkers (3)
  "Blood eosinophils (log10)",
  "FeNO (log10)",
  "Total IgE (log10)"
)

# === HANDLE NaN ROW BEFORE PLOTTING ===
results_nb_plot <- Univariate_deltaattack_results_nb

# Replace NaN CI with NA so forplo skips the point cleanly
results_nb_plot["Number of severe exacerbations or hospitalisations (past 12 months)", 
                c("2.5 %", "97.5 %")] <- NA

# Build beta column with header as first element
beta_col <- c(
  "Beta",                                          # header row
  sprintf("%.3f", results_nb_plot$estimate)        # 18 data rows
)
beta_col[is.nan(beta_col) | beta_col == "NaN"] <- ""

pdf(
  "/Users/joel/ORACLE - Placebo - Exacerbations/Project/FigureS7A_Forplo_Attacks_Univar_NB.pdf",
  width  = 22,
  height = 14
)
par(mar = c(4, 2, 2, 2))

SUP_Forplo_Attacks_Univar_NB <- forplo(
  as.data.frame(results_nb_plot[, c("estimate", "2.5 %", "97.5 %")]),
  
  linreg = TRUE,
  xlim   = c(-0.06, 0.06),
  
  em     = "",
  
  row.labels = row.names(results_nb_plot),
  left.align = FALSE,
  
  add.arrow.right = FALSE,
  add.arrow.left  = FALSE,
  left.bar        = FALSE,
  right.bar       = FALSE,
  
  shade.every = 1,
  shade.col   = "grey",
  shade.alpha = 0.15,
  
  margin.left    = 18,
  margin.right   = 8,
  column.spacing = 1,
  
  add.columns  = list(beta_col),
  add.colnames = "",           # ← blank since header is in beta_col[1]
  
  groups    = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 3)),
  grouplabs = c(
    "Demographics",
    "Comorbidities",
    "Asthma history",
    "Baseline lung function",
    "Inflammatory biomarkers"
  ),
  
  char = 20,
  size = 1.5,
  
  col = c(
    rep("darkgreen",  4),
    rep("purple",     6),
    rep("darkred",    3),
    rep("darkblue",   2),
    rep("darkorange", 3)
  )
)
dev.off()

#======================================================================================================================================
# Figure S6.B - Complete model (with covariates)
#======================================================================================================================================
# === Multivariate MODEL (MODEL 2 - Supplementary) - ASTHMA ATTACKS ===
Multivariate_deltaattack_model_nb <- MI_estimates(
  data        = Data_Oracle_attack_shifted,
  outcome_var = "Delta_attack_shifted",
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
    "FEV1_preBD_L_0W_per_10",             # 14 - Baseline lung function
    "FEV1_reversibility_0W_per_10",        # 15 - Baseline lung function
    "BEC_log10",                           # 16 - Inflammatory biomarkers
    "FeNO_log10",                          # 17 - Inflammatory biomarkers
    "IgE_log10"                            # 18 - Inflammatory biomarkers
  ),
  covariables = c(
    "Sex",
    "BMI_per_5",
    "Smoking_history_yes_no",
    "CRSwNP",
    "Psychiatric_disease",
    "GINA_step_numeric",
    "ACQ_score_0W", 
    "Attack_history_num",
    "FEV1_preBD_L_0W_per_10",
    "BEC_log10",
    "FeNO_log10"
  ),
  imp_col              = ".imp",
  followup_offset      = "Yes",
  followup_col         = "Follow_up_duration_days",
  random_intercept_var = "Enrolled_Trial_name",
  model_type           = "nb"
)

# === EXTRACT & LABEL RESULTS ===
Multivariate_deltaattack_results_nb <- attr(Multivariate_deltaattack_model_nb, "combined_results")

row.names(Multivariate_deltaattack_results_nb) <- c(
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
  "FEV1 pre-BD (L., per 0.1 decrease)",
  "FEV1 reversibility (per 10% increase)",
  # Inflammatory biomarkers (3)
  "Blood eosinophils (log10)",
  "FeNO (log10)",
  "Total IgE (log10)"
)

# === HANDLE NaN ROW BEFORE PLOTTING ===
results_multivar_nb_plot <- Multivariate_deltaattack_results_nb
results_multivar_nb_plot["Exacerbation history (count)",
                         c("2.5 %", "97.5 %")] <- 0   # ← 0 instead of NA

ci_col <- sprintf("%.3f [%.3f; %.3f]",
                  results_multivar_nb_plot$estimate,
                  results_multivar_nb_plot$`2.5 %`,
                  results_multivar_nb_plot$`97.5 %`)
ci_col[results_multivar_nb_plot$`2.5 %` == 0] <- ""   # ← blank text for that row

pdf(
  "/Users/joel/ORACLE - Placebo - Exacerbations/Project/FigureS7B_Forplo_Attacks_Multivar_NB.pdf",
  width  = 22,
  height = 14
)
par(mar = c(4, 2, 2, 2))

SUP_Forplo_Attacks_Multivar_NB <- forplo(
  as.data.frame(results_multivar_nb_plot[, c("estimate", "2.5 %", "97.5 %")]),
  
  linreg = TRUE,
  xlim   = c(-0.06, 0.06),
  
  em     = "",               # ← suppress default CI column
  
  row.labels = row.names(results_multivar_nb_plot),
  left.align = FALSE,
  
  add.arrow.right = FALSE,
  add.arrow.left  = FALSE,
  left.bar        = FALSE,
  
  shade.every = 1,
  shade.col   = "grey",
  shade.alpha = 0.15,
  
  margin.left  = 18,
  margin.right = 18,
  
  add.columns  = list(ci_col),
  add.colnames = "Beta [95% CI]",
  
  groups    = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 3)),
  grouplabs = c(
    "Demographic",
    "Comorbidities",
    "Asthma history",
    "Baseline lung function",
    "Inflammatory biomarkers"
  ),
  
  char = 20,
  size = 1.5,
  
  col = c(
    rep("darkgreen",  4),
    rep("purple",     6),
    rep("darkred",    3),
    rep("darkblue",   2),
    rep("darkorange", 3)
  )
)
dev.off()



#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S7. Asthma attack - Sensitivity analysis – CAPTAIN + PACT 
# A) All trials (with CAPTAIN + PACT)
# B) Only placebo trials (without CAPTAIN + PACT) - MAIN ANALYSIS IN THE MANUSCRIPT
# c) Only CAPTAIN + PACT
#======================================================================================================================================

#Sensitivity analysis : Trials excluded if no "true placebo arm"
#======================================================================================================================================
#Preparation of dataset
#======================================================================================================================================
Data_Oracle %>%
  filter(.imp == 0) %>%
  filter(Enrolled_Trial_name %in% c("CAPTAIN", "PACT")) %>%
  group_by(Enrolled_Trial_name) %>%
  summarise(
    N = n(),
    # --- Attack variables ---
    NA_Attack_history_num              = sum(is.na(Attack_history_num)),
    NA_Attack_number_during_followup   = sum(is.na(Attack_number_during_followup)),
    # --- FEV1 ---
    NA_FEV1_preBD_L_0W                 = sum(is.na(FEV1_preBD_L_0W)),
    NA_FEV1_preBD_L_52W                = sum(is.na(FEV1_preBD_L_52W)),
    # --- ACQ ---
    NA_ACQ_score_0W                    = sum(is.na(ACQ_score_0W)),
    NA_ACQ_score_24W                   = sum(is.na(ACQ_score_24W)),
    NA_ACQ_score_26W                   = sum(is.na(ACQ_score_26W)),
    .groups = "drop"
  ) %>%
  print()
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
# Placebo response on asthma attack - Data_Oracle_attack_withPACT_CAPTAIN 
# Selection of the studies with attack history in numeric

# Sensitivity: adding PACT and CAPTAIN
Studies_included_attack_withPACT_CAPTAIN <- c("AZISAST", "BENRAP2B", "DREAM", "DRI12544", "EXTRA",
                                              "LAVOLTA_1", "LAVOLTA_2", "LUSTER_1", "LUSTER_2",
                                              "NAVIGATOR", "PATHWAY", "QUEST", "STRATOS_1", "STRATOS_2", "PACT", "CAPTAIN")

Data_Oracle_attack_withPACT_CAPTAIN <- ORACLE_after_imputation %>%
  filter(Enrolled_Trial_name %in% Studies_included_attack_withPACT_CAPTAIN)

unique(Data_Oracle_attack_withPACT_CAPTAIN$Enrolled_Trial_name)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Preparation for delta asthma attack analysis - Data_Oracle_attack_withPACT_CAPTAIN 
sum(is.na(Data_Oracle_attack_withPACT_CAPTAIN $Attack_number_during_followup))
sum(is.na(Data_Oracle_attack_withPACT_CAPTAIN $Attack_history_num))
sum(is.na(Data_Oracle_attack_withPACT_CAPTAIN $Placebo_response_Attack))

# Delta attack
Data_Oracle_attack_withPACT_CAPTAIN  <- Data_Oracle_attack_withPACT_CAPTAIN  %>%
  mutate(Delta_attack = Attack_history_num - Attack_number_during_followup)

## Calculate Rate Attack in trial
Data_Oracle_attack_withPACT_CAPTAIN  <- Data_Oracle_attack_withPACT_CAPTAIN  %>%
  mutate(
    Attack_rate_before_trial      = Attack_history_num / 1,
    Attack_rate_during_trial      = Attack_number_during_followup / Follow_up_duration_days,
    Delta_attack_rate             = Attack_rate_before_trial - Attack_rate_during_trial,
    Attack_rate_before_trial_log  = log(Attack_rate_before_trial + 1)
  )

#======================================================================================================================================
# UNIVARIATE & MULTIVARIATE ANALYSIS
#======================================================================================================================================

#======================================================================
# === CREATE DATASETS ===
#======================================================================

# A) Original / main analysis — without PACT & CAPTAIN (post-imputation)
Data_Oracle_attack_main <- ORACLE_Placebo_Attack_imputed

# B) Sensitivity — with PACT & CAPTAIN included (post-imputation)
Data_Oracle_attack_withPACT_CAPTAIN <- Data_Oracle_attack_withPACT_CAPTAIN %>%
  # Identify participants where Delta_attack is missing in .imp == 0
  group_by(Sequential_number) %>%
  filter(!any(.imp == 0 & is.na(Placebo_response_Attack))) %>%
  ungroup() %>%
  # Keep only .imp = 1 to 10
  filter(.imp != 0)
unique(Data_Oracle_attack_withPACT_CAPTAIN$.imp)
nrow(Data_Oracle_attack_withPACT_CAPTAIN)

## C) Preparation of dataset : Only PACT and CAPTAIN
Data_Oracle_attack_onlyPACT_CAPTAIN <- Data_Oracle_attack_withPACT_CAPTAIN %>%
  filter(Enrolled_Trial_name %in% c("PACT","CAPTAIN"))
unique(Data_Oracle_attack_onlyPACT_CAPTAIN$Enrolled_Trial_name)
sum(is.na(Data_Oracle_attack_onlyPACT_CAPTAIN$Placebo_response_Attack))

# === SANITY CHECK
cat("\nDataset sizes (imp = 1):\n")
cat("Main (without PACT & CAPTAIN) :", nrow(Data_Oracle_attack_main             %>% filter(.imp == 1)), "\n")
cat("With PACT & CAPTAIN           :", nrow(Data_Oracle_attack_withPACT_CAPTAIN %>% filter(.imp == 1)), "\n")
cat("With PACT & CAPTAIN           :", nrow(Data_Oracle_attack_onlyPACT_CAPTAIN %>% filter(.imp == 1)), "\n")

#======================================================================================================================================
# === UNIVARIATE MODEL - FUNCTION TO AVOID REPETITION ===

run_univariate_attack <- function(data, label) {
  
  cat("\nRunning univariate model:", label, "\n")
  
  # === RUN MODEL ===
  model <- MI_estimates(
    data           = data,
    outcome_var    = "Placebo_response_Attack",
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
      "FEV1_preBD_L_0W_per_10",             # 14 - Baseline lung function
      "FEV1_reversibility_0W_per_10",        # 15 - Baseline lung function
      "BEC_log10",                           # 16 - Inflammatory biomarkers
      "FeNO_log10",                          # 17 - Inflammatory biomarkers
      "IgE_log10"                            # 18 - Inflammatory biomarkers
    ),
    covariables          = c("Attack_history_num"),
    imp_col              = ".imp",
    followup_offset      = "Yes",
    followup_col         = "Follow_up_duration_days",
    random_intercept_var = "Enrolled_Trial_name",
    model_type           = "lm"
  )
  
  # === EXTRACT RESULTS ===
  results <- attr(model, "combined_results")
  
  # === LABEL ROWS ===
  row.names(results) <- c(
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
    "FEV1 pre-BD (L., per 0.1 decrease)",
    "FEV1 reversibility (per 10% increase)",
    # Inflammatory biomarkers (3)
    "Blood eosinophils (log10)",
    "FeNO (log10)",
    "Total IgE (log10)"
  )
  
  return(results)
}

#======================================================================
# === RUN UNIVARIATE MODELS ===
#======================================================================

# 1. Main analysis (without PACT & CAPTAIN)
Univariate_deltaattack_results_main <- run_univariate_attack(
  data  = Data_Oracle_attack_main,
  label = "Main analysis (without PACT & CAPTAIN)"
)

# 2. Sensitivity — with PACT & CAPTAIN
Univariate_deltaattack_results_withPACT_CAPTAIN <- run_univariate_attack(
  data  = Data_Oracle_attack_withPACT_CAPTAIN,
  label = "Sensitivity: with PACT & CAPTAIN"
)

# 3. Sensitivity — only PACT & CAPTAIN
Univariate_deltaattack_results_onlyPACT_CAPTAIN <- run_univariate_attack(
  data  = Data_Oracle_attack_onlyPACT_CAPTAIN,
  label = "Sensitivity: only PACT & CAPTAIN"
)

#======================================================================================================================================
# === FOREST PLOT - FUNCTION TO AVOID REPETITION ===

plot_univariate_attack <- function(results, label, filename) {
  pdf(filename, width = 18, height = 12)
  par(mar = c(4, 2, 2, 2))
  forplo(
    as.data.frame(results[, c("estimate", "2.5 %", "97.5 %")]),
    xlim            = c(-1.5, 1.5),
    em              = "aRC",
    linreg          = TRUE,
    row.labels      = row.names(results),
    left.align      = FALSE,
    add.arrow.right = FALSE,
    arrow.right.length = 20,
    add.arrow.left  = FALSE,
    arrow.left.length  = 20,
    left.bar        = FALSE,
    shade.every     = 1,
    shade.col       = "grey",
    shade.alpha     = 0.2,
    margin.left     = 18,
    margin.right    = 18,   # ↑ increased from 8 to pull values left
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
  cat("✓ Saved:", filename, "\n")
}

#======================================================================
# === FOREST PLOTS — UNIVARIATE ===
#======================================================================

base_path <- "/Users/joel/ORACLE - Placebo - Exacerbations/Project/"

plot_univariate_attack(
  results  = Univariate_deltaattack_results_main,
  label    = "Main analysis (without PACT & CAPTAIN)",
  filename = paste0(base_path, "Forrest_plot_Attacks_Univar_Main.pdf")
)

plot_univariate_attack(
  results  = Univariate_deltaattack_results_withPACT_CAPTAIN,
  label    = "Sensitivity: with PACT & CAPTAIN",
  filename = paste0(base_path, "Forrest_plot_Attacks_Univar_WithPACT_CAPTAIN.pdf")
)


plot_univariate_attack(
  results  = Univariate_deltaattack_results_onlyPACT_CAPTAIN,
  label    = "Sensitivity: only PACT & CAPTAIN",
  filename = paste0(base_path, "Forrest_plot_Attacks_Univar_onlyPACT_CAPTAIN.pdf")
)

cat("\n✓ PACT & CAPTAIN sensitivity univariate models completed.\n")

#======================================================================
# === COMBINED FOREST PLOT — 3 PANELS SIDE BY SIDE (base R / forplo) ===
#======================================================================

# --- 1. Save the 3 individual PDFs as usual ---------------------------

plot_univariate_attack(
  results  = Univariate_deltaattack_results_withPACT_CAPTAIN,
  label    = "A) All trials (Placebo + Control)",
  filename = paste0(base_path, "tmp_panel_A.pdf")
)

plot_univariate_attack(
  results  = Univariate_deltaattack_results_main,
  label    = "B) Only placebo (main analysis)",
  filename = paste0(base_path, "tmp_panel_B.pdf")
)

plot_univariate_attack(
  results  = Univariate_deltaattack_results_onlyPACT_CAPTAIN,
  label    = "C) Only control (PACT + CAPTAIN)",
  filename = paste0(base_path, "tmp_panel_C.pdf")
)

# --- 2. Combine the 3 PDFs side by side in one page -------------------
install.packages("pdftools")
library(pdftools)   # install.packages("pdftools") if needed

# Rasterize each PDF page to bitmap, then stitch horizontally
bitmaps <- lapply(
  c(paste0(base_path, "tmp_panel_A.pdf"),
    paste0(base_path, "tmp_panel_B.pdf"),
    paste0(base_path, "tmp_panel_C.pdf")),
  function(f) pdf_render_page(f, page = 1, dpi = 200)
)

# --- 3. Write combined PNG (side by side) ------------------------------
library(png)   # install.packages("png") if needed

# Check all panels have same height; use max for safety
heights <- sapply(bitmaps, function(b) dim(b)[1])
widths  <- sapply(bitmaps, function(b) dim(b)[2])
max_h   <- max(heights)

# Pad shorter panels to same height (white padding)
bitmaps_padded <- lapply(bitmaps, function(b) {
  h_diff <- max_h - dim(b)[1]
  if (h_diff > 0) {
    pad <- array(1, dim = c(h_diff, dim(b)[2], dim(b)[3]))
    b   <- abind::abind(b, pad, along = 1)
  }
  b
})

combined_bitmap <- abind::abind(bitmaps_padded[[1]],
                                bitmaps_padded[[2]],
                                bitmaps_padded[[3]], along = 2)

writePNG(combined_bitmap,
         paste0(base_path, "Forrest_plot_Attacks_Univar_Combined_3panels.png"))

cat("\n✓ Combined 3-panel PNG saved.\n")

# --- 4. Optional: convert PNG → PDF -----------------------------------
install.packages("magick") 
library(magick)

magick::image_read(paste0(base_path, "Forrest_plot_Attacks_Univar_Combined_3panels.png")) %>%
  magick::image_convert("pdf") %>%
  magick::image_write(paste0(base_path, "Forrest_plot_Attacks_Univar_Combined_3panels.pdf"))

cat("✓ Combined 3-panel PDF saved.\n")

# --- 5. Clean up temp files -------------------------------------------
file.remove(paste0(base_path, c("tmp_panel_A.pdf",
                                "tmp_panel_B.pdf",
                                "tmp_panel_C.pdf")))

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# MULTIVARIATE ANALYSIS
#======================================================================================================================================
#======================================================================================================================================
# === MULTIVARIATE MODEL - FUNCTION TO AVOID REPETITION ===
run_multivariate_attack <- function(data, label) {
  
  cat("\nRunning multivariate model:", label, "\n")
  
  # === RUN MODEL ===
  model <- MI_estimates(
    data           = data,
    outcome_var    = "Placebo_response_Attack",
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
      "FEV1_preBD_L_0W_per_10",             # 14 - Baseline lung function
      "FEV1_reversibility_0W_per_10",        # 15 - Baseline lung function
      "BEC_log10",                           # 16 - Inflammatory biomarkers
      "FeNO_log10",                          # 17 - Inflammatory biomarkers
      "IgE_log10"                            # 18 - Inflammatory biomarkers
    ),
    covariables = c(
      "Sex",
      "BMI_per_5",
      "Smoking_history_yes_no",
      "CRSwNP",
      "Psychiatric_disease",
      "GINA_step_numeric",
      "ACQ_score_0W",
      "Attack_history_num",
      "FEV1_preBD_L_0W_per_10",
      "BEC_log10",
      "FeNO_log10"
    ),
    imp_col              = ".imp",
    followup_offset      = "Yes",
    followup_col         = "Follow_up_duration_days",
    random_intercept_var = "Enrolled_Trial_name",
    model_type           = "lm"
  )
  
  # === EXTRACT RESULTS ===
  results <- attr(model, "combined_results")
  
  # === LABEL ROWS ===
  row.names(results) <- c(
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
    "FEV1 pre-BD (L., per 0.1 decrease)",
    "FEV1 reversibility (per 10% increase)",
    # Inflammatory biomarkers (3)
    "Blood eosinophils (log10)",
    "FeNO (log10)",
    "Total IgE (log10)"
  )
  
  return(results)
}

#======================================================================================================================================
# === FOREST PLOT - FUNCTION (MULTIVARIATE) ===

plot_multivariate_attack <- function(results, label, filename) {
  
  pdf(filename, width = 18, height = 12)
  par(mar = c(4, 2, 2, 2))
  
  forplo(
    as.data.frame(results[, c("estimate", "2.5 %", "97.5 %")]),
    xlim            = c(-1.5, 1.5),
    em              = "aRC",
    linreg          = TRUE,
    row.labels      = row.names(results),
    left.align      = FALSE,
    add.arrow.right = FALSE,
    arrow.right.length = 20,
    add.arrow.left  = FALSE,
    arrow.left.length  = 20,
    left.bar        = FALSE,
    shade.every     = 1,
    shade.col       = "grey",
    shade.alpha     = 0.2,
    margin.left     = 18,
    margin.right    = 18,
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
  cat("✓ Saved:", filename, "\n")
}

#======================================================================
# === RUN MULTIVARIATE MODELS ===
#======================================================================

# 1. Main analysis (without PACT & CAPTAIN)
Multivariate_deltaattack_results_main <- run_multivariate_attack(
  data  = Data_Oracle_attack_main,
  label = "Main analysis (without PACT & CAPTAIN)"
)

# 2. Sensitivity — with PACT & CAPTAIN
Multivariate_deltaattack_results_withPACT_CAPTAIN <- run_multivariate_attack(
  data  = Data_Oracle_attack_withPACT_CAPTAIN,
  label = "Sensitivity: with PACT & CAPTAIN"
)

# 3. Sensitivity — only PACT & CAPTAIN
Multivariate_deltaattack_results_onlyPACT_CAPTAIN <- run_multivariate_attack(
  data  = Data_Oracle_attack_onlyPACT_CAPTAIN,
  label = "Sensitivity: only PACT & CAPTAIN"
)

#======================================================================
# === FOREST PLOTS — MULTIVARIATE ===
#======================================================================

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_main,
  label    = "Main analysis (without PACT & CAPTAIN)",
  filename = paste0(base_path, "Forrest_plot_Attacks_Multivar_Main.pdf")
)

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_withPACT_CAPTAIN,
  label    = "Sensitivity: with PACT & CAPTAIN",
  filename = paste0(base_path, "Forrest_plot_Attacks_Multivar_WithPACT_CAPTAIN.pdf")
)

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_onlyPACT_CAPTAIN,
  label    = "Sensitivity: only PACT & CAPTAIN",
  filename = paste0(base_path, "Forrest_plot_Attacks_Multivar_onlyPACT_CAPTAIN.pdf")
)

cat("\n✓ PACT & CAPTAIN sensitivity multivariate models completed.\n")

#=================================================================================================
# === stitch 3 PDFs side by side → one combined PDF === for Univariable and Multivariable
#=================================================================================================
library(pdftools)  
library(png)
library(magick)

stitch_3panels <- function(file_A, file_B, file_C,
                           title_A, title_B, title_C,
                           filename_out,
                           dpi = 200) {
  
  library(magick)
  
  # --- Read & annotate each panel with its title ---------------------
  add_title <- function(pdf_path, title_text) {
    img <- magick::image_read_pdf(pdf_path, density = dpi)
    magick::image_annotate(
      img,
      text     = title_text,
      gravity  = "North",
      location = "+0+20",
      size     = 40,
      font     = "Arial",
      color    = "black",
      weight   = 700
    )
  }
  
  panel_A <- add_title(file_A, title_A)
  panel_B <- add_title(file_B, title_B)
  panel_C <- add_title(file_C, title_C)
  
  # --- Normalise heights (pad shorter panels with white) -------------
  h <- max(magick::image_info(panel_A)$height,
           magick::image_info(panel_B)$height,
           magick::image_info(panel_C)$height)
  
  pad <- function(img) {
    info <- magick::image_info(img)
    if (info$height < h) {
      magick::image_extent(img,
                           geometry = paste0(info$width, "x", h),
                           gravity  = "North",
                           color    = "white")
    } else { img }
  }
  
  panel_A <- pad(panel_A)
  panel_B <- pad(panel_B)
  panel_C <- pad(panel_C)
  
  # --- Stitch horizontally -------------------------------------------
  combined <- magick::image_append(c(panel_A, panel_B, panel_C),
                                   stack = FALSE)   # FALSE = side by side
  
  # --- Save as PDF -----------------------------------------------------
  magick::image_write(combined, path = filename_out, format = "pdf")
  cat("✓ Combined 3-panel PDF saved:", filename_out, "\n")
}

#======================================================================
# === UNIVARIATE — save 3 individual PDFs then stitch ================
#======================================================================

plot_univariate_attack(
  results  = Univariate_deltaattack_results_withPACT_CAPTAIN,
  label    = "A) All trials (Placebo + Control)",
  filename = paste0(base_path, "tmp_Univar_A.pdf")
)
plot_univariate_attack(
  results  = Univariate_deltaattack_results_main,
  label    = "B) Only placebo (main analysis)",
  filename = paste0(base_path, "tmp_Univar_B.pdf")
)
plot_univariate_attack(
  results  = Univariate_deltaattack_results_onlyPACT_CAPTAIN,
  label    = "C) Only control (PACT + CAPTAIN)",
  filename = paste0(base_path, "tmp_Univar_C.pdf")
)

stitch_3panels(
  file_A       = paste0(base_path, "tmp_Univar_A.pdf"),
  file_B       = paste0(base_path, "tmp_Univar_B.pdf"),
  file_C       = paste0(base_path, "tmp_Univar_C.pdf"),
  title_A      = "A) All trials (Placebo + Control)",
  title_B      = "B) Only placebo (main analysis)",
  title_C      = "C) Only control (PACT + CAPTAIN)",
  filename_out = paste0(base_path, "Forrest_plot_Attacks_Univar_Combined_3panels.pdf")
)

# Clean up temp files
file.remove(paste0(base_path, c("tmp_Univar_A.pdf",
                                "tmp_Univar_B.pdf",
                                "tmp_Univar_C.pdf")))

cat("\n✓ Univariate combined 3-panel plot completed.\n")

#======================================================================
# === MULTIVARIATE — save 3 individual PDFs then stitch ==============
#======================================================================

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_withPACT_CAPTAIN,
  label    = "A) All trials (Placebo + Control)",
  filename = paste0(base_path, "tmp_Multivar_A.pdf")
)
plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_main,
  label    = "B) Only placebo (main analysis)",
  filename = paste0(base_path, "tmp_Multivar_B.pdf")
)
plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_onlyPACT_CAPTAIN,
  label    = "C) Only control (PACT + CAPTAIN)",
  filename = paste0(base_path, "tmp_Multivar_C.pdf")
)

stitch_3panels(
  file_A       = paste0(base_path, "tmp_Multivar_A.pdf"),
  file_B       = paste0(base_path, "tmp_Multivar_B.pdf"),
  file_C       = paste0(base_path, "tmp_Multivar_C.pdf"),
  title_A      = "A) All trials (Placebo + Control)",
  title_B      = "B) Only placebo (main analysis)",
  title_C      = "C) Only control (PACT + CAPTAIN)",
  filename_out = paste0(base_path, "Forrest_plot_Attacks_Multivar_Combined_3panels.pdf")
)

# Clean up temp files
file.remove(paste0(base_path, c("tmp_Multivar_A.pdf",
                                "tmp_Multivar_B.pdf",
                                "tmp_Multivar_C.pdf")))

cat("\n✓ Multivariate combined 3-panel plot completed.\n")


#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S9. Asthma attacks - Sensitivity analysis : Attacks inclusion criteria
# >2 ; >1 ; no criterion
#======================================================================================================================================
#Sensitvity analysis : Inclusion criteria - Asthma attacks
#======================================================================================================================================
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

# === JOIN CRITERIA TO IMPUTED DATA ===
Data_Oracle_attack_imputated <- ORACLE_Placebo_Attack_imputed %>%
  dplyr::left_join(trial_criteria_map_attack, by = "Enrolled_Trial_name")

# === CREATE SUBSETS BY EXPLICIT TRIAL NAMES ===

# 1. With ≥2 asthma attacks
Data_Oracle_attack_2plus <- Data_Oracle_attack_imputated %>%
  dplyr::filter(Enrolled_Trial_name %in% c(
    "STRATOS_1", "STRATOS_2",
    "PATHWAY",
    "NAVIGATOR",
    "LUSTER_1",  "LUSTER_2",
    "DREAM",
    "BENRAP2B",
    "AZISAST"
  ))

# 2. With ≥1 asthma attack
Data_Oracle_attack_1plus <- Data_Oracle_attack_imputated %>%
  dplyr::filter(Enrolled_Trial_name %in% c(
    "QUEST",
    "EXTRA",
    "DRI12544"
  ))

# 3. No inclusion criteria
Data_Oracle_attack_nocriteria <- Data_Oracle_attack_imputated %>%
  dplyr::filter(Enrolled_Trial_name %in% c(
    "LAVOLTA_1",
    "LAVOLTA_2"
  ))

# 4. With any inclusion criteria (≥1 + ≥2 combined)
Data_Oracle_attack_withcriteria <- Data_Oracle_attack_imputated %>%
  dplyr::filter(Enrolled_Trial_name %in% c(
    # ≥1 asthma attack
    "QUEST",
    "EXTRA",
    "DRI12544",
    # ≥2 asthma attacks
    "STRATOS_1", "STRATOS_2",
    "PATHWAY",
    "NAVIGATOR",
    "LUSTER_1",  "LUSTER_2",
    "DREAM",
    "BENRAP2B",
    "AZISAST"
  ))


#======================================================================================================================================
# === RUN ALL 5 MODELS ===

# 1. All trials (primary analysis — already run, keeping for consistency)
Univariate_deltaattack_results <- run_univariate_attack(
  data  = Data_Oracle_attack_imputated,
  label = "All trials"
)

# 2. With ≥2 asthma attacks
Univariate_deltaattack_results_2plus <- run_univariate_attack(
  data  = Data_Oracle_attack_2plus,
  label = "With >=2 asthma attacks"
)

# 3. With ≥1 asthma attack
Univariate_deltaattack_results_1plus <- run_univariate_attack(
  data  = Data_Oracle_attack_1plus,
  label = "With >=1 asthma attack"
)

# 4. No inclusion criteria
Univariate_deltaattack_results_nocriteria <- run_univariate_attack(
  data  = Data_Oracle_attack_nocriteria,
  label = "No inclusion criteria"
)

# 5. Any inclusion criteria (≥1 + ≥2 combined)
Univariate_deltaattack_results_withcriteria <- run_univariate_attack(
  data  = Data_Oracle_attack_withcriteria,
  label = "Any inclusion criteria"
)

#======================================================================================================================================
# === GENERATE ALL 5 FOREST PLOTS ===

base_path <- "/Users/joel/ORACLE - Placebo - Exacerbations/Project/"

plot_univariate_attack(
  results  = Univariate_deltaattack_results,
  label    = "All trials",
  filename = paste0(base_path, "Forrest_plot_Attacks_Univar_AllTrials.pdf")
)

plot_univariate_attack(
  results  = Univariate_deltaattack_results_2plus,
  label    = "With >=2 asthma attacks",
  filename = paste0(base_path, "Forrest_plot_Attacks_Univar_2plus.pdf")
)

plot_univariate_attack(
  results  = Univariate_deltaattack_results_1plus,
  label    = "With >=1 asthma attack",
  filename = paste0(base_path, "Forrest_plot_Attacks_Univar_1plus.pdf")
)

plot_univariate_attack(
  results  = Univariate_deltaattack_results_nocriteria,
  label    = "No inclusion criteria",
  filename = paste0(base_path, "Forrest_plot_Attacks_Univar_NoCriteria.pdf")
)

plot_univariate_attack(
  results  = Univariate_deltaattack_results_withcriteria,
  label    = "Any inclusion criteria",
  filename = paste0(base_path, "Forrest_plot_Attacks_Univar_WithCriteria.pdf")
)

cat("\n✓ All 5 models and forest plots completed.\n")


#======================================================================================================================================
# === SENSITIVITY 1 — All trials vs. With any criteria ===

Sensitivity_1 <- compare_results(
  results_list = list(
    Univariate_deltaattack_results,
    Univariate_deltaattack_results_withcriteria
  ),
  names_list = c(
    "All trials (N=17)",
    "With any criteria — ≥1 or ≥2 (N=12)"
  )
)

#======================================================================================================================================
# === SENSITIVITY 2 — ≥1 vs. ≥2 ===

Sensitivity_2 <- compare_results(
  results_list = list(
    Univariate_deltaattack_results_1plus,
    Univariate_deltaattack_results_2plus
  ),
  names_list = c(
    "With ≥1 attack (N=3)",
    "With ≥2 attacks (N=9)"
  )
)

#======================================================================================================================================
# === SENSITIVITY 3 — Gradient: No criteria vs. ≥1 vs. ≥2 ===

Sensitivity_3 <- compare_results(
  results_list = list(
    Univariate_deltaattack_results_nocriteria,
    Univariate_deltaattack_results_1plus,
    Univariate_deltaattack_results_2plus
  ),
  names_list = c(
    "No criteria (N=2)",
    "With ≥1 attack (N=3)",
    "With ≥2 attacks (N=9)"
  )
)

#======================================================================================================================================
# === SENSITIVITY 4 — With any criteria vs. ≥2 only ===

Sensitivity_4 <- compare_results(
  results_list = list(
    Univariate_deltaattack_results_withcriteria,
    Univariate_deltaattack_results_2plus
  ),
  names_list = c(
    "With any criteria — ≥1 or ≥2 (N=12)",
    "With ≥2 attacks only (N=9)"
  )
)

#======================================================================================================================================
# === SENSITIVITY 5 — All trials vs. ≥2 only ===

Sensitivity_5 <- compare_results(
  results_list = list(
    Univariate_deltaattack_results,
    Univariate_deltaattack_results_2plus
  ),
  names_list = c(
    "All trials (N=17)",
    "With ≥2 attacks only (N=9)"
  )
)

#======================================================================================================================================
# === ADD CATEGORY COLUMN TO ALL ===

add_category <- function(df) {
  df %>%
    mutate(
      Category = case_when(
        Variable %in% c(
          "Age (per 10 years)", "Sex (male)",
          "Body Mass Index (per 5 kg/m\u00b2)", "Smoking status")        ~ "Demographics",
        Variable %in% c(
          "Allergy testing positive", "Allergic rhinitis", "Eczema",
          "Chronic Rhinosinusitis without Nasal Polyps",
          "Chronic Rhinosinusitis with Nasal Polyps")                     ~ "Comorbidities",
        Variable %in% c(
          "Treatment step (per step increase)", "ACQ-5",
          "Number of severe exacerbation or hospitalisation in the past 12 months") ~ "Asthma History",
        Variable %in% c(
          "Baseline FEV1 pre-BD (per 10% predicted)",
          "Baseline FEV1 reversibility to BD (per 10%)")                  ~ "Baseline Lung Function",
        Variable %in% c(
          "Blood eosinophils (10\u2079 cells/L) (log)",
          "FeNO (ppb) (log10)",
          "Total IgE (ng/mL) (log10)")                                    ~ "Inflammatory Biomarkers"
      )
    ) %>%
    dplyr::select(Category, Variable, everything())  # <-- explicit dplyr::
}

# Re-run all
Sensitivity_1 <- add_category(Sensitivity_1)
Sensitivity_2 <- add_category(Sensitivity_2)
Sensitivity_3 <- add_category(Sensitivity_3)
Sensitivity_4 <- add_category(Sensitivity_4)
Sensitivity_5 <- add_category(Sensitivity_5)

#======================================================================================================================================
# === FLEXTABLE FUNCTION ===

make_sensitivity_ft <- function(df, title_label) {
  
  ft <- flextable(df) %>%
    theme_vanilla() %>%
    merge_v(j = "Category") %>%
    bold(part = "header") %>%
    bold(j = "Category") %>%
    align(j = "Category", align = "center", part = "body") %>%
    align(j = seq(3, ncol(df)), align = "center", part = "all") %>%
    fontsize(size = 9,  part = "all") %>%
    fontsize(size = 10, part = "header") %>%
    width(j = "Category", width = 1.3) %>%
    width(j = "Variable",  width = 2.8) %>%
    add_footer_lines(
      paste0(
        title_label,
        ". * = 95% CI excludes 0 (significant association).
         Beta [95% CI] shown for each subset.
         Linear mixed-effects model (LMM) with offset for follow-up duration,
         random intercept for trial, adjusted for attack history.
         Pooled across multiple imputations using Rubin's rules."
      )
    ) %>%
    fontsize(size = 8, part = "footer") %>%
    color(color = "grey40", part = "footer") %>%
    autofit()
  
  return(ft)
}

#======================================================================================================================================
# === SAVE ALL TO ONE WORD DOCUMENT ===

base_path <- "/Users/joel/ORACLE - Placebo - Exacerbations/Project/"

doc <- read_docx() %>%
  
  # Title page
  body_add_par(
    "Sensitivity Analyses — Inclusion Criteria for Asthma Attacks",
    style = "heading 1"
  ) %>%
  body_add_par(
    "All comparisons use univariate linear mixed-effects models (LMM),
     adjusted for attack history and follow-up duration,
     with random intercept for trial, pooled across imputations.",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  
  # Sensitivity 1
  body_add_par(
    "Sensitivity 1 — All trials vs. With any inclusion criterion (≥1 or ≥2)",
    style = "heading 2"
  ) %>%
  body_add_par(
    "Does restricting to trials requiring any attack history change results?",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(make_sensitivity_ft(
    Sensitivity_1,
    "Sensitivity 1: All trials (N=17) vs. With any criteria (N=12)"
  )) %>%
  body_add_par("", style = "Normal") %>%
  
  # Sensitivity 2
  body_add_par(
    "Sensitivity 2 — With ≥1 attack vs. With ≥2 attacks",
    style = "heading 2"
  ) %>%
  body_add_par(
    "Does the stringency of the attack history criterion affect results?",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(make_sensitivity_ft(
    Sensitivity_2,
    "Sensitivity 2: With ≥1 attack (N=3) vs. With ≥2 attacks (N=9)"
  )) %>%
  body_add_par("", style = "Normal") %>%
  
  # Sensitivity 3
  body_add_par(
    "Sensitivity 3 — Gradient: No criteria vs. ≥1 vs. ≥2",
    style = "heading 2"
  ) %>%
  body_add_par(
    "Is there a dose-response pattern across enrollment stringency?",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(make_sensitivity_ft(
    Sensitivity_3,
    "Sensitivity 3: No criteria (N=2) vs. ≥1 attack (N=3) vs. ≥2 attacks (N=9)"
  )) %>%
  body_add_par("", style = "Normal") %>%
  
  # Sensitivity 4
  body_add_par(
    "Sensitivity 4 — With any criteria vs. With ≥2 attacks only",
    style = "heading 2"
  ) %>%
  body_add_par(
    "Are ≥1 trials diluting the signal seen in ≥2 trials?",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(make_sensitivity_ft(
    Sensitivity_4,
    "Sensitivity 4: With any criteria (N=12) vs. ≥2 attacks only (N=9)"
  )) %>%
  body_add_par("", style = "Normal") %>%
  
  # Sensitivity 5
  body_add_par(
    "Sensitivity 5 — All trials vs. With ≥2 attacks only",
    style = "heading 2"
  ) %>%
  body_add_par(
    "Does the primary result hold in the most homogeneous population?",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(make_sensitivity_ft(
    Sensitivity_5,
    "Sensitivity 5: All trials (N=17) vs. ≥2 attacks only (N=9)"
  )) %>%
  body_add_par("", style = "Normal") %>%
  
  # Abbreviations
  body_add_par("Abbreviations", style = "heading 2") %>%
  body_add_par(
    "ACQ = Asthma Control Questionnaire; BD = bronchodilator; BMI = body mass
     index; BEC = blood eosinophil count; CI = confidence interval; CRSsNP =
     chronic rhinosinusitis without nasal polyps; CRSwNP = chronic rhinosinusitis
     with nasal polyps; FEV1 = forced expiratory volume in 1 second; FeNO =
     fractional exhaled nitric oxide; IgE = immunoglobulin E; LMM = linear
     mixed-effects model; MI = multiple imputation.",
    style = "Normal"
  )

print(doc,
      target = paste0(base_path, "Sensitivity_Analyses_Criteria_Attacks.docx"))

cat("✓ All 5 sensitivity analyses saved to one Word document.\n")

#======================================================================================================================================
# === RUN ALL 5 MULTIVARIATE MODELS ===

base_path <- "/Users/joel/ORACLE - Placebo - Exacerbations/Project/"

# 1. All trials (primary)
Multivariate_deltaattack_results <- run_multivariate_attack(
  data  = Data_Oracle_attack_imputated,
  label = "All trials"
)

# 2. With ≥2 asthma attacks
Multivariate_deltaattack_results_2plus <- run_multivariate_attack(
  data  = Data_Oracle_attack_2plus,
  label = "With >=2 asthma attacks"
)

# 3. With ≥1 asthma attack
Multivariate_deltaattack_results_1plus <- run_multivariate_attack(
  data  = Data_Oracle_attack_1plus,
  label = "With >=1 asthma attack"
)

# 4. No inclusion criteria
Multivariate_deltaattack_results_nocriteria <- run_multivariate_attack(
  data  = Data_Oracle_attack_nocriteria,
  label = "No inclusion criteria"
)

# 5. Any inclusion criteria (≥1 + ≥2 combined)
Multivariate_deltaattack_results_withcriteria <- run_multivariate_attack(
  data  = Data_Oracle_attack_withcriteria,
  label = "Any inclusion criteria"
)

#======================================================================================================================================
# === GENERATE ALL 5 FOREST PLOTS ===

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results,
  label    = "All trials",
  filename = paste0(base_path, "Forrest_plot_Attacks_Multivar_AllTrials.pdf")
)

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_2plus,
  label    = "With >=2 asthma attacks",
  filename = paste0(base_path, "Forrest_plot_Attacks_Multivar_2plus.pdf")
)

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_1plus,
  label    = "With >=1 asthma attack",
  filename = paste0(base_path, "Forrest_plot_Attacks_Multivar_1plus.pdf")
)

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_nocriteria,
  label    = "No inclusion criteria",
  filename = paste0(base_path, "Forrest_plot_Attacks_Multivar_NoCriteria.pdf")
)

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_withcriteria,
  label    = "Any inclusion criteria",
  filename = paste0(base_path, "Forrest_plot_Attacks_Multivar_WithCriteria.pdf")
)

cat("\n✓ All 5 multivariate models and forest plots completed.\n")

#======================================================================================================================================
# === SENSITIVITY COMPARISONS (MULTIVARIATE) ===

# Sensitivity 1 — All trials vs. With any criteria
Sensitivity_Multi_1 <- compare_results(
  results_list = list(
    Multivariate_deltaattack_results,
    Multivariate_deltaattack_results_withcriteria
  ),
  names_list = c(
    "All trials (N=17)",
    "With any criteria — ≥1 or ≥2 (N=12)"
  )
)

# Sensitivity 2 — ≥1 vs. ≥2
Sensitivity_Multi_2 <- compare_results(
  results_list = list(
    Multivariate_deltaattack_results_1plus,
    Multivariate_deltaattack_results_2plus
  ),
  names_list = c(
    "With ≥1 attack (N=3)",
    "With ≥2 attacks (N=9)"
  )
)

# Sensitivity 3 — Gradient: No criteria vs. ≥1 vs. ≥2
Sensitivity_Multi_3 <- compare_results(
  results_list = list(
    Multivariate_deltaattack_results_nocriteria,
    Multivariate_deltaattack_results_1plus,
    Multivariate_deltaattack_results_2plus
  ),
  names_list = c(
    "No criteria (N=2)",
    "With ≥1 attack (N=3)",
    "With ≥2 attacks (N=9)"
  )
)

# Sensitivity 4 — With any criteria vs. ≥2 only
Sensitivity_Multi_4 <- compare_results(
  results_list = list(
    Multivariate_deltaattack_results_withcriteria,
    Multivariate_deltaattack_results_2plus
  ),
  names_list = c(
    "With any criteria — ≥1 or ≥2 (N=12)",
    "With ≥2 attacks only (N=9)"
  )
)

# Sensitivity 5 — All trials vs. ≥2 only
Sensitivity_Multi_5 <- compare_results(
  results_list = list(
    Multivariate_deltaattack_results,
    Multivariate_deltaattack_results_2plus
  ),
  names_list = c(
    "All trials (N=17)",
    "With ≥2 attacks only (N=9)"
  )
)

#======================================================================================================================================
# === ADD CATEGORIES ===

Sensitivity_Multi_1 <- add_category(Sensitivity_Multi_1)
Sensitivity_Multi_2 <- add_category(Sensitivity_Multi_2)
Sensitivity_Multi_3 <- add_category(Sensitivity_Multi_3)
Sensitivity_Multi_4 <- add_category(Sensitivity_Multi_4)
Sensitivity_Multi_5 <- add_category(Sensitivity_Multi_5)

#======================================================================================================================================
# === FLEXTABLE FUNCTION (MULTIVARIATE) ===

make_sensitivity_multi_ft <- function(df, title_label) {
  
  ft <- flextable(df) %>%
    theme_vanilla() %>%
    merge_v(j = "Category") %>%
    bold(part = "header") %>%
    bold(j = "Category") %>%
    align(j = "Category", align = "center", part = "body") %>%
    align(j = seq(3, ncol(df)), align = "center", part = "all") %>%
    fontsize(size = 9,  part = "all") %>%
    fontsize(size = 10, part = "header") %>%
    width(j = "Category", width = 1.3) %>%
    width(j = "Variable",  width = 2.8) %>%
    add_footer_lines(
      paste0(
        title_label,
        ". * = 95% CI excludes 0 (significant association).
         Beta [95% CI] shown for each subset.
         Multivariable linear mixed-effects model (LMM) with offset for
         follow-up duration, random intercept for trial.
         Covariates: attack history, sex, smoking, allergic rhinitis, CRSsNP,
         treatment step, ACQ-5, baseline FEV1, BEC, FEV1 reversibility, FeNO.
         Pooled across multiple imputations using Rubin's rules."
      )
    ) %>%
    fontsize(size = 8, part = "footer") %>%
    color(color = "grey40", part = "footer") %>%
    autofit()
  
  return(ft)
}

#======================================================================================================================================
# === SAVE ALL TO ONE WORD DOCUMENT ===

doc <- read_docx() %>%
  
  body_add_par(
    "Sensitivity Analyses — Inclusion Criteria for Asthma Attacks (Multivariable)",
    style = "heading 1"
  ) %>%
  body_add_par(
    "All comparisons use multivariable linear mixed-effects models (LMM),
     with offset for follow-up duration and random intercept for trial.
     Covariates: attack history, sex, smoking status, allergic rhinitis,
     CRSsNP, treatment step, ACQ-5, baseline FEV1, BEC,
     FEV1 reversibility, and FeNO. Pooled across imputations.",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  
  # Sensitivity 1
  body_add_par(
    "Sensitivity 1 — All trials vs. With any inclusion criterion (≥1 or ≥2)",
    style = "heading 2"
  ) %>%
  body_add_par(
    "Does restricting to trials requiring any attack history change results?",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(make_sensitivity_multi_ft(
    Sensitivity_Multi_1,
    "Sensitivity 1: All trials (N=17) vs. With any criteria (N=12)"
  )) %>%
  body_add_par("", style = "Normal") %>%
  
  # Sensitivity 2
  body_add_par(
    "Sensitivity 2 — With ≥1 attack vs. With ≥2 attacks",
    style = "heading 2"
  ) %>%
  body_add_par(
    "Does the stringency of the attack history criterion affect results?",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(make_sensitivity_multi_ft(
    Sensitivity_Multi_2,
    "Sensitivity 2: With ≥1 attack (N=3) vs. With ≥2 attacks (N=9)"
  )) %>%
  body_add_par("", style = "Normal") %>%
  
  # Sensitivity 3
  body_add_par(
    "Sensitivity 3 — Gradient: No criteria vs. ≥1 vs. ≥2",
    style = "heading 2"
  ) %>%
  body_add_par(
    "Is there a dose-response pattern across enrollment stringency?",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(make_sensitivity_multi_ft(
    Sensitivity_Multi_3,
    "Sensitivity 3: No criteria (N=2) vs. ≥1 attack (N=3) vs. ≥2 attacks (N=9)"
  )) %>%
  body_add_par("", style = "Normal") %>%
  
  # Sensitivity 4
  body_add_par(
    "Sensitivity 4 — With any criteria vs. With ≥2 attacks only",
    style = "heading 2"
  ) %>%
  body_add_par(
    "Are ≥1 trials diluting the signal seen in ≥2 trials?",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(make_sensitivity_multi_ft(
    Sensitivity_Multi_4,
    "Sensitivity 4: With any criteria (N=12) vs. ≥2 attacks only (N=9)"
  )) %>%
  body_add_par("", style = "Normal") %>%
  
  # Sensitivity 5
  body_add_par(
    "Sensitivity 5 — All trials vs. With ≥2 attacks only",
    style = "heading 2"
  ) %>%
  body_add_par(
    "Does the primary result hold in the most homogeneous population?",
    style = "Normal"
  ) %>%
  body_add_par("", style = "Normal") %>%
  body_add_flextable(make_sensitivity_multi_ft(
    Sensitivity_Multi_5,
    "Sensitivity 5: All trials (N=17) vs. ≥2 attacks only (N=9)"
  )) %>%
  body_add_par("", style = "Normal") %>%
  
  # Abbreviations
  body_add_par("Abbreviations", style = "heading 2") %>%
  body_add_par(
    "ACQ = Asthma Control Questionnaire; BD = bronchodilator; BEC = blood
     eosinophil count; BMI = body mass index; CI = confidence interval;
     CRSsNP = chronic rhinosinusitis without nasal polyps; CRSwNP = chronic
     rhinosinusitis with nasal polyps; FEV1 = forced expiratory volume in 1
     second; FeNO = fractional exhaled nitric oxide; IgE = immunoglobulin E;
     LMM = linear mixed-effects model; MI = multiple imputation.",
    style = "Normal"
  )

print(doc,
      target = paste0(base_path,
                      "Sensitivity_Analyses_Criteria_Attacks_Multivariate.docx"))

cat("✓ All 5 multivariate sensitivity analyses saved.\n")


#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S9. Asthma attacks - Distribution (histogram) of follow-up duration by participants
#======================================================================================================================================
subset_1 <- Data_Oracle_attack %>%
  filter(.imp == 0, !is.na(Delta_attack))

##Distribution of follow-up duration

pdf("Figure_S9_Histogram_FollowupDuration_Attacks_noshort.pdf",
    width = 10, height = 8)
ggplot(subset_1, aes(x = Follow_up_duration_days)) +
  geom_histogram(
    binwidth = 1,
    fill = "steelblue",
    color = "black",
    alpha = 0.7
  ) +
  scale_x_continuous(breaks = seq(0, 370, 50)) +
  labs(
    title = "Distribution of Follow-up duration",
    x = "Days",
    y = "Count"
  )
dev.off()


#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S10. Asthma attacks - Sensitivity analysis for short follow-up duration 
# with (main analysis) Vs. without <26 weeks’ trials (AZISTAT and DRI12544)
#======================================================================================================================================

#======================================================================
# === CREATE DATASETS: ALL TRIALS vs WITHOUT SHORT FOLLOW-UP TRIALS ===
#======================================================================

# A) All trials (primary analysis)
Data_Oracle_attack_all <- ORACLE_Placebo_Attack_imputed

# B) Excluding short follow-up trials (AZISTAT + DRI12544)
Data_Oracle_attack_noshort <- ORACLE_Placebo_Attack_imputed %>%
  filter(!Enrolled_Trial_name %in% c("AZISAST", "DRI12544"))
unique(Data_Oracle_attack_noshort$Enrolled_Trial_name)
# === SANITY CHECK
cat("\nDataset sizes (imp = 1):\n")
cat("All trials                :", nrow(Data_Oracle_attack_all     %>% filter(.imp == 1)), "\n")
cat("Exclude short follow-up   :", nrow(Data_Oracle_attack_noshort %>% filter(.imp == 1)), "\n")

#======================================================================================================================================
# UNIVARIATE ANALYSIS
#======================================================================================================================================

#======================================================================
# === RUN UNIVARIATE MODELS ===
#======================================================================

# 1. All trials
Univariate_deltaattack_results_all <- run_univariate_attack(
  data  = Data_Oracle_attack_all,
  label = "All trials"
)

# 2. Excluding short follow-up trials (AZISTAT + DRI12544)
Univariate_deltaattack_results_noshort <- run_univariate_attack(
  data  = Data_Oracle_attack_noshort,
  label = "Exclude short follow-up trials (AZISTAT + DRI12544)"
)

#======================================================================
# === FOREST PLOTS — SHORT FOLLOW-UP SENSITIVITY ===
#======================================================================

base_path <- "/Users/joel/ORACLE - Placebo - Exacerbations/Project/"

plot_univariate_attack(
  results  = Univariate_deltaattack_results_all,
  label    = "All trials",
  filename = paste0(base_path,
                    "Forrest_plot_Attacks_Univar_AllTrials.pdf")
)

plot_univariate_attack(
  results  = Univariate_deltaattack_results_noshort,
  label    = "Exclude short follow-up trials (AZISTAT + DRI12544)",
  filename = paste0(base_path,
                    "Forrest_plot_Attacks_Univar_NoShortFollowup.pdf")
)

cat("\n✓ Short follow-up sensitivity univariate models completed.\n")

#======================================================================================================================================
# MULTIVARIATE ANALYSIS
#======================================================================================================================================

#======================================================================
# === RUN MULTIVARIATE MODELS ===
#======================================================================

# 1. All trials
Multivariate_deltaattack_results_all <- run_multivariate_attack(
  data  = Data_Oracle_attack_all,
  label = "All trials"
)

# 2. Excluding short follow-up trials (AZISTAT + DRI12544)
Multivariate_deltaattack_results_noshort <- run_multivariate_attack(
  data  = Data_Oracle_attack_noshort,
  label = "Exclude short follow-up trials (AZISTAT + DRI12544)"
)

#======================================================================
# === FOREST PLOTS — SHORT FOLLOW-UP SENSITIVITY ===
#======================================================================

base_path <- "/Users/joel/ORACLE - Placebo - Exacerbations/Project/"

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_all,
  label    = "All trials",
  filename = paste0(base_path,
                    "Forrest_plot_Attacks_Multivar_AllTrials.pdf")
)

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_noshort,
  label    = "Exclude short follow-up trials (AZISTAT + DRI12544)",
  filename = paste0(base_path,
                    "Forrest_plot_Attacks_Multivar_NoShortFollowup.pdf")
)

cat("\n✓ Short follow-up sensitivity multivariate models completed.\n")

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S11. Asthma attacks - Sensitivity analysis for zero-inflation pre-trial exacerbations 
# without 0 exacerbation at baseline
#======================================================================================================================================
#======================================================================
# === CREATE DATASETS: INCLUDE vs EXCLUDE ZERO BASELINE ATTACKS ===
#======================================================================

# A) All patients (primary analysis – includes zero baseline)
Data_Oracle_attack_all <- ORACLE_Placebo_Attack_imputed

# B) Excluding zero baseline attacks (sensitivity – floor effect removed)
Data_Oracle_attack_nozero <- ORACLE_Placebo_Attack_imputed %>%
  dplyr::filter(Attack_history_num > 0)

# === SANITY CHECK (optional but recommended)
cat("\nDataset sizes (imp = 1):\n")
cat("All patients              :", nrow(Data_Oracle_attack_all  %>% filter(.imp == 1)), "\n")
cat("Exclude zero baseline     :", nrow(Data_Oracle_attack_nozero %>% filter(.imp == 1)), "\n")


#======================================================================
# === RUN UNIVARIATE MODELS — ZERO BASELINE SENSITIVITY ===
#======================================================================

# 1. All patients (includes zero baseline)
Univariate_deltaattack_results_all <- run_univariate_attack(
  data  = Data_Oracle_attack_all,
  label = "All patients (including zero baseline)"
)

# 2. Excluding zero baseline attacks
Univariate_deltaattack_results_nozero <- run_univariate_attack(
  data  = Data_Oracle_attack_nozero,
  label = "Exclude zero baseline attacks"
)

#======================================================================
# === FOREST PLOTS — ZERO BASELINE SENSITIVITY ===
#======================================================================

base_path <- "/Users/joel/ORACLE - Placebo - Exacerbations/Project/"

plot_univariate_attack(
  results  = Univariate_deltaattack_results_all,
  label    = "All patients (incl. zero baseline)",
  filename = paste0(base_path,
                    "Forrest_plot_Attacks_Univar_AllPatients.pdf")
)

plot_univariate_attack(
  results  = Univariate_deltaattack_results_nozero,
  label    = "Exclude zero baseline attacks",
  filename = paste0(base_path,
                    "Forrest_plot_Attacks_Univar_NoZeroBaseline.pdf")
)

cat("\n✓ Zero baseline sensitivity univariate models completed.\n")


#======================================================================================================================================
# MULTIVARIATE ANALYSIS
#======================================================================================================================================

#======================================================================
# === RUN MULTIVARIATE MODELS — ZERO BASELINE SENSITIVITY ===
#======================================================================

# 1. All patients (includes zero baseline)
Multivariate_deltaattack_results_all <- run_multivariate_attack(
  data  = Data_Oracle_attack_all,
  label = "All patients (including zero baseline)"
)

# 2. Excluding zero baseline attacks
Multivariate_deltaattack_results_nozero <- run_multivariate_attack(
  data  = Data_Oracle_attack_nozero,
  label = "Exclude zero baseline attacks"
)

#======================================================================
# === FOREST PLOTS — ZERO BASELINE SENSITIVITY ===
#======================================================================

base_path <- "/Users/joel/ORACLE - Placebo - Exacerbations/Project/"

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_all,
  label    = "All patients (incl. zero baseline)",
  filename = paste0(base_path,
                    "Forrest_plot_Attacks_Multivar_AllPatients.pdf")
)

plot_multivariate_attack(
  results  = Multivariate_deltaattack_results_nozero,
  label    = "Exclude zero baseline attacks",
  filename = paste0(base_path,
                    "Forrest_plot_Attacks_Multivar_NoZeroBaseline.pdf")
)

cat("\n✓ Zero baseline sensitivity univariate models completed.\n")


#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S12 - SEE "Simple model with only FEV1_preBD_L_0W_per_10 as covariable" (PART B) in the PlaceboAsthma_MAIN_Script_clean
#======================================================================================================================================

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S13 - Lung function - Forest plot’s results per trial + I2 (for FeNO)
#======================================================================================================================================
# =============================================================================
# Per-trial forest plot: FeNO_log10 → Delta_FEV1
# Uses MI_estimates() filtered per trial + full-dataset pooled overall
# Adjusted for: FEV1_preBD_L_0W | Random intercept: No (per-trial)
# =============================================================================

library(dplyr)
library(ggplot2)
library(metafor)


# =============================================================================
# 0. SETTINGS
# =============================================================================

data         <- ORACLE_Placebo_Lung_imputed
outcome      <- "Placebo_response_Lung"
predictor    <- "FeNO_log10"
covariate    <- "FEV1_preBD_L_0W"
trial_col    <- "Enrolled_Trial_name"
imp_col      <- ".imp"
followup_col <- "Follow_up_duration_days"


# =============================================================================
# 1. GET TRIAL NAMES
# =============================================================================

trials <- data[[trial_col]][data[[imp_col]] == min(data[[imp_col]])] |>
  unique() |> sort()
cat("Trials:", paste(trials, collapse = ", "), "\n\n")


# =============================================================================
# 2. RUN MI_estimates PER TRIAL
# =============================================================================

run_trial <- function(trial_name) {
  cat("Fitting:", trial_name, "... ")
  
  df <- data[data[[trial_col]] == trial_name, ]
  n  <- nrow(df[df[[imp_col]] == min(df[[imp_col]]), ])
  
  res <- tryCatch(
    MI_estimates(
      data             = df,
      outcome_var      = outcome,
      predictor_vars   = predictor,
      covariables      = covariate,
      imp_col          = imp_col,
      followup_offset  = "No",
      followup_col     = followup_col,
      model_type       = "lm"
    ),
    error = function(e) { message("ERROR: ", e$message); NULL }
  )
  if (is.null(res)) return(NULL)
  
  # MI_estimates returns a data.frame when predictor_vars has 1 element
  # Extract the FeNO row
  row <- res[grepl(predictor, res$term), ]
  if (nrow(row) == 0) { message("  predictor row not found"); return(NULL) }
  
  cat("OK (n =", n, ")\n")
  data.frame(
    trial   = trial_name,
    n       = n,
    est     = row$estimate[1],
    se      = row$std.error[1],
    lo95    = row$`2.5 %`[1],
    hi95    = row$`97.5 %`[1],
    p.value = row$p.value[1],
    stringsAsFactors = FALSE
  )
}

results_list <- lapply(trials, run_trial)
results      <- do.call(rbind, Filter(Negate(is.null), results_list))


# =============================================================================
# 3. POOLED OVERALL — full dataset, random intercept on trial
# =============================================================================

cat("\nFitting overall pooled model...\n")

overall_fit <- MI_estimates(
  data                 = data,
  outcome_var          = outcome,
  predictor_vars       = predictor,
  covariables          = covariate,
  imp_col              = imp_col,
  followup_offset      = "No",
  followup_col         = followup_col,
  random_intercept_var = trial_col,
  model_type           = "lm"
)

overall_row <- overall_fit[grepl(predictor, overall_fit$term), ]

overall <- data.frame(
  trial   = "Overall",
  n       = sum(results$n),
  est     = overall_row$estimate[1],
  se      = overall_row$std.error[1],
  lo95    = overall_row$`2.5 %`[1],
  hi95    = overall_row$`97.5 %`[1],
  p.value = overall_row$p.value[1],
  stringsAsFactors = FALSE
)

cat(sprintf("Overall: \u03b2 = %.3f [%.3f, %.3f]  p = %s\n",
            overall$est, overall$lo95, overall$hi95,
            ifelse(overall$p.value < 0.001, "<0.001",
                   sprintf("%.3f", overall$p.value))))


# =============================================================================
# 4. HETEROGENEITY (Q-based I², Borenstein H method)
# =============================================================================

meta   <- rma(yi = est, sei = se, data = results, method = "REML")
k_meta <- nrow(results)
w_meta <- 1 / results$se^2
Q_meta <- sum(w_meta * (results$est - weighted.mean(results$est, w_meta))^2)
df     <- k_meta - 1
i2_est <- max(0, (Q_meta - df) / Q_meta)

H      <- sqrt(max(Q_meta, df) / df)
lnH    <- log(H)
se_lnH <- if (Q_meta > df + 1) (log(Q_meta) - log(df)) / (2*(Q_meta - df)) else sqrt(1/(2*df))
H_lo   <- exp(lnH - 1.96 * se_lnH)
H_hi   <- exp(lnH + 1.96 * se_lnH)
i2_lo  <- max(0,        1 - 1/H_lo^2)
i2_hi  <- max(0, min(1, 1 - 1/H_hi^2))

i2_label <- if (Q_meta <= df) {
  sprintf("Heterogeneity: I\u00b2 = 0 (Q(%d) = %.2f, p = %s) \u2014 no between-trial variance detected",
          df, Q_meta,
          ifelse(meta$QEp < 0.001, "<0.001", sprintf("%.3f", meta$QEp)))
} else {
  sprintf("Heterogeneity: I\u00b2 = %.2f [%.2f, %.2f]  |  Q(%d) = %.2f, p = %s",
          i2_est, i2_lo, i2_hi, df, Q_meta,
          ifelse(meta$QEp < 0.001, "<0.001", sprintf("%.3f", meta$QEp)))
}
cat("\n", i2_label, "\n")


# =============================================================================
# 5. PREPARE PLOT DATA
# =============================================================================

# Alphabetical, A at top → reverse for ggplot y-axis (bottom to top)
results <- results |>
  arrange(trial) |>
  mutate(
    weight_pct = n / sum(n) * 100,
    sig        = p.value < 0.05,
    p_label    = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    label      = sprintf("%.2f [%.2f, %.2f]   p=%s   W=%.1f%%",
                         est, lo95, hi95, p_label, weight_pct),
    trial      = factor(trial, levels = rev(trial))
  )

K         <- nrow(results)
trial_y   <- seq(K, 1)
overall_y <- -1

plot_trials <- results |> mutate(y = trial_y)

# Dot size proportional to weight (manual, no size legend)
min_pt <- 2; max_pt <- 7
w_min  <- min(results$weight_pct); w_max <- max(results$weight_pct)
plot_trials <- plot_trials |>
  mutate(pt_size = min_pt + (weight_pct - w_min) / (w_max - w_min) * (max_pt - min_pt))

overall_label <- sprintf("%.2f [%.2f, %.2f]   p=%s",
                         overall$est, overall$lo95, overall$hi95,
                         ifelse(overall$p.value < 0.001, "<0.001",
                                sprintf("%.3f", overall$p.value)))

x_lo    <- min(c(results$lo95, overall$lo95)) - 0.1
x_hi    <- max(c(results$hi95, overall$hi95)) + 0.1
label_x <- x_hi + (x_hi - x_lo) * 0.02

y_breaks <- c(trial_y, overall_y)
y_labels <- c(as.character(results$trial), "Overall")


# =============================================================================
# 6. FOREST PLOT
# =============================================================================

forest_plot <- ggplot() +
  
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  
  geom_hline(yintercept = (overall_y + 1) / 2,
             colour = "grey60", linewidth = 0.35) +
  
  # CI lines
  geom_errorbarh(
    data = plot_trials,
    aes(y = y, xmin = lo95, xmax = hi95),
    height = 0.18, linewidth = 0.5, colour = "grey30"
  ) +
  
  # Points (manual size, no legend)
  geom_point(
    data  = plot_trials,
    aes(x = est, y = y, colour = sig, fill = sig),
    size  = plot_trials$pt_size,
    shape = 21, stroke = 0.55
  ) +
  
  # Overall diamond
  geom_polygon(
    data = data.frame(
      x = c(overall$lo95, overall$est, overall$hi95, overall$est),
      y = c(overall_y, overall_y + 0.42, overall_y, overall_y - 0.42)
    ),
    aes(x = x, y = y),
    fill = "#185FA5", colour = "#0D3D6E", linewidth = 0.5
  ) +
  
  # Trial labels (right)
  geom_text(
    data = plot_trials,
    aes(x = label_x, y = y, label = label, colour = sig),
    hjust = 0, size = 2.7, show.legend = FALSE
  ) +
  
  # Overall label (right)
  annotate("text", x = label_x, y = overall_y, label = overall_label,
           hjust = 0, size = 2.7, colour = "#185FA5", fontface = "bold") +
  
  # I² annotation
  annotate("text",
           x = (x_lo + x_hi) / 2, y = overall_y - 0.85,
           label = i2_label,
           size = 2.5, colour = "grey35", hjust = 0.5, fontface = "italic") +
  
  # Column header
  annotate("text", x = label_x, y = K + 0.85,
           label = "\u03b2 [95% CI]   p-value   Weight",
           hjust = 0, size = 2.7, colour = "grey30", fontface = "bold") +
  
  scale_colour_manual(
    values = c("TRUE" = "#A32D2D", "FALSE" = "#4A4A47"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "p \u2265 0.05"),
    name   = NULL
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#D96060", "FALSE" = "#AAAAAA"),
    guide  = "none"
  ) +
  scale_x_continuous(
    limits = c(x_lo, x_hi + (x_hi - x_lo) * 0.52),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = y_labels,
    expand = expansion(add = c(1.4, 1.2))
  ) +
  
  labs(
    title    = "Effect of FeNO (log10) on Delta FEV1 (mL) - Per trial",
    subtitle = paste0(
      "Univariate linear regression, adjusted for ", covariate
    ),
    x = "Adjusted coefficient for FeNO (log10)  (95% CI)",
    y = NULL
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    panel.border       = element_rect(colour = "grey70", linewidth = 0.5),
    axis.text.y        = element_text(size = 12),
    axis.text.x        = element_text(size = 12),
    plot.title         = element_text(face = "bold", size = 20),
    plot.subtitle      = element_text(size = 8.5, colour = "grey40"),
    legend.position    = "bottom",
    legend.text        = element_text(size = 12),
    plot.margin        = margin(10, 10, 10, 10)
  )

print(forest_plot)


# =============================================================================
# 7. SAVE
# =============================================================================

plot_height <- max(6, K * 0.55 + 3.5)

ggsave(file.path(getwd(), "Forest_plot_FeNO_per_trial_FEV1_MIestimates.pdf"),
       plot = forest_plot, width = 18, height = plot_height, units = "in")

ggsave(file.path(getwd(), "Forest_plot_FeNO_per_trial_FEV1_MIestimates.png"),
       plot = forest_plot, width = 18, height = plot_height,
       dpi = 300, units = "in")

cat("\nSaved to:", getwd(), "\n")


# =============================================================================
# 8. RETURN
# =============================================================================

list(results = results, overall = overall,
     I2 = i2_est, I2_ci = c(i2_lo, i2_hi), i2_label = i2_label)

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S14. LUNG FUNCTION - Sensitivity analysis using ∆FEV₁ as % of predicted
#======================================================================================================================================
#======================================================================================================================================
# Figure S14.A - Simple model
#======================================================================================================================================
# Sensitivity analysis with Delta_FEV1 with Delta_%ofpredicted(PCT)
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
    "FEV1_preBD_L_0W_per_10",             # 14 - Baseline lung function
    "FEV1_reversibility_0W_per_10",        # 15 - Baseline lung function
    "BEC_log10",                           # 16 - Inflammatory biomarkers
    "FeNO_log10",                          # 17 - Inflammatory biomarkers
    "IgE_log10"                            # 18 - Inflammatory biomarkers
  ),
  covariables          = c("FEV1_preBD_L_0W_per_10"),
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
  "FEV1 pre-BD (L., per 0.1 decrease)",
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
  margin.left     = 18,
  margin.right    = 18,
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
    "FEV1_preBD_L_0W_per_10",             # 14 - Baseline lung function
    "FEV1_reversibility_0W_per_10",        # 15 - Baseline lung function
    "BEC_log10",                           # 16 - Inflammatory biomarkers
    "FeNO_log10",                          # 17 - Inflammatory biomarkers
    "IgE_log10"                            # 18 - Inflammatory biomarkers
  ),
  covariables = c(
    "Age_per_10",
    "Sex",
    "BMI_per_5",
    "FEV1_preBD_L_0W_per_10",
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
  "FEV1 pre-BD (L., per 0.1 decrease)",
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
  margin.left     = 18,
  margin.right    = 18,
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
#======================================================================================================================================
#======================================================================================================================================
# Figure S15. LUNG FUNCTION - Sensitivity analysis – CAPTAIN
# With CAPTAIN
#======================================================================================================================================
#SENSITIVITY ANALYSIS FOR LUNG FUNCTION - WITH CAPTAIN
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## For lung function change
#Selection of the studies with the total ∆FEV1
Studies_included_lung_withCAPTAIN<- c("BENRAP2B","LAVOLTA_1","LAVOLTA_2","NAVIGATOR","PATHWAY","QUEST","STRATOS_1","STRATOS_2", "CAPTAIN")

Data_Oracle_lung_withCAPTAIN <- ORACLE_after_imputation %>%
  filter(Enrolled_Trial_name%in%Studies_included_lung_withCAPTAIN)

nrow(Data_Oracle_lung_withCAPTAIN %>%
       filter(.imp==0))

Data_Oracle_lung_withCAPTAIN_imp0 <-Data_Oracle_lung_withCAPTAIN %>%
  filter(.imp==0)
sum(is.na(Data_Oracle_lung_withCAPTAIN_imp0$FEV1_preBD_L_0W))
sum(is.na(Data_Oracle_lung_withCAPTAIN_imp0$FEV1_preBD_L_52W))
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##CALCULATION OF DELTA FEV1 (calculated without imputated data)
Data_Oracle_lung_withCAPTAIN <- Data_Oracle_lung_withCAPTAIN %>%
  mutate(delta_FEV1_preBD_L= FEV1_preBD_L_52W-FEV1_preBD_L_0W)%>%
  mutate(delta_FEV1_preBD_mL= (FEV1_preBD_L_52W-FEV1_preBD_L_0W)*1000) %>%
  # Calculated with imputated data
  mutate(delta_FEV1_preBD_PCT= FEV1_preBD_PCT_52W-FEV1_preBD_PCT_0W)

colnames(Data_Oracle_lung)

#Put Treatment step in ordinal
Data_Oracle_lung_withCAPTAIN$Treatment_step_ordinal <- as.numeric(gsub("Step ", "", as.character(Data_Oracle_lung_withCAPTAIN$Treatment_step)))
summary(Data_Oracle_lung_withCAPTAIN$Treatment_step_ordinal)

# Sensitivity — with CAPTAIN 
Data_Oracle_lung_withCAPTAIN <- Data_Oracle_lung_withCAPTAIN%>%
  # Identify participants where Delta_FEV1 is missing in .imp == 0
  group_by(Sequential_number) %>%
  filter(!any(.imp == 0 & is.na(delta_FEV1_preBD_mL))) %>%
  ungroup() %>%
  # Keep only .imp = 1 to 10
  filter(.imp != 0)
unique(Data_Oracle_lung_withCAPTAIN$.imp)
sum(is.na(Data_Oracle_lung_withCAPTAIN$delta_FEV1_preBD_mL))
sum(is.na(Data_Oracle_lung_withCAPTAIN$Placebo_response_Lung))
Data_Oracle_lung_withCAPTAIN$Placebo_response_Lung

# === SANITY CHECK
cat("\nDataset sizes (imp = 1):\n")
cat("Main (without PACT & CAPTAIN) :", nrow(Data_Oracle_lung_imputated  %>% filter(.imp == 1)), "\n")
cat("With CAPTAIN           :", nrow(Data_Oracle_lung_withCAPTAIN %>% filter(.imp == 1)), "\n")

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# === UNIVARIATE MODEL (MODEL 1) - LUNG FUNCTION - WITH CAPTAIN ===

Univariate_deltalung_model_withCAPTAIN <- MI_estimates(
  data        = Data_Oracle_lung_withCAPTAIN,
  outcome_var = "Placebo_response_Lung",
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
    "FEV1_preBD_L_0W_per_10",           # 14 - Baseline lung function
    "FEV1_reversibility_0W_per_10",        # 15 - Baseline lung function
    "BEC_log10",                           # 16 - Inflammatory biomarkers
    "FeNO_log10",                          # 17 - Inflammatory biomarkers
    "IgE_log10"                            # 18 - Inflammatory biomarkers
  ),
  covariables          = c("FEV1_preBD_L_0W_per_10"),
  imp_col              = ".imp",
  followup_offset      = "No",
  followup_col         = "Follow_up_duration_days",
  random_intercept_var = "Enrolled_Trial_name",
  model_type           = "lm"
)

# === EXTRACT & LABEL RESULTS ===
Univariate_deltalung_results_withCAPTAIN <- attr(Univariate_deltalung_model_withCAPTAIN, "combined_results")

row.names(Univariate_deltalung_results_withCAPTAIN) <- c(
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
  "FEV1 pre-BD (L., per 0.1 decrease)",
  "FEV1 reversibility (per 10% increase)",
  # Inflammatory biomarkers (3)
  "Blood eosinophils (log10)",
  "FeNO (log10)",
  "Total IgE (log10)"
)

# === FOREST PLOT - UNIVARIATE WITH CAPTAIN ===
pdf(
  "/Users/joel/ORACLE - Placebo - Exacerbations/Project/Forrest_plot_Lung_Univar_WithCAPTAIN.pdf",
  width = 18, height = 12
)
par(mar = c(4, 2, 2, 2))

Forrest_plot_Lung_Univar_WithCAPTAIN <- forplo(
  as.data.frame(Univariate_deltalung_results_withCAPTAIN[, c("estimate", "2.5 %", "97.5 %")]),
  xlim            = c(-180, 180),
  em              = "aRC",
  linreg          = TRUE,
  row.labels      = row.names(Univariate_deltalung_results_withCAPTAIN),
  left.align      = FALSE,
  add.arrow.right = FALSE,
  arrow.right.length = 20,
  add.arrow.left  = FALSE,
  arrow.left.length  = 20,
  left.bar        = FALSE,
  shade.every     = 1,
  shade.col       = "grey",
  shade.alpha     = 0.2,
  margin.left     = 18,
  margin.right    = 18,
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

#--------------------------------------------------------------------------------------------------
# === MULTIVARIATE MODEL (MODEL 2) - LUNG FUNCTION - WITH CAPTAIN ===

Multivariate_deltalung_model_withCAPTAIN <- MI_estimates(
  data        = Data_Oracle_lung_withCAPTAIN,
  outcome_var = "Placebo_response_Lung",
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
    "FEV1_preBD_L_0W_per_10",           # 14 - Baseline lung function
    "FEV1_reversibility_0W_per_10",        # 15 - Baseline lung function
    "BEC_log10",                           # 16 - Inflammatory biomarkers
    "FeNO_log10",                          # 17 - Inflammatory biomarkers
    "IgE_log10"                            # 18 - Inflammatory biomarkers
  ),
  covariables = c(
    "Age_per_10",
    "Sex",
    "BMI_per_5",
    "FEV1_preBD_L_0W_per_10",
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
Multivariate_deltalung_results_withCAPTAIN <- attr(Multivariate_deltalung_model_withCAPTAIN, "combined_results")

row.names(Multivariate_deltalung_results_withCAPTAIN) <- c(
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
  "FEV1 pre-BD (L., per 0.1 decrease)",
  "FEV1 reversibility (per 10% increase)",
  # Inflammatory biomarkers (3)
  "Blood eosinophils (log10)",
  "FeNO (log10)",
  "Total IgE (log10)"
)

# === FOREST PLOT - MULTIVARIATE WITH CAPTAIN ===
pdf(
  "/Users/joel/ORACLE - Placebo - Exacerbations/Project/Forrest_plot_Lung_Multivar_WithCAPTAIN.pdf",
  width = 18, height = 12
)
par(mar = c(4, 2, 2, 2))

Forrest_plot_Lung_Multivar_WithCAPTAIN <- forplo(
  as.data.frame(Multivariate_deltalung_results_withCAPTAIN[, c("estimate", "2.5 %", "97.5 %")]),
  xlim            = c(-100, 150),
  em              = "aRC",
  linreg          = TRUE,
  row.labels      = row.names(Multivariate_deltalung_results_withCAPTAIN),
  left.align      = FALSE,
  add.arrow.right = FALSE,
  arrow.right.length = 20,
  add.arrow.left  = FALSE,
  arrow.left.length  = 20,
  left.bar        = FALSE,
  shade.every     = 1,
  shade.col       = "grey",
  shade.alpha     = 0.2,
  margin.left     = 18,
  margin.right    = 18,
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
#======================================================================================================================================
#======================================================================================================================================
# Figure S16 - SEE "Simple model with only ACQ_score_0W as covariable" (PART C) in the PlaceboAsthma_MAIN_Script_clean
#======================================================================================================================================

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S17 - ACQ-5 - Forest plot’s results per trial + I2 (for FeNO)
#======================================================================================================================================

# =============================================================================
# Per-trial forest plot: FeNO_log10 → Delta_ACQ
# Uses MI_estimates() filtered per trial + full-dataset pooled overall
# Adjusted for: ACQ_score_0W | Random intercept: No (per-trial)
# =============================================================================

library(dplyr)
library(ggplot2)
library(metafor)


# =============================================================================
# 0. SETTINGS
# =============================================================================

data         <- ORACLE_Placebo_ACQ_imputed
outcome      <- "Placebo_response_ACQ"
predictor    <- "FeNO_log10"
covariate    <- "ACQ_score_0W"
trial_col    <- "Enrolled_Trial_name"
imp_col      <- ".imp"
followup_col <- "Follow_up_duration_days"


# =============================================================================
# 1. GET TRIAL NAMES
# =============================================================================

trials <- data[[trial_col]][data[[imp_col]] == min(data[[imp_col]])] |>
  unique() |> sort()
cat("Trials:", paste(trials, collapse = ", "), "\n\n")


# =============================================================================
# 2. RUN MI_estimates PER TRIAL
# =============================================================================

run_trial <- function(trial_name) {
  cat("Fitting:", trial_name, "... ")
  
  df <- data[data[[trial_col]] == trial_name, ]
  n  <- nrow(df[df[[imp_col]] == min(df[[imp_col]]), ])
  
  res <- tryCatch(
    MI_estimates(
      data             = df,
      outcome_var      = outcome,
      predictor_vars   = predictor,
      covariables      = covariate,
      imp_col          = imp_col,
      followup_offset  = "No",
      followup_col     = followup_col,
      model_type       = "lm"
    ),
    error = function(e) { message("ERROR: ", e$message); NULL }
  )
  if (is.null(res)) return(NULL)
  
  # MI_estimates returns a data.frame when predictor_vars has 1 element
  # Extract the FeNO row
  row <- res[grepl(predictor, res$term), ]
  if (nrow(row) == 0) { message("  predictor row not found"); return(NULL) }
  
  cat("OK (n =", n, ")\n")
  data.frame(
    trial   = trial_name,
    n       = n,
    est     = row$estimate[1],
    se      = row$std.error[1],
    lo95    = row$`2.5 %`[1],
    hi95    = row$`97.5 %`[1],
    p.value = row$p.value[1],
    stringsAsFactors = FALSE
  )
}

results_list <- lapply(trials, run_trial)
results      <- do.call(rbind, Filter(Negate(is.null), results_list))


# =============================================================================
# 3. POOLED OVERALL — full dataset, random intercept on trial
# =============================================================================

cat("\nFitting overall pooled model...\n")

overall_fit <- MI_estimates(
  data                 = data,
  outcome_var          = outcome,
  predictor_vars       = predictor,
  covariables          = covariate,
  imp_col              = imp_col,
  followup_offset      = "No",
  followup_col         = followup_col,
  random_intercept_var = trial_col,
  model_type           = "lm"
)

overall_row <- overall_fit[grepl(predictor, overall_fit$term), ]

overall <- data.frame(
  trial   = "Overall",
  n       = sum(results$n),
  est     = overall_row$estimate[1],
  se      = overall_row$std.error[1],
  lo95    = overall_row$`2.5 %`[1],
  hi95    = overall_row$`97.5 %`[1],
  p.value = overall_row$p.value[1],
  stringsAsFactors = FALSE
)

cat(sprintf("Overall: \u03b2 = %.3f [%.3f, %.3f]  p = %s\n",
            overall$est, overall$lo95, overall$hi95,
            ifelse(overall$p.value < 0.001, "<0.001",
                   sprintf("%.3f", overall$p.value))))


# =============================================================================
# 4. HETEROGENEITY (Q-based I², Borenstein H method)
# =============================================================================

meta   <- rma(yi = est, sei = se, data = results, method = "REML")
k_meta <- nrow(results)
w_meta <- 1 / results$se^2
Q_meta <- sum(w_meta * (results$est - weighted.mean(results$est, w_meta))^2)
df     <- k_meta - 1
i2_est <- max(0, (Q_meta - df) / Q_meta)

H      <- sqrt(max(Q_meta, df) / df)
lnH    <- log(H)
se_lnH <- if (Q_meta > df + 1) (log(Q_meta) - log(df)) / (2*(Q_meta - df)) else sqrt(1/(2*df))
H_lo   <- exp(lnH - 1.96 * se_lnH)
H_hi   <- exp(lnH + 1.96 * se_lnH)
i2_lo  <- max(0,        1 - 1/H_lo^2)
i2_hi  <- max(0, min(1, 1 - 1/H_hi^2))

i2_label <- if (Q_meta <= df) {
  sprintf("Heterogeneity: I\u00b2 = 0 (Q(%d) = %.2f, p = %s) \u2014 no between-trial variance detected",
          df, Q_meta,
          ifelse(meta$QEp < 0.001, "<0.001", sprintf("%.3f", meta$QEp)))
} else {
  sprintf("Heterogeneity: I\u00b2 = %.2f [%.2f, %.2f]  |  Q(%d) = %.2f, p = %s",
          i2_est, i2_lo, i2_hi, df, Q_meta,
          ifelse(meta$QEp < 0.001, "<0.001", sprintf("%.3f", meta$QEp)))
}
cat("\n", i2_label, "\n")


# =============================================================================
# 5. PREPARE PLOT DATA
# =============================================================================

# Alphabetical, A at top → reverse for ggplot y-axis (bottom to top)
results <- results |>
  arrange(trial) |>
  mutate(
    weight_pct = n / sum(n) * 100,
    sig        = p.value < 0.05,
    p_label    = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    label      = sprintf("%.2f [%.2f, %.2f]   p=%s   W=%.1f%%",
                         est, lo95, hi95, p_label, weight_pct),
    trial      = factor(trial, levels = rev(trial))
  )

K         <- nrow(results)
trial_y   <- seq(K, 1)
overall_y <- -1

plot_trials <- results |> mutate(y = trial_y)

# Dot size proportional to weight (manual, no size legend)
min_pt <- 2; max_pt <- 7
w_min  <- min(results$weight_pct); w_max <- max(results$weight_pct)
plot_trials <- plot_trials |>
  mutate(pt_size = min_pt + (weight_pct - w_min) / (w_max - w_min) * (max_pt - min_pt))

overall_label <- sprintf("%.2f [%.2f, %.2f]   p=%s",
                         overall$est, overall$lo95, overall$hi95,
                         ifelse(overall$p.value < 0.001, "<0.001",
                                sprintf("%.3f", overall$p.value)))

x_lo    <- min(c(results$lo95, overall$lo95)) - 0.1
x_hi    <- max(c(results$hi95, overall$hi95)) + 0.1
label_x <- x_hi + (x_hi - x_lo) * 0.02

y_breaks <- c(trial_y, overall_y)
y_labels <- c(as.character(results$trial), "Overall")


# =============================================================================
# 6. FOREST PLOT
# =============================================================================

forest_plot <- ggplot() +
  
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  
  geom_hline(yintercept = (overall_y + 1) / 2,
             colour = "grey60", linewidth = 0.35) +
  
  # CI lines
  geom_errorbarh(
    data = plot_trials,
    aes(y = y, xmin = lo95, xmax = hi95),
    height = 0.18, linewidth = 0.5, colour = "grey30"
  ) +
  
  # Points (manual size, no legend)
  geom_point(
    data  = plot_trials,
    aes(x = est, y = y, colour = sig, fill = sig),
    size  = plot_trials$pt_size,
    shape = 21, stroke = 0.55
  ) +
  
  # Overall diamond
  geom_polygon(
    data = data.frame(
      x = c(overall$lo95, overall$est, overall$hi95, overall$est),
      y = c(overall_y, overall_y + 0.42, overall_y, overall_y - 0.42)
    ),
    aes(x = x, y = y),
    fill = "#185FA5", colour = "#0D3D6E", linewidth = 0.5
  ) +
  
  # Trial labels (right)
  geom_text(
    data = plot_trials,
    aes(x = label_x, y = y, label = label, colour = sig),
    hjust = 0, size = 2.7, show.legend = FALSE
  ) +
  
  # Overall label (right)
  annotate("text", x = label_x, y = overall_y, label = overall_label,
           hjust = 0, size = 2.7, colour = "#185FA5", fontface = "bold") +
  
  # I² annotation
  annotate("text",
           x = (x_lo + x_hi) / 2, y = overall_y - 0.85,
           label = i2_label,
           size = 2.5, colour = "grey35", hjust = 0.5, fontface = "italic") +
  
  # Column header
  annotate("text", x = label_x, y = K + 0.85,
           label = "\u03b2 [95% CI]   p-value   Weight",
           hjust = 0, size = 2.7, colour = "grey30", fontface = "bold") +
  
  scale_colour_manual(
    values = c("TRUE" = "#A32D2D", "FALSE" = "#4A4A47"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "p \u2265 0.05"),
    name   = NULL
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#D96060", "FALSE" = "#AAAAAA"),
    guide  = "none"
  ) +
  scale_x_continuous(
    limits = c(x_lo, x_hi + (x_hi - x_lo) * 0.52),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = y_labels,
    expand = expansion(add = c(1.4, 1.2))
  ) +
  
  labs(
    title    = "Effect of FeNO (log10) on Delta ACQ-5 (score) - Per trial",
    subtitle = paste0(
      "Univariate linear regression, adjusted for ", covariate
    ),
    x = "Adjusted coefficient for FeNO (log10)  (95% CI)",
    y = NULL
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    panel.border       = element_rect(colour = "grey70", linewidth = 0.5),
    axis.text.y        = element_text(size = 12),
    axis.text.x        = element_text(size = 12),
    plot.title         = element_text(face = "bold", size = 20),
    plot.subtitle      = element_text(size = 12, colour = "grey40"),
    legend.position    = "bottom",
    legend.text        = element_text(size = 12),
    plot.margin        = margin(10, 10, 10, 10)
  )

print(forest_plot)


# =============================================================================
# 7. SAVE
# =============================================================================

plot_height <- max(6, K * 0.55 + 3.5)

ggsave(file.path(getwd(), "Forest_plot_FeNO_per_trial_ACQ5_MIestimates.pdf"),
       plot = forest_plot, width = 18, height = plot_height, units = "in")

ggsave(file.path(getwd(), "Forest_plot_FeNO_per_trial_ACQ5_MIestimates.png"),
       plot = forest_plot, width = 18, height = plot_height,
       dpi = 300, units = "in")

cat("\nSaved to:", getwd(), "\n")


# =============================================================================
# 8. RETURN
# =============================================================================

list(results = results, overall = overall,
     I2 = i2_est, I2_ci = c(i2_lo, i2_hi), i2_label = i2_label)



#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S18. Comparison with Results from Original imputation (Meulmeester et al., Lancet Res, 2025)
# Old imputation vs. NEW imputation (by SML)
#======================================================================================================================================
## IT HAS BEEN GENERATED PREVIOUSLY WITH MI_estimates
# TO DO IT AGAIN THE METHOD IS DESCRIBED BELOW
#======================================================================================================================================
#Import the DATA from the file
setwd("/Users/joel/ORACLE - Placebo - Exacerbations/Raw_data")
#setwd(paste0("/Users/",user_name,"/Library/CloudStorage/OneDrive-USherbrooke/Recherche/Projets/Pneumologie/ORACLE/ORACLE_RAW"))

## Importing the original data (will be used to create table 1)
col_types <- c("guess", "guess", "text", "text", "numeric", "text", "numeric", "numeric", "numeric", "text", "numeric", "text", "numeric", "numeric", "numeric", "text", "text", "text", "text", "text", "text", "text", "numeric", "text", "numeric", "text", "numeric", "numeric", "numeric", "text", "numeric", "text", "numeric", "numeric", "text", "text", "text", "text", "numeric", "numeric", "text", "numeric", "numeric", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "text", "numeric", "text", "numeric", "numeric", "text", "numeric", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "text", "numeric", "numeric", "numeric", "numeric")
data_original_imported <- suppressWarnings(read_excel("data_ORACLE_original_20240429.xlsx", col_types = col_types))

## Importing the imputated data without the systematically missing data (d)= Only imputation of of all missing data
col_types <-c("guess", "text" ,"text","text", "text", "text", "numeric", "numeric", "numeric", "text", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "text", "numeric", "text", "numeric", "numeric", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric","numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "text" )
Data_imputated_all_imported <- suppressWarnings(read_excel("data_ORACLE_imp_20250102_joined.xlsx", col_types = col_types))

## Importing the imputated data with the systematically missing data not replaced (NR) = Only imputation of non-systematically missing
col_types <-c("guess","text","text","text","text","text","numeric","text","numeric","text","text","text","text","numeric","numeric","text","text","numeric","text","text","text","text","text","text","text","text","numeric","numeric","numeric","text","numeric","text","text","text","text","numeric","text","text","text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","text","numeric","text","numeric","numeric","text","numeric","numeric","text","text","text","text","text","text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","text")
Data_imputated_without_systematically_missing_imported<-suppressWarnings(read_excel("data_ORACLE_imp_sysREMOVED_20250102_joined.xlsx", col_types = col_types))

data_original <- data_original_imported
data_imputated<- Data_imputated_all_imported
data_imputated_without_systematically_missing<- Data_imputated_without_systematically_missing_imported
colnames(data_imputated_without_systematically_missing)

# Create a vector with "numeric" for all columns
col_types <- rep("numeric", 260)
# Set the first and third columns to "guess"
col_types[c(1, 3)] <- "guess"
Raw_ACQ_FEV1 <- suppressWarnings(read_excel("ACQ_FEV1.DataSet.xlsx", col_types = col_types))
Raw_ACQ_FEV1<-as.data.frame(Raw_ACQ_FEV1)
#======================================================================================================================================
#Selecting the current working project for the project
##Selecting the working directory
setwd("/Users/joel/ORACLE - Placebo - Exacerbations/Project")
#setwd(paste0("/Users/",user_name,"/Library/CloudStorage/OneDrive-USherbrooke/Recherche/Projets/Pneumologie/ORACLE/ORACLE_Placebo"))
#======================================================================================================================================
#TRANSFORMATION of Values
#Putting NA for the missing value
data_original[data_original == "NA"] <- NA
data_imputated[data_imputated== "NA"] <- NA
data_imputated_without_systematically_missing[data_imputated_without_systematically_missing == "NA"] <- NA

#Transform values to put them in the the good units
#Modify the values of the stratos trials to put the reversibility in % for all the studies instead of decimal
data_original<- data_original %>%
  mutate(FEV1_reversibility_percent_postBD_real = case_when(
    Enrolled_Trial_name=="STRATOS_1" ~ FEV1_PCT_reversibility_postBD*100,
    Enrolled_Trial_name=="STRATOS_2" ~ FEV1_PCT_reversibility_postBD*100,
    TRUE ~ FEV1_PCT_reversibility_postBD))

#Modify the values of the Captain study to put adherence in trial in percentage everywhere
data_original<- data_original %>%
  mutate(Adherence_InTrial_quantity_real = case_when(
    Enrolled_Trial_name=="CAPTAIN" ~ Adherence_InTrial_quantity*100,
    TRUE ~ Adherence_InTrial_quantity))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#======================================================================================================================================
###Creation of a unique dataframe which contain original value, imputated and imputated_only_not_systematically_missing
###VÉRIFIER AVEC SML???

colnames(data_original)
unique(data_original$Subject_ID)

data_original_main<-data_original[,c("Sequential_number","Enrolled_Trial_name","Treatment_arm","Age","Gender_0Female_1Male","BMI","Ethnicity","Country","Region","Treatment_step","Any_severe_attack_previous_12m_0no_1yes","Any_attack_or_hospitalization_previous_12_months","Number_severe_attack_previous_12m","Number_hospitalisations_for_asthma_previous_12_months","Number_hospitalizations_previous_12m_1Yes_0No","Previous_ICU_0no_1yes_9999notknown","Previous_Intubation_0no_1yes_9999notknown","Previous_ICU_or_intubation_0no_1yes","Smoking_0never_1ex_2current","Pack_years","Psychiatric_disease_0no_1yes_9999notknown","Atopy_history_0no_1yes_9999notknown","Eczema_0no_1yes_9999notknown","AllergicRhinitis__0no_1yes_9999notknown","Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown","Chronic_Rhinosinusitis_0no_1yes_9999notknown","Nasal_polyposis_0no_1yes_9999notknown","Previous_nasal_polypectomy_0no_1yes_9999notknown","ICS_DOSE_CLASS","LABA_prescribed_0no_1yes","LAMA_prescribed__0no_1yes","maintenance_OCS_prescribed__0no_1yes","Theophylline_prescribed__0no_1yes","Intranasal_seroid_prescribed__0no_1yes","FEV1_predicted_L","FVC_predicted_L","FEV1_preBD_L_Baseline","FEV1_preBD_PCT_Baseline","FVC_preBD_L_Baseline","FEV1_postBD_L_Baseline","FEV1_postBD_PCT_Baseline","FVC_postBD_L_Baseline","FVC_postBD_PCT_Baseline","FEV1_PCT_reversibility_postBD","FEV1_FVC_ratio","ACQ_baseline_score_mean","ACT_baseline_score","Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced","FeNO_baseline_ppb","Total_IgE","Follow_up_duration_days_nozero","Number_severe_asthma_attacks_during_followup","Time_to_First_attack","Time_to_2n_attack","Time_to_3n_attack","Time_to_4n_attack","Time_to_5n_attack","End_FollowUp_Reason","FEV1PREBD_L_52W","FEV1PREBD_PCT_52W","FEV1POSTBD_L_52W","FEV1POSTBD_PCT_52W","FEV1_reversibility_percent_postBD_real")]

data_imputated_all<-data_imputated[,c("Sequential_number",".imp","Age","Gender_0Female_1Male","BMI","Any_severe_attack_previous_12m_0no_1yes","Number_severe_attack_previous_12m_con","Number_hospitalisations_for_asthma_previous_12_months_con","Previous_ICU_or_intubation_0no_1yes","Smoking_0never_1ex_2current","Pack_years","Atopy_history_0no_1yes_9999notknown","Eczema_0no_1yes_9999notknown","AllergicRhinitis__0no_1yes_9999notknown","Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown","Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown","Chronic_Rhinosinusitis_0no_1yes_9999notknown","Nasal_polyposis_0no_1yes_9999notknown","Previous_nasal_polypectomy_0no_1yes_9999notknown","FEV1_preBD_L_Baseline","FEV1_preBD_PCT_Baseline","FVC_preBD_L_Baseline","FEV1_postBD_L_Baseline","FEV1_postBD_PCT_Baseline","FVC_postBD_L_Baseline","FEV1_PCT_reversibility_postBD","FEV1_FVC_ratio","ACQ_baseline_score_mean","Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced","FeNO_baseline_ppb","Total_IgE")]
colnames(data_imputated_all)[-c(1:2)] <- paste0(colnames(data_imputated_all)[-c(1:2)], "_imputated")
colnames(data_imputated_all) <- make.unique(colnames(data_imputated_all))
colnames(data_imputated_all)

data_imputated_no_systematically_missing<-data_imputated_without_systematically_missing[,c("Sequential_number",".imp","Age","Gender_0Female_1Male","BMI","Any_severe_attack_previous_12m_0no_1yes","Number_severe_attack_previous_12m_con","Number_hospitalisations_for_asthma_previous_12_months_con","Previous_ICU_or_intubation_0no_1yes","Smoking_0never_1ex_2current","Pack_years","Atopy_history_0no_1yes_9999notknown","Eczema_0no_1yes_9999notknown","AllergicRhinitis__0no_1yes_9999notknown","Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown","Chronic_Rhinosinusitis_0no_1yes_9999notknown","Nasal_polyposis_0no_1yes_9999notknown","Previous_nasal_polypectomy_0no_1yes_9999notknown","FEV1_preBD_L_Baseline","FEV1_preBD_PCT_Baseline","FVC_preBD_L_Baseline","FEV1_postBD_L_Baseline","FEV1_postBD_PCT_Baseline","FVC_postBD_L_Baseline","FEV1_PCT_reversibility_postBD","FEV1_FVC_ratio","ACQ_baseline_score_mean","Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced","FeNO_baseline_ppb","Total_IgE")]
data_imputated_no_systematically_missing
colnames(data_imputated_no_systematically_missing)[-c(1:2)] <- paste0(colnames(data_imputated_no_systematically_missing)[-c(1:2)], "_imputated_no_systematically_missing")
colnames(data_imputated_no_systematically_missing)

unique(data_imputated_all$.imp)

#selection of the 8 first imputation for both dataset
#data_imputated_all<-data_imputated_all %>%
#filter(.imp<9)
#data_imputated_no_systematically_missing<- data_imputated_no_systematically_missing %>%
#filter(.imp<9)
nrow(data_imputated_all)
nrow(data_imputated_no_systematically_missing)

#combining the datasets
merged_data_imputated <- merge(data_imputated_all, data_imputated_no_systematically_missing,
                               by = c("Sequential_number", ".imp"),
                               all = TRUE)

All_data <- merge(data_original, merged_data_imputated,
                  by = c("Sequential_number"),
                  all = TRUE)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Removing the trials open label
Excluded_trials <- c("Novel_START", "PRACTICAL")
unique (All_data$Enrolled_Trial_name)

All_data__without_open_label<- All_data %>%
  filter(!Enrolled_Trial_name %in% Excluded_trials )
nrow(All_data__without_open_label)
unique (All_data__without_open_label$Enrolled_Trial_name)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Identify the categorical data
Data_Oracle<-All_data__without_open_label%>%
  #Gender
  mutate(Gender_0Female_1Male= case_when(Gender_0Female_1Male ==0 ~ "Female", Gender_0Female_1Male ==1 ~ "Male",TRUE ~ NA)) %>%
  mutate(Gender_0Female_1Male_imputated= case_when(Gender_0Female_1Male_imputated ==0 ~ "Female", Gender_0Female_1Male_imputated ==1 ~ "Male",TRUE ~ NA)) %>%
  mutate(Gender_0Female_1Male_imputated_no_systematically_missing= case_when(Gender_0Female_1Male_imputated_no_systematically_missing ==0 ~ "Female", Gender_0Female_1Male_imputated_no_systematically_missing ==1 ~ "Male",TRUE ~ NA)) %>%
  #Smoking
  mutate(Smoking_0never_1ex_2current= case_when(Smoking_0never_1ex_2current ==0 ~ "Never", Smoking_0never_1ex_2current ==1 ~ "Yes (current or ex)",Smoking_0never_1ex_2current ==2 ~ "Yes (current or ex)",TRUE ~ NA)) %>%
  mutate(Smoking_0never_1ex_2current_imputated= case_when(Smoking_0never_1ex_2current_imputated ==0 ~ "Never", Smoking_0never_1ex_2current_imputated ==1 ~ "Yes (current or ex)",Smoking_0never_1ex_2current_imputated ==2 ~ "Yes (current or ex)",TRUE ~ NA)) %>%
  mutate(Smoking_0never_1ex_2current_imputated_no_systematically_missing= case_when(Smoking_0never_1ex_2current_imputated_no_systematically_missing ==0 ~ "Never", Smoking_0never_1ex_2current_imputated_no_systematically_missing ==1 ~ "Yes (current or ex)",Smoking_0never_1ex_2current_imputated_no_systematically_missing ==2 ~ "Yes (current or ex)",TRUE ~ NA)) %>%
  #Atopy
  mutate(Atopy_history_0no_1yes_9999notknown= case_when(Atopy_history_0no_1yes_9999notknown ==0 ~ "No", Atopy_history_0no_1yes_9999notknown ==1 ~ "Yes",TRUE ~ NA)) %>%
  mutate(Atopy_history_0no_1yes_9999notknown_imputated= case_when(Atopy_history_0no_1yes_9999notknown_imputated ==0 ~ "No", Atopy_history_0no_1yes_9999notknown_imputated ==1 ~ "Yes",TRUE ~ NA)) %>%
  mutate(Atopy_history_0no_1yes_9999notknown_imputated_no_systematically_missing= case_when(Atopy_history_0no_1yes_9999notknown_imputated_no_systematically_missing ==0 ~ "No", Atopy_history_0no_1yes_9999notknown_imputated_no_systematically_missing ==1 ~ "Yes",TRUE ~ NA)) %>%
  #Airborne_allergen_sensitisation
  mutate(Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown= case_when(Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown ==0 ~ "No", Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown ==1 ~ "Yes",TRUE ~ NA)) %>%
  mutate(Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown_imputated= case_when(Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown_imputated ==0 ~ "No", Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown_imputated ==1 ~ "Yes",TRUE ~ NA)) %>%
  mutate(Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown_imputated_no_systematically_missing= case_when(Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown_imputated_no_systematically_missing ==0 ~ "No", Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown_imputated_no_systematically_missing ==1 ~ "Yes",TRUE ~ NA)) %>%
  #Eczema
  mutate(Eczema_0no_1yes_9999notknown= case_when(Eczema_0no_1yes_9999notknown ==0 ~ "No", Eczema_0no_1yes_9999notknown ==1 ~ "Yes",TRUE ~ NA)) %>%
  mutate(Eczema_0no_1yes_9999notknown_imputated= case_when(Eczema_0no_1yes_9999notknown_imputated ==0 ~ "No", Eczema_0no_1yes_9999notknown_imputated ==1 ~ "Yes",TRUE ~ NA)) %>%
  mutate(Eczema_0no_1yes_9999notknown_imputated_no_systematically_missing= case_when(Eczema_0no_1yes_9999notknown_imputated_no_systematically_missing ==0 ~ "No", Eczema_0no_1yes_9999notknown_imputated_no_systematically_missing ==1 ~ "Yes",TRUE ~ NA)) %>%
  #Allergic rhinitis
  mutate(AllergicRhinitis__0no_1yes_9999notknown= case_when(AllergicRhinitis__0no_1yes_9999notknown ==0 ~ "No", AllergicRhinitis__0no_1yes_9999notknown ==1 ~ "Yes",TRUE ~ NA)) %>%
  mutate(AllergicRhinitis__0no_1yes_9999notknown_imputated= case_when(AllergicRhinitis__0no_1yes_9999notknown_imputated ==0 ~ "No", AllergicRhinitis__0no_1yes_9999notknown_imputated ==1 ~ "Yes",TRUE ~ NA)) %>%
  mutate(AllergicRhinitis__0no_1yes_9999notknown_imputated_no_systematically_missing= case_when(AllergicRhinitis__0no_1yes_9999notknown_imputated_no_systematically_missing ==0 ~ "No", AllergicRhinitis__0no_1yes_9999notknown_imputated_no_systematically_missing ==1 ~ "Yes",TRUE ~ NA)) %>%
  #Chronic rhinosinusitis
  mutate(Chronic_Rhinosinusitis_0no_1yes_9999notknown= case_when(Chronic_Rhinosinusitis_0no_1yes_9999notknown ==0 ~ "No", Chronic_Rhinosinusitis_0no_1yes_9999notknown ==1 ~ "Yes",TRUE ~ NA)) %>%
  mutate(Chronic_Rhinosinusitis_0no_1yes_9999notknown_imputated= case_when(Chronic_Rhinosinusitis_0no_1yes_9999notknown_imputated ==0 ~ "No", Chronic_Rhinosinusitis_0no_1yes_9999notknown_imputated ==1 ~ "Yes",TRUE ~ NA)) %>%
  mutate(Chronic_Rhinosinusitis_0no_1yes_9999notknown_imputated_no_systematically_missing= case_when(Chronic_Rhinosinusitis_0no_1yes_9999notknown_imputated_no_systematically_missing ==0 ~ "No", Chronic_Rhinosinusitis_0no_1yes_9999notknown_imputated_no_systematically_missing ==1 ~ "Yes",TRUE ~ NA)) %>%
  #Nasal polyposis
  mutate(Nasal_polyposis_0no_1yes_9999notknown= case_when(Nasal_polyposis_0no_1yes_9999notknown ==0 ~ "No", Nasal_polyposis_0no_1yes_9999notknown ==1 ~ "Yes",TRUE ~ NA)) %>%
  mutate(Nasal_polyposis_0no_1yes_9999notknown_imputated= case_when(Nasal_polyposis_0no_1yes_9999notknown_imputated ==0 ~ "No", Nasal_polyposis_0no_1yes_9999notknown_imputated ==1 ~ "Yes",TRUE ~ NA)) %>%
  mutate(Nasal_polyposis_0no_1yes_9999notknown_imputated_no_systematically_missing= case_when(Nasal_polyposis_0no_1yes_9999notknown_imputated_no_systematically_missing ==0 ~ "No", Nasal_polyposis_0no_1yes_9999notknown_imputated_no_systematically_missing ==1 ~ "Yes",TRUE ~ NA))


colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Gender_0Female_1Male" )]<- "Gender"
Data_Oracle$Gender<-as.factor(Data_Oracle$Gender)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Gender_0Female_1Male_imputated" )]<- "Gender_imputated"
Data_Oracle$Gender_imputated<-as.factor(Data_Oracle$Gender_imputated)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Gender_0Female_1Male_imputated_no_systematically_missing" )]<- "Gender_imputated_no_systematically_missing"
Data_Oracle$Gender_imputated_no_systematically_missing<-as.factor(Data_Oracle$Gender_imputated_no_systematically_missing)

colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Smoking_0never_1ex_2current" )]<- "Smoking_Statut"
Data_Oracle$Smoking_Statut<-as.factor(Data_Oracle$Smoking_Statut)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Smoking_0never_1ex_2current_imputated" )]<- "Smoking_Statut_imputated"
Data_Oracle$Smoking_Statut_imputated<-as.factor(Data_Oracle$Smoking_Statut_imputated)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Smoking_0never_1ex_2current_imputated_no_systematically_missing" )]<- "Smoking_Statut_imputated_no_systematically_missing"
Data_Oracle$Smoking_Statut_imputated_no_systematically_missing<-as.factor(Data_Oracle$Smoking_Statut_imputated_no_systematically_missing)

colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Atopy_history_0no_1yes_9999notknown" )]<- "Atopy_history"
Data_Oracle$Atopy_history<-as.factor(Data_Oracle$Atopy_history)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Atopy_history_0no_1yes_9999notknown_imputated" )]<- "Atopy_history_imputated"
Data_Oracle$Atopy_history_imputated<-as.factor(Data_Oracle$Atopy_history_imputated)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Atopy_history_0no_1yes_9999notknown_imputated_no_systematically_missing" )]<- "Atopy_history_imputated_no_systematically_missing"
Data_Oracle$Atopy_history_imputated_no_systematically_missing<-as.factor(Data_Oracle$Atopy_history_imputated_no_systematically_missing)

colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown" )]<- "Airborne_allergen_sensibilisation"
Data_Oracle$Airborne_allergen_sensibilisation<-as.factor(Data_Oracle$Airborne_allergen_sensibilisation)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown_imputated" )]<- "Airborne_allergen_sensibilisation_imputated"
Data_Oracle$Airborne_allergen_sensibilisation_imputated<-as.factor(Data_Oracle$Airborne_allergen_sensibilisation_imputated)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Airborne_allergen_sensitisation_on_testing_0no_1yes_9999notknown_imputated_no_systematically_missing" )]<- "Airborne_allergen_sensibilisation_imputated_no_systematically_missing"
Data_Oracle$Airborne_allergen_sensibilisation_imputated_no_systematically_missing<-as.factor(Data_Oracle$Airborne_allergen_sensibilisation_imputated_no_systematically_missing)

colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Eczema_0no_1yes_9999notknown" )]<- "Eczema"
Data_Oracle$Eczema<-as.factor(Data_Oracle$Eczema)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Eczema_0no_1yes_9999notknown_imputated" )]<- "Eczema_imputated"
Data_Oracle$Eczema_imputated<-as.factor(Data_Oracle$Eczema_imputated)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Eczema_0no_1yes_9999notknown_imputated_no_systematically_missing" )]<- "Eczema_imputated_no_systematically_missing"
Data_Oracle$Eczema_imputated_no_systematically_missing<-as.factor(Data_Oracle$Eczema_imputated_no_systematically_missing)

colnames(Data_Oracle)[which(colnames(Data_Oracle)=="AllergicRhinitis__0no_1yes_9999notknown" )]<- "Allergic_rhinitis"
Data_Oracle$Allergic_rhinitis<-as.factor(Data_Oracle$Allergic_rhinitis)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="AllergicRhinitis__0no_1yes_9999notknown_imputated" )]<- "Allergic_rhinitis_imputated"
Data_Oracle$Allergic_rhinitis_imputated<-as.factor(Data_Oracle$Allergic_rhinitis_imputated)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="AllergicRhinitis__0no_1yes_9999notknown_imputated_no_systematically_missing" )]<- "Allergic_rhinitis_imputated_no_systematically_missing"
Data_Oracle$Allergic_rhinitis_imputated_no_systematically_missing<-as.factor(Data_Oracle$Allergic_rhinitis_imputated_no_systematically_missing)

colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Chronic_Rhinosinusitis_0no_1yes_9999notknown" )]<- "Chronic_rhinosinusitis"
Data_Oracle$Chronic_rhinosinusitis<-as.factor(Data_Oracle$Chronic_rhinosinusitis)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Chronic_Rhinosinusitis_0no_1yes_9999notknown_imputated" )]<- "Chronic_rhinosinusitis_imputated"
Data_Oracle$Chronic_rhinosinusitis_imputated<-as.factor(Data_Oracle$Chronic_rhinosinusitis_imputated)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Chronic_Rhinosinusitis_0no_1yes_9999notknown_imputated_no_systematically_missing" )]<- "Chronic_rhinosinusitis_imputated_no_systematically_missing"
Data_Oracle$Chronic_rhinosinusitis_imputated_no_systematically_missing<-as.factor(Data_Oracle$Chronic_rhinosinusitis_imputated_no_systematically_missing)

colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Nasal_polyposis_0no_1yes_9999notknown" )]<- "Nasal_polyposis"
Data_Oracle$Nasal_polyposis<-as.factor(Data_Oracle$Nasal_polyposis)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Nasal_polyposis_0no_1yes_9999notknown_imputated" )]<- "Nasal_polyposis_imputated"
Data_Oracle$Nasal_polyposis_imputated<-as.factor(Data_Oracle$Nasal_polyposis_imputated)
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Nasal_polyposis_0no_1yes_9999notknown_imputated_no_systematically_missing" )]<- "Nasal_polyposis_imputated_no_systematically_missing"
Data_Oracle$Nasal_polyposis_imputated_no_systematically_missing<-as.factor(Data_Oracle$Nasal_polyposis_imputated_no_systematically_missing)

## Modify the name of columns for clear names
Data_Oracle<-Data_Oracle %>%
  mutate(Eosinophils_Log=log10(Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced))
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced_imputated" )]<- "Eosinophils_Log_imputated"
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced_imputated_no_systematically_missing" )]<- "Eosinophils_Log_imputated_no_systematically_missing"

Data_Oracle<-Data_Oracle %>%
  mutate(FeNO_Log=log10(FeNO_baseline_ppb))
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="FeNO_baseline_ppb_imputated" )]<- "FeNO_Log_imputated"
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="FeNO_baseline_ppb_imputated_no_systematically_missing" )]<- "FeNO_Log_imputated_no_systematically_missing"

Data_Oracle<-Data_Oracle %>%
  mutate(IgE_Log=log10(Total_IgE))
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Total_IgE_imputated")]<- "IgE_Log_imputated"
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Total_IgE_imputated_no_systematically_missing")]<- "IgE_Log_imputated_no_systematically_missing"

colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Number_severe_attack_previous_12m_con"  )]<- "Attack_12mo_Nb"
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Number_severe_attack_previous_12m_con_imputated"  )]<- "Attack_12mo_Nb_imputated"
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Number_severe_attack_previous_12m_con_imputated_no_systematically_missing"  )]<- "Attack_12mo_Nb_imputated_no_systematically_missing"

colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Number_hospitalisations_for_asthma_previous_12_months_con"  )]<- "Hospitalisations_12mo_Nb"
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Number_hospitalisations_for_asthma_previous_12_months_con_imputated"  )]<- "Hospitalisations_12mo_Nb_imputated"
colnames(Data_Oracle)[which( colnames(Data_Oracle)=="Number_hospitalisations_for_asthma_previous_12_months_con_imputated_no_systematically_missing"  )]<- "Hospitalisations_12mo_Nb_imputated_no_systematically_missing"
#======================================================================================================================================
#Calculation of the continuous variable

##Calculation of predicted spirometric parameters according to %FEV1 or %FVC when the FEV1_predicted_L was not given
Data_Oracle<-Data_Oracle%>%
  mutate(FEV1_predicted_L= case_when(
    is.na(FEV1_predicted_L) & !is.na(FEV1_preBD_L_Baseline) & !is.na(FEV1_preBD_PCT_Baseline) ~ (100*FEV1_preBD_L_Baseline)/FEV1_preBD_PCT_Baseline,
    is.na(FEV1_predicted_L) & !is.na(FEV1_postBD_L_Baseline) & !is.na(FEV1_postBD_PCT_Baseline) ~ (100*FEV1_postBD_L_Baseline)/FEV1_postBD_PCT_Baseline,
    TRUE ~ FEV1_predicted_L
  )) %>%
  mutate(FVC_predicted_L= case_when(
    is.na(FVC_predicted_L) & !is.na(FVC_preBD_L_Baseline) & !is.na(FVC_preBD_PCT_Baseline) ~ (100*FVC_preBD_L_Baseline)/FVC_preBD_PCT_Baseline,
    is.na(FVC_predicted_L) & !is.na(FVC_postBD_L_Baseline) & !is.na(FVC_postBD_PCT_Baseline) ~ (100*FVC_postBD_L_Baseline)/FVC_postBD_PCT_Baseline,
    TRUE ~ FVC_predicted_L
  ))
colnames(Data_Oracle)

## Calculation FEV1/FVC

Data_Oracle<-Data_Oracle%>%
  mutate(Tiffeneau=FEV1_preBD_L_Baseline/FVC_preBD_L_Baseline)%>%
  mutate(Tiffeneau_imputated=FEV1_preBD_L_Baseline_imputated/FVC_preBD_L_Baseline_imputated)%>%
  mutate(Tiffeneau_imputated_no_systematically_missing=FEV1_preBD_L_Baseline_imputated_no_systematically_missing/FVC_preBD_L_Baseline_imputated_no_systematically_missing)

## Calculate the value in absolute for eosinophils, FeNO and IgE
Data_Oracle$IgE_Log_imputated
Data_Oracle<-Data_Oracle %>%
  mutate(Blood_Eos_baseline_x10_9_cells_per_L_imputated=10^Eosinophils_Log_imputated) %>%
  mutate(Blood_Eos_baseline_x10_9_cells_per_L_imputated_no_systematically_missing=10^Eosinophils_Log_imputated_no_systematically_missing) %>%
  
  mutate(FeNO_baseline_ppb_imputated=10^FeNO_Log_imputated)%>%
  mutate(FeNO_baseline_ppb_imputated_no_systematically_missing=10^FeNO_Log_imputated_no_systematically_missing)%>%
  
  mutate(Total_IgE_imputated=10^IgE_Log_imputated) %>%
  mutate(Total_IgE_imputated_no_systematically_missing=10^IgE_Log_imputated_no_systematically_missing)

##Put the Follow up duration in days
Data_Oracle$Follow_up_duration_days_nozero
Data_Oracle<-Data_Oracle %>%
  mutate(Follow_up_duration_days=Data_Oracle$Follow_up_duration_days_nozero)
Data_Oracle$Follow_up_duration_days
#======================================================================================================================================
## Calculation of the parameters by a definite change
Data_Oracle<-Data_Oracle%>%
  #Age per 10 year increase
  mutate(Age_per_10=Age/10) %>%
  mutate(Age_per_10_imputated=Age_imputated/10) %>%
  mutate(Age_per_10_imputated_no_systematically_missing=Age_imputated_no_systematically_missing/10) %>%
  #BMI per 5 increase
  mutate(BMI_per_5=BMI/5) %>%
  mutate(BMI_per_5_imputated=BMI_imputated/5) %>%
  mutate(BMI_per_5_imputated_no_systematically_missing=BMI_imputated_no_systematically_missing/5) %>%
  #FEV1 per 10% decrease
  mutate(FEV1_preBD_per10_Baseline=-FEV1_preBD_PCT_Baseline/10) %>%
  mutate(FEV1_preBD_per10_Baseline_imputated=-FEV1_preBD_PCT_Baseline_imputated/10) %>%
  mutate(FEV1_preBD_per10_Baseline_imputated_no_systematically_missing=-FEV1_preBD_PCT_Baseline_imputated_no_systematically_missing/10) %>%
  #Reversibility per 10%
  mutate(FEV1_per10_reversibilityBD= FEV1_PCT_reversibility_postBD/10) %>%
  mutate(FEV1_per10_reversibilityBD_imputated= FEV1_PCT_reversibility_postBD_imputated/10) %>%
  mutate(FEV1_per10_reversibilityBD_imputated_no_systematically_missing= FEV1_PCT_reversibility_postBD_imputated_no_systematically_missing/10)
#======================================================================================================================================
#Creation of new categoricals data
##Create the CRsNP factor (CRS without NP)
#CREATING of CRSwNP

Data_Oracle <- Data_Oracle  %>%
  mutate(CRSwNP=case_when(
    #Identifying as "Yes" for patients with NP
    Nasal_polyposis=="Yes"~ "Yes",
    #Identifying as "Yes" for patients with polypectomy history
    Previous_nasal_polypectomy_0no_1yes_9999notknown =="Yes"~ "Yes",
    #Identifying as "No" for patients without NP
    Nasal_polyposis=="No" ~ "No",
    #Identifying as NA for patients that NP was not assessed
    is.na(Nasal_polyposis) ~ NA,
    TRUE ~ NA
  )) %>%
  mutate(CRSwNP_imputated=case_when(
    #Identifying as "Yes" for patients with NP
    Nasal_polyposis_imputated=="Yes"~ "Yes",
    #Identifying as "Yes" for patients with polypectomy history
    Previous_nasal_polypectomy_0no_1yes_9999notknown_imputated =="Yes"~ "Yes",
    #Identifying as "No" for patients without NP
    Nasal_polyposis_imputated=="No" ~ "No",
    #Identifying as NA for patients that NP was not assessed
    is.na(Nasal_polyposis_imputated) ~ NA,
    TRUE ~ NA
  )) %>%
  mutate(CRSwNP_imputated_no_systematically_missing=case_when(
    #Identifying as "Yes" for patients with NP
    Nasal_polyposis_imputated_no_systematically_missing=="Yes"~ "Yes",
    #Identifying as "Yes" for patients with polypectomy history
    Previous_nasal_polypectomy_0no_1yes_9999notknown_imputated_no_systematically_missing =="Yes"~ "Yes",
    #Identifying as "No" for patients without NP
    Nasal_polyposis_imputated_no_systematically_missing=="No" ~ "No",
    #Identifying as NA for patients that NP was not assessed
    is.na(Nasal_polyposis_imputated_no_systematically_missing) ~ NA,
    TRUE ~ NA
  ))

Data_Oracle$CRSwNP <- factor(Data_Oracle$CRSwNP, levels=c("No","Yes"), labels=c("No","Yes"))
Data_Oracle$CRSwNP_imputated <- factor(Data_Oracle$CRSwNP_imputated, levels=c("No","Yes"), labels=c("No","Yes"))
Data_Oracle$CRSwNP_imputated_no_systematically_missing <- factor(Data_Oracle$CRSwNP_imputated_no_systematically_missing, levels=c("No","Yes"), labels=c("No","Yes"))
##CRSsNP (No CRSwNP + No history of nasal polypectomy + CRS=="Yes")
Data_Oracle<- Data_Oracle %>%
  mutate(CRSsNP=case_when(
    #Identifying as NA for patients that nasal polyposis was not assessed
    is.na(CRSwNP)~ NA,
    #Identifying as NA for patients that Chronic_Rhinosinusitis was not assessed
    is.na(Chronic_rhinosinusitis) ~ NA,
    #Identifying as "Yes" for patients with CRSwNP
    CRSwNP=="Yes"~ "No",
    #Identifying as "No" for patients without Chronic rhinosinusitis
    Chronic_rhinosinusitis=="No" ~ "No",
    #Identifying as "Yes" for patients with Chronic rhinosinusitis
    Chronic_rhinosinusitis=="Yes" ~ "Yes",
    TRUE ~ NA
  ))%>%
  mutate(CRSsNP_imputated=case_when(
    #Identifying as NA for patients that nasal polyposis was not assessed
    is.na(CRSwNP_imputated)~ NA,
    #Identifying as NA for patients that Chronic_Rhinosinusitis was not assessed
    is.na(Chronic_rhinosinusitis_imputated) ~ NA,
    #Identifying as "Yes" for patients with CRSwNP
    CRSwNP_imputated=="Yes"~ "No",
    #Identifying as "No" for patients without Chronic rhinosinusitis
    Chronic_rhinosinusitis_imputated=="No" ~ "No",
    #Identifying as "Yes" for patients with Chronic rhinosinusitis
    Chronic_rhinosinusitis_imputated=="Yes" ~ "Yes",
    TRUE ~ NA
  ))%>%
  mutate(CRSsNP_imputated_no_systematically_missing=case_when(
    #Identifying as NA for patients that nasal polyposis was not assessed
    is.na(CRSwNP_imputated_no_systematically_missing)~ NA,
    #Identifying as NA for patients that Chronic_Rhinosinusitis was not assessed
    is.na(Chronic_rhinosinusitis_imputated_no_systematically_missing) ~ NA,
    #Identifying as "Yes" for patients with CRSwNP
    CRSwNP_imputated_no_systematically_missing=="Yes"~ "No",
    #Identifying as "No" for patients without Chronic rhinosinusitis
    Chronic_rhinosinusitis_imputated_no_systematically_missing=="No" ~ "No",
    #Identifying as "Yes" for patients with Chronic rhinosinusitis
    Chronic_rhinosinusitis_imputated_no_systematically_missing=="Yes" ~ "Yes",
    TRUE ~ NA
  ))
Data_Oracle$CRSsNP <- factor(Data_Oracle$CRSsNP, levels=c("No","Yes"), labels=c("No","Yes"))
Data_Oracle$CRSsNP_imputated <- factor(Data_Oracle$CRSsNP_imputated, levels=c("No","Yes"), labels=c("No","Yes"))
Data_Oracle$CRSsNP_imputated_no_systematically_missing <- factor(Data_Oracle$CRSsNP_imputated_no_systematically_missing, levels=c("No","Yes"), labels=c("No","Yes"))

##Category by inflammatory marker
Data_Oracle$Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced
Data_Oracle$Eosinophils_by_group<-cut(Data_Oracle$Blood_Eos_baseline_x10_9_cells_per_L_zeroreplaced,breaks = c(0,0.15, 0.3,100000),labels=c('<0.15', '0.15-0.3', '>0.3'))
Data_Oracle$Eosinophils_by_group_imputated<-cut(Data_Oracle$Blood_Eos_baseline_x10_9_cells_per_L_imputated,breaks = c(0,0.15, 0.3,100000),labels=c('<0.15', '0.15-0.3', '>0.3'))
Data_Oracle$Eosinophils_by_group_imputated_no_systematically_missing<-cut(Data_Oracle$Blood_Eos_baseline_x10_9_cells_per_L_imputated_no_systematically_missing,breaks = c(0,0.15, 0.3,100000),labels=c('<0.15', '0.15-0.3', '>0.3'))

Data_Oracle$FeNO_baseline_by_group<-cut(Data_Oracle$FeNO_baseline_ppb,breaks = c(0,25, 50,100000),labels=c('<25', '25-50', '>50'))
Data_Oracle$FeNO_baseline_by_group_imputated<-cut(Data_Oracle$FeNO_baseline_ppb_imputated,breaks = c(0,25, 50,100000),labels=c('<25', '25-50', '>50'))
Data_Oracle$FeNO_baseline_by_group_imputated_no_systematically_missing<-cut(Data_Oracle$FeNO_baseline_ppb_imputated_no_systematically_missing,breaks = c(0,25, 50,100000),labels=c('<25', '25-50', '>50'))


Data_Oracle$IgE_by_group<-cut(Data_Oracle$Total_IgE,breaks = c(0,150, 600,100000),labels=c('<150', '150-600', '>600'))
Data_Oracle$IgE_by_group_imputated<-cut(Data_Oracle$Total_IgE_imputated,breaks = c(0,150, 600,100000),labels=c('<150', '150-600', '>600'))
Data_Oracle$IgE_by_group_imputated_no_systematically_missing<-cut(Data_Oracle$Total_IgE_imputated_no_systematically_missing,breaks = c(0,150, 600,100000),labels=c('<150', '150-600', '>600'))


##Category by lung function
Data_Oracle$FEV1_preBD_Baseline_by_group<-cut(Data_Oracle$FEV1_preBD_PCT_Baseline,breaks = c(0,50,60,70,100000),labels=c('<50%',"50-60%",'60-70%',">70%"))
Data_Oracle$FEV1_preBD_Baseline_by_group_imputated<-cut(Data_Oracle$FEV1_preBD_PCT_Baseline_imputated,breaks = c(0,50,60,70,100000),labels=c('<50%',"50-60%",'60-70%',">70%"))
Data_Oracle$FEV1_preBD_Baseline_by_group_imputated_no_systematically_missing<-cut(Data_Oracle$FEV1_preBD_PCT_Baseline_imputated_no_systematically_missing,breaks = c(0,50,60,70,100000),labels=c('<50%',"50-60%",'60-70%',">70%"))

Data_Oracle$FEV1_preBD_Baseline_by_group <- factor(Data_Oracle$FEV1_preBD_Baseline_by_group, levels = c(">70%","60-70%","50-60%","<50%"))
Data_Oracle$FEV1_preBD_Baseline_by_group_imputated <- factor(Data_Oracle$FEV1_preBD_Baseline_by_group_imputated, levels = c(">70%","60-70%","50-60%","<50%"))
Data_Oracle$FEV1_preBD_Baseline_by_group_imputated_no_systematically_missing <- factor(Data_Oracle$FEV1_preBD_Baseline_by_group_imputated_no_systematically_missing, levels = c(">70%","60-70%","50-60%","<50%"))


##Category by GINA_treatment_step
Data_Oracle<-Data_Oracle %>%
  mutate(Treatment_step= case_when(Treatment_step=="1"~"Step 1",Treatment_step=="2"~"Step 2",Treatment_step=="3"~"Step 3",Treatment_step=="4"~"Step 4",Treatment_step=="5"~"Step 5",TRUE~Treatment_step))
Data_Oracle$Treatment_step
Data_Oracle$Treatment_step <- factor(Data_Oracle$Treatment_step, levels = c("Step 1","Step 2","Step 3","Step 4","Step 5"))

Data_Oracle<-Data_Oracle %>%
  mutate(Treatment_step_combine= case_when(Treatment_step=="Step 1"~"Step 1-2",Treatment_step=="Step 2"~"Step 1-2",TRUE~Treatment_step))
Data_Oracle$Treatment_step_combine <- factor(Data_Oracle$Treatment_step_combine, levels = c("Step 3","Step 1-2","Step 4","Step 5"))

##Create a treatment step_1-2_vs_3-4-5
Data_Oracle <- Data_Oracle %>%
  mutate(Treatment_step_1_2vs3_5 = factor(case_when(
    Treatment_step %in% c("Step 1", "Step 2") ~ "Step 1-2",
    Treatment_step %in% c("Step 3", "Step 4", "Step 5") ~ "Step 3-5",
    TRUE ~ Treatment_step  # Keep other values as they are
  )))

##Create a treatment step_1-2_vs_3-4_vs_5
Data_Oracle <- Data_Oracle %>%
  mutate(Treatment_step_1_2vs3_4vs5 = factor(case_when(
    Treatment_step %in% c("Step 1", "Step 2") ~ "Step 1-2",
    Treatment_step %in% c("Step 3", "Step 4") ~ "Step 3-4",
    Treatment_step %in% c("Step 5") ~ "Step 5",
    TRUE ~ Treatment_step  # Keep other values as they are
  )))

Data_Oracle <- Data_Oracle %>%
  mutate(Treatment_step_1_2vs3_4vs5 = relevel(Treatment_step_1_2vs3_4vs5, ref = "Step 3-4"))

##Category by ACQ-5
Data_Oracle$ACQ5_by_group<-cut(Data_Oracle$ACQ_baseline_score_mean,breaks = c(-10,1.5,3,100000),labels=c("<1.5","1.5-3",">3"))
Data_Oracle$ACQ5_by_group_imputated<-cut(Data_Oracle$ACQ_baseline_score_mean_imputated,breaks = c(-10,1.5,3,100000),labels=c("<1.5","1.5-3",">3"))
Data_Oracle$ACQ5_by_group_imputated_no_systematically_missing<-cut(Data_Oracle$ACQ_baseline_score_mean_imputated_no_systematically_missing,breaks = c(-10,1.5,3,100000),labels=c("<1.5","1.5-3",">3"))

##Category by BMI
Data_Oracle$BMI_by_group<-cut(Data_Oracle$BMI,breaks = c(-10,25,30,35,100000),labels=c("<25","25-30","30-35",">35"))
Data_Oracle$BMI_by_group_imputated<-cut(Data_Oracle$BMI_imputated,breaks = c(-10,25,30,35,100000),labels=c("<25","25-30","30-35",">35"))
Data_Oracle$BMI_by_group_imputated_no_systematically_missing<-cut(Data_Oracle$BMI_imputated_no_systematically_missing,breaks = c(-10,25,30,35,100000),labels=c("<25","25-30","30-35",">35"))

##Category by Age group
Data_Oracle$Age_by_group_imputated<-cut(Data_Oracle$Age_imputated,breaks = c(-10,40,50,60,100000),labels=c("<40","40-50","50-60",">60"))

#Category of ICS DOSE_CLASS
Data_Oracle$ICS_DOSE_CLASS <- factor(Data_Oracle$ICS_DOSE_CLASS,
                                     levels = c("0", "Low", "Medium", "High"))
Data_Oracle$ICS_DOSE_CLASS <- relevel(Data_Oracle$ICS_DOSE_CLASS, ref = "High")
Data_Oracle$ICS_DOSE_NUMERIC <- as.numeric(factor(Data_Oracle$ICS_DOSE_CLASS,
                                                  levels = c("0", "Low", "Medium", "High"))) - 1
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Extrapolation of the missing value of lung function

##For the FVC predicted
Data_Oracle<- Data_Oracle %>%
  mutate(FVC_predicted_L = case_when(
    is.na(FVC_predicted_L)&!is.na(FVC_preBD_L_Baseline)&!is.na(FVC_preBD_PCT_Baseline) ~ FVC_preBD_L_Baseline*100/FVC_preBD_PCT_Baseline,
    is.na(FVC_predicted_L)&!is.na(FVC_postBD_L_Baseline)&!is.na(FVC_postBD_PCT_Baseline) ~ FVC_postBD_L_Baseline*100/FVC_postBD_PCT_Baseline,
    TRUE ~ FVC_predicted_L))


##For the FEV1/FVC preBD
Data_Oracle<- Data_Oracle %>%
  mutate(FEV1_FVC_ratio = FEV1_preBD_L_Baseline/FVC_preBD_L_Baseline)

##For the FEV1/FVC postBD
Data_Oracle<- Data_Oracle %>%
  mutate(FEV1_FVC_ratio_postBD = FEV1_postBD_L_Baseline/FVC_postBD_L_Baseline)

##For the FEV1 post BD in L
Data_Oracle<- Data_Oracle %>%
  mutate(FEV1_postBD_L_Baseline = case_when(
    is.na(FEV1_postBD_L_Baseline)&!is.na(FEV1_postBD_PCT_Baseline)&!is.na(FEV1_predicted_L) ~ (FEV1_postBD_PCT_Baseline/100)*FEV1_predicted_L,
    TRUE ~ FEV1_postBD_L_Baseline))

##For the reversibility
Data_Oracle<- Data_Oracle %>%
  mutate(FEV1_PCT_reversibility_postBD = case_when(
    is.na(FEV1_PCT_reversibility_postBD)&!is.na(FEV1_preBD_L_Baseline)&!is.na(FEV1_postBD_L_Baseline) ~ ((FEV1_postBD_L_Baseline-FEV1_preBD_L_Baseline)/FEV1_preBD_L_Baseline)*100,
    TRUE ~ FEV1_PCT_reversibility_postBD))

##For the FEV1 post BD at 52week in %
Data_Oracle<- Data_Oracle %>%
  mutate(FEV1POSTBD_PCT_52W = case_when(
    is.na(FEV1POSTBD_PCT_52W)&!is.na(FEV1POSTBD_L_52W)&!is.na(FEV1_predicted_L) ~ (FEV1POSTBD_L_52W/FEV1_predicted_L)*100,
    TRUE ~ FEV1POSTBD_PCT_52W))

##Calculation of delta_FEV1_preBD
Data_Oracle<- Data_Oracle %>%
  mutate(delta_FEV1_preBD_L= FEV1PREBD_L_52W-FEV1_preBD_L_Baseline) %>%
  mutate(delta_FEV1_preBD_mL= (FEV1PREBD_L_52W-FEV1_preBD_L_Baseline)*1000) %>%
  mutate(delta_FEV1_preBD_PCT= FEV1PREBD_PCT_52W-FEV1_preBD_PCT_Baseline)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Data preparation for country
Data_Oracle$Country <- factor(Data_Oracle$Country)
levels(Data_Oracle$Country)

Data_Oracle <- Data_Oracle %>%
  mutate(Country = case_when(
    Country %in% c("United_States", "United States") ~ "United States",
    Country %in% c("Hungary", "Hungria") ~ "Hungary",
    Country %in% c("Russia", "Russian Federation") ~ "Russia",
    Country %in% c("South_Korea", "Korea","Korea, Republic of") ~ "South Korea",
    Country %in% c("South_Africa", "South Africa") ~ "South Africa",
    Country %in% c("UK - CMD", "United Kingdom") ~ "United Kingdom",
    TRUE ~ Country
  ))
Data_Oracle$Country <- factor(Data_Oracle$Country)
Data_Oracle$Country <- relevel(Data_Oracle$Country, ref = "United States")

Data_Oracle <- Data_Oracle %>%
  mutate(Country_per_region_1 = case_when(
    Country %in% c("United States", "Canada", "Mexico", "North America") ~ "North America",
    Country %in% c("Argentina", "Chile", "Peru", "Brazil", "South America") ~ "South America",
    Country %in% c("Australia", "New Zealand", "Oceania") ~ "Oceania",
    Country == "South Africa" ~ "Africa South",
    Country %in% c("Italy", "United Kingdom", "Belgium", "Poland", "Romania", "Hungary",
                   "Czech Republic", "Belarus", "Ukraine", "Spain", "Slovakia", "Russia",
                   "France", "Serbia", "Bulgaria", "Germany", "Netherlands", "Latvia",
                   "Lithuania", "Europe") ~ "Europe",
    Country %in% c("Japan", "Korea", "Israel", "Turkey", "Vietnam", "South Korea", "Asia") ~ "Asia",
    TRUE ~ NA
  ))

Data_Oracle <- Data_Oracle %>%
  mutate(Country_per_region_2 = case_when(
    Country %in% c("United States", "Canada", "Mexico") ~ "North America",
    Country %in% c("Argentina", "Chile", "Peru", "Brazil") ~ "South America",
    Country %in% c("Australia", "New_Zealand") ~ "Oceania",
    Country == "South_Africa" ~ "Africa_South",
    Country %in% c("Israel", "Turkey") ~ "Middle East",
    Country %in% c("Japan", "South Korea", "Vietnam") ~ "Asia",
    Country %in% c("United Kingdom", "France", "Germany", "Netherlands") ~ "Europe_Western",
    Country %in% c("Spain", "Italy") ~ "Europe_Southern",
    Country %in% c("Poland", "Hungary", "Czech_Republic", "Slovakia") ~ "Europe_Central Eastern",
    Country %in% c("Belarus", "Ukraine", "Russia", "Lithuania", "Latvia") ~ "Europe_Eastern",
    Country %in% c("Romania", "Bulgaria", "Serbia") ~ "Europe_South Eastern",
    TRUE ~ NA
  ))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#PREPARATION OF DATA FOR THE DIFFERENT ANALYSIS
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Placebo response on asthma attack
#Selection of the studies with attack history in numeric
Studies_included_attack<- c("AZISAST","BENRAP2B","DREAM","DRI12544","EXTRA", "LAVOLTA_1","LAVOLTA_2","LUSTER_1","LUSTER_2","NAVIGATOR","PATHWAY","QUEST","STRATOS_1","STRATOS_2")

Data_Oracle_attack<-Data_Oracle %>%
  filter(Enrolled_Trial_name%in%Studies_included_attack)
nrow(Data_Oracle_attack%>%
       filter(.imp==1))
#======================================================================================================================================
##Preparation for delta asthma attack analysis

colnames(Data_Oracle_attack$Attack_12mo_Nb_imputated)

Data_Oracle_attack<- Data_Oracle_attack %>%
  mutate(Delta_attack= Attack_12mo_Nb_imputated-Number_severe_asthma_attacks_during_followup)%>%
  mutate(Delta_attack_notimputated= Attack_12mo_Nb-Number_severe_asthma_attacks_during_followup) 

Data_Oracle_attack$Treatment_step_ordinal <- as.numeric(gsub("Step ", "", as.character(Data_Oracle_attack$Treatment_step)))
summary(Data_Oracle_attack$Treatment_step_ordinal)

##Calculate Rate Attack in trial
Data_Oracle_attack <- Data_Oracle_attack %>%
  mutate(
    Attack_rate_before_trial = Attack_12mo_Nb_imputated / 1,
    Attack_rate_during_trial = Number_severe_asthma_attacks_during_followup / Follow_up_duration_years,
    Delta_Attack_rate = Attack_rate_before_trial - Attack_rate_during_trial,
    Attack_rate_before_trial_log = log(Attack_rate_before_trial + 1)
  )

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## For lung function change
#Selection of the studies with the total ∆FEV1
Studies_included_lung<- c("BENRAP2B","LAVOLTA_1","LAVOLTA_2","NAVIGATOR","PATHWAY","QUEST","STRATOS_1","STRATOS_2")
Data_Oracle_lung<-Data_Oracle %>%
  filter(Enrolled_Trial_name%in%Studies_included_lung)
Data_Oracle_lung$delta_FEV1_preBD_mL
colnames(Data_Oracle_lung)

#Put Treatment step in ordinal
Data_Oracle_lung$Treatment_step_ordinal <- as.numeric(gsub("Step ", "", as.character(Data_Oracle_lung$Treatment_step)))
summary(Data_Oracle_lung$Treatment_step_ordinal)

nrow(Data_Oracle_lung%>%
       filter(.imp==1))

#Selection of patients with delta FEV1 is not NA
Data_Oracle_lung<-Data_Oracle_lung %>%
  filter(!is.na(delta_FEV1_preBD_mL))

#======================================================================================================================================
#======================================================================================================================================
#Preparation of the data with change over time

# Rename Day_*w to Lung_Day_*W
colnames(Raw_ACQ_FEV1) <- gsub("^Day_(\\d+)w$", "Lung_day_\\1W", colnames(Raw_ACQ_FEV1))

## Exclude MILLY from the analysis
#Raw_ACQ_FEV1<- Raw_ACQ_FEV1 %>%
#filter(Enrolled_Trial_name != "MILLY")

# Identify the patient(s) with FEV1PREBD_PCT_24W == 249
patients_to_update <- Raw_ACQ_FEV1 %>%
  filter(FEV1PREBD_PCT_24W == 249) %>%
  pull(Subject_ID)  # Extract Subject_ID(s)

# Replace the selected variables with NA for these patients
Raw_ACQ_FEV1 <- Raw_ACQ_FEV1 %>%
  mutate(
    FEV1PREBD_L_24W = if_else(Subject_ID %in% patients_to_update, NA_real_, FEV1PREBD_L_24W),
    FEV1PREBD_PCT_24W = if_else(Subject_ID %in% patients_to_update, NA_real_, FEV1PREBD_PCT_24W),
    FVCPREBD_L_24W = if_else(Subject_ID %in% patients_to_update, NA_real_, FVCPREBD_L_24W),
    FVCPREBD_PCT_24W = if_else(Subject_ID %in% patients_to_update, NA_real_, FVCPREBD_PCT_24W),
    FEV1POSTBD_L_24W = if_else(Subject_ID %in% patients_to_update, NA_real_, FEV1POSTBD_L_24W),
    Lung_day_24W= if_else(Subject_ID %in% patients_to_update, NA_real_, Lung_day_24W)
  )

# Check if the values were updated correctly
detach("package:Matrix", unload = TRUE)
Raw_ACQ_FEV1 %>%
  as.data.frame() %>%
  dplyr::filter(Subject_ID %in% patients_to_update) %>%
  dplyr::select(
    Subject_ID,
    FEV1PREBD_L_24W,
    FEV1PREBD_PCT_24W,
    FVCPREBD_L_24W,
    FVCPREBD_PCT_24W,
    FEV1POSTBD_L_24W
  )

#Calculation of predicted value for each subject
Raw_ACQ_FEV1<- Raw_ACQ_FEV1 %>%
  mutate(FEV1_predicted_1W= FEV1PREBD_L_1W/(FEV1PREBD_PCT_1W/100)) %>%
  mutate(FEV1_predicted_2W= FEV1PREBD_L_2W/(FEV1PREBD_PCT_2W/100))%>%
  mutate(FEV1_predicted_4W= FEV1PREBD_L_4W/(FEV1PREBD_PCT_4W/100))%>%
  mutate(FEV1_predicted_8W= FEV1PREBD_L_8W/(FEV1PREBD_PCT_8W/100))%>%
  mutate(FEV1_predicted_12W= FEV1PREBD_L_12W/(FEV1PREBD_PCT_12W/100))%>%
  mutate(FEV1_predicted_16W= FEV1PREBD_L_16W/(FEV1PREBD_PCT_16W/100))%>%
  mutate(FEV1_predicted_20W= FEV1PREBD_L_20W/(FEV1PREBD_PCT_20W/100))%>%
  mutate(FEV1_predicted_21W= FEV1PREBD_L_21W/(FEV1PREBD_PCT_21W/100))%>%
  mutate(FEV1_predicted_24W= FEV1PREBD_L_24W/(FEV1PREBD_PCT_24W/100))%>%
  mutate(FEV1_predicted_26W= FEV1PREBD_L_26W/(FEV1PREBD_PCT_26W/100))%>%
  mutate(FEV1_predicted_28W= FEV1PREBD_L_28W/(FEV1PREBD_PCT_28W/100))%>%
  mutate(FEV1_predicted_32W= FEV1PREBD_L_32W/(FEV1PREBD_PCT_32W/100))%>%
  mutate(FEV1_predicted_36W= FEV1PREBD_L_36W/(FEV1PREBD_PCT_36W/100))%>%
  mutate(FEV1_predicted_40W= FEV1PREBD_L_40W/(FEV1PREBD_PCT_40W/100))%>%
  mutate(FEV1_predicted_44W= FEV1PREBD_L_44W/(FEV1PREBD_PCT_44W/100))%>%
  mutate(FEV1_predicted_48W= FEV1PREBD_L_48W/(FEV1PREBD_PCT_48W/100))%>%
  mutate(FEV1_predicted_52W= FEV1PREBD_L_52W/(FEV1PREBD_PCT_52W/100)) %>%
  
  mutate(FVC_predicted_1W= FVCPREBD_L_1W/(FVCPREBD_PCT_1W/100)) %>%
  mutate(FVC_predicted_2W= FVCPREBD_L_2W/(FVCPREBD_PCT_2W/100))%>%
  mutate(FVC_predicted_4W= FVCPREBD_L_4W/(FVCPREBD_PCT_4W/100))%>%
  mutate(FVC_predicted_8W= FVCPREBD_L_8W/(FVCPREBD_PCT_8W/100))%>%
  mutate(FVC_predicted_12W= FVCPREBD_L_12W/(FVCPREBD_PCT_12W/100))%>%
  mutate(FVC_predicted_16W= FVCPREBD_L_16W/(FVCPREBD_PCT_16W/100))%>%
  mutate(FVC_predicted_20W= FVCPREBD_L_20W/(FVCPREBD_PCT_20W/100))%>%
  mutate(FVC_predicted_21W= FVCPREBD_L_21W/(FVCPREBD_PCT_21W/100))%>%
  mutate(FVC_predicted_24W= FVCPREBD_L_24W/(FVCPREBD_PCT_24W/100))%>%
  mutate(FVC_predicted_26W= FVCPREBD_L_26W/(FVCPREBD_PCT_26W/100))%>%
  mutate(FVC_predicted_28W= FVCPREBD_L_28W/(FVCPREBD_PCT_28W/100))%>%
  mutate(FVC_predicted_32W= FVCPREBD_L_32W/(FVCPREBD_PCT_32W/100))%>%
  mutate(FVC_predicted_36W= FVCPREBD_L_36W/(FVCPREBD_PCT_36W/100))%>%
  mutate(FVC_predicted_40W= FVCPREBD_L_40W/(FVCPREBD_PCT_40W/100))%>%
  mutate(FVC_predicted_44W= FVCPREBD_L_44W/(FVCPREBD_PCT_44W/100))%>%
  mutate(FVC_predicted_48W= FVCPREBD_L_48W/(FVCPREBD_PCT_48W/100))%>%
  mutate(FVC_predicted_52W= FVCPREBD_L_52W/(FVCPREBD_PCT_52W/100)) %>%
  
  mutate(Tif_predicted_1W= FEV1_predicted_1W/FVC_predicted_1W) %>%
  mutate(Tif_predicted_2W= FEV1_predicted_2W/FVC_predicted_2W) %>%
  mutate(Tif_predicted_4W= FEV1_predicted_4W/FVC_predicted_4W) %>%
  mutate(Tif_predicted_8W= FEV1_predicted_8W/FVC_predicted_8W) %>%
  mutate(Tif_predicted_12W= FEV1_predicted_12W/FVC_predicted_12W) %>%
  mutate(Tif_predicted_16W=FEV1_predicted_16W/FVC_predicted_16W) %>%
  mutate(Tif_predicted_20W= FEV1_predicted_20W/FVC_predicted_20W) %>%
  mutate(Tif_predicted_21W= FEV1_predicted_21W/FVC_predicted_21W) %>%
  mutate(Tif_predicted_24W= FEV1_predicted_24W/FVC_predicted_24W) %>%
  mutate(Tif_predicted_26W= FEV1_predicted_26W/FVC_predicted_26W) %>%
  mutate(Tif_predicted_28W= FEV1_predicted_28W/FVC_predicted_28W) %>%
  mutate(Tif_predicted_32W= FEV1_predicted_32W/FVC_predicted_32W) %>%
  mutate(Tif_predicted_36W= FEV1_predicted_36W/FVC_predicted_36W) %>%
  mutate(Tif_predicted_40W= FEV1_predicted_40W/FVC_predicted_40W) %>%
  mutate(Tif_predicted_44W= FEV1_predicted_44W/FVC_predicted_44W) %>%
  mutate(Tif_predicted_48W= FEV1_predicted_48W/FVC_predicted_48W) %>%
  mutate(Tif_predicted_52W= FEV1_predicted_52W/FVC_predicted_52W)


#Calculation of Tif for each time
Raw_ACQ_FEV1<- Raw_ACQ_FEV1 %>%
  mutate(TifPREBD_1W= FEV1PREBD_L_1W/FVCPREBD_L_1W) %>%
  mutate(TifPREBD_2W= FEV1PREBD_L_2W/FVCPREBD_L_2W) %>%
  mutate(TifPREBD_4W= FEV1PREBD_L_4W/FVCPREBD_L_4W) %>%
  mutate(TifPREBD_8W= FEV1PREBD_L_8W/FVCPREBD_L_8W) %>%
  mutate(TifPREBD_12W= FEV1PREBD_L_12W/FVCPREBD_L_12W) %>%
  mutate(TifPREBD_16W= FEV1PREBD_L_16W/FVCPREBD_L_16W) %>%
  mutate(TifPREBD_20W= FEV1PREBD_L_20W/FVCPREBD_L_20W) %>%
  mutate(TifPREBD_21W= FEV1PREBD_L_21W/FVCPREBD_L_21W) %>%
  mutate(TifPREBD_24W= FEV1PREBD_L_24W/FVCPREBD_L_24W) %>%
  mutate(TifPREBD_26W= FEV1PREBD_L_26W/FVCPREBD_L_26W) %>%
  mutate(TifPREBD_28W= FEV1PREBD_L_28W/FVCPREBD_L_28W) %>%
  mutate(TifPREBD_32W= FEV1PREBD_L_32W/FVCPREBD_L_32W) %>%
  mutate(TifPREBD_36W= FEV1PREBD_L_36W/FVCPREBD_L_36W) %>%
  mutate(TifPREBD_40W= FEV1PREBD_L_40W/FVCPREBD_L_40W) %>%
  mutate(TifPREBD_44W= FEV1PREBD_L_44W/FVCPREBD_L_44W) %>%
  mutate(TifPREBD_48W= FEV1PREBD_L_48W/FVCPREBD_L_48W) %>%
  mutate(TifPREBD_52W= FEV1PREBD_L_52W/FVCPREBD_L_52W) %>%
  
  
  mutate(TifPOSTBD_1W= TifPREBD_1W/FVCPOSTBD_L_1W) %>%
  mutate(TifPOSTBD_2W= FEV1POSTBD_L_2W/FVCPOSTBD_L_2W) %>%
  mutate(TifPOSTBD_4W= FEV1POSTBD_L_4W/FVCPOSTBD_L_4W) %>%
  mutate(TifPOSTBD_8W= FEV1POSTBD_L_8W/FVCPOSTBD_L_8W) %>%
  mutate(TifPOSTBD_12W= FEV1POSTBD_L_12W/FVCPOSTBD_L_12W) %>%
  mutate(TifPOSTBD_16W= FEV1POSTBD_L_16W/FVCPOSTBD_L_16W) %>%
  mutate(TifPOSTBD_20W= FEV1POSTBD_L_20W/FVCPOSTBD_L_20W) %>%
  mutate(TifPOSTBD_21W= FEV1POSTBD_L_21W/FVCPOSTBD_L_21W) %>%
  mutate(TifPOSTBD_24W= FEV1POSTBD_L_24W/FVCPOSTBD_L_24W) %>%
  mutate(TifPOSTBD_26W= FEV1POSTBD_L_26W/FVCPOSTBD_L_26W) %>%
  mutate(TifPOSTBD_28W= FEV1POSTBD_L_28W/FVCPOSTBD_L_28W) %>%
  mutate(TifPOSTBD_32W= FEV1POSTBD_L_32W/FVCPOSTBD_L_32W) %>%
  mutate(TifPOSTBD_36W= FEV1POSTBD_L_36W/FVCPOSTBD_L_36W) %>%
  mutate(TifPOSTBD_40W= FEV1POSTBD_L_40W/FVCPOSTBD_L_40W) %>%
  mutate(TifPOSTBD_44W= FEV1POSTBD_L_44W/FVCPOSTBD_L_44W) %>%
  mutate(TifPOSTBD_48W= FEV1POSTBD_L_48W/FVCPOSTBD_L_48W) %>%
  mutate(TifPOSTBD_52W= FEV1POSTBD_L_52W/FVCPOSTBD_L_52W)

#Calculation of TifPCT for each time
Raw_ACQ_FEV1<- Raw_ACQ_FEV1 %>%
  mutate(TifPREBD_PCT_1W= (TifPREBD_1W/Tif_predicted_1W)*100) %>%
  mutate(TifPREBD_PCT_2W= (TifPREBD_2W/Tif_predicted_2W)*100) %>%
  mutate(TifPREBD_PCT_4W= (TifPREBD_4W/Tif_predicted_4W)*100) %>%
  mutate(TifPREBD_PCT_8W= (TifPREBD_8W/Tif_predicted_8W)*100) %>%
  mutate(TifPREBD_PCT_12W= (TifPREBD_12W/Tif_predicted_12W)*100) %>%
  mutate(TifPREBD_PCT_16W= (TifPREBD_16W/Tif_predicted_16W)*100) %>%
  mutate(TifPREBD_PCT_20W= (TifPREBD_20W/Tif_predicted_20W)*100) %>%
  mutate(TifPREBD_PCT_21W= (TifPREBD_21W/Tif_predicted_21W)*100) %>%
  mutate(TifPREBD_PCT_24W= (TifPREBD_24W/Tif_predicted_24W)*100) %>%
  mutate(TifPREBD_PCT_26W= (TifPREBD_26W/Tif_predicted_26W)*100) %>%
  mutate(TifPREBD_PCT_28W= (TifPREBD_28W/Tif_predicted_28W)*100) %>%
  mutate(TifPREBD_PCT_32W= (TifPREBD_32W/Tif_predicted_32W)*100) %>%
  mutate(TifPREBD_PCT_36W= (TifPREBD_36W/Tif_predicted_36W)*100) %>%
  mutate(TifPREBD_PCT_40W= (TifPREBD_40W/Tif_predicted_40W)*100) %>%
  mutate(TifPREBD_PCT_44W= (TifPREBD_44W/Tif_predicted_44W)*100) %>%
  mutate(TifPREBD_PCT_48W= (TifPREBD_48W/Tif_predicted_48W)*100) %>%
  mutate(TifPREBD_PCT_52W= (TifPREBD_52W/Tif_predicted_52W)*100) %>%
  
  mutate(TifPOSTBD_PCT_1W= (TifPOSTBD_1W/Tif_predicted_1W)*100) %>%
  mutate(TifPOSTBD_PCT_2W= (TifPOSTBD_2W/Tif_predicted_2W)*100) %>%
  mutate(TifPOSTBD_PCT_4W= (TifPOSTBD_4W/Tif_predicted_4W)*100) %>%
  mutate(TifPOSTBD_PCT_8W= (TifPOSTBD_8W/Tif_predicted_8W)*100) %>%
  mutate(TifPOSTBD_PCT_12W= (TifPOSTBD_12W/Tif_predicted_12W)*100) %>%
  mutate(TifPOSTBD_PCT_16W= (TifPOSTBD_16W/Tif_predicted_16W)*100) %>%
  mutate(TifPOSTBD_PCT_20W= (TifPOSTBD_20W/Tif_predicted_20W)*100) %>%
  mutate(TifPOSTBD_PCT_21W= (TifPOSTBD_21W/Tif_predicted_21W)*100) %>%
  mutate(TifPOSTBD_PCT_24W= (TifPOSTBD_24W/Tif_predicted_24W)*100) %>%
  mutate(TifPOSTBD_PCT_26W= (TifPOSTBD_26W/Tif_predicted_26W)*100) %>%
  mutate(TifPOSTBD_PCT_28W= (TifPOSTBD_28W/Tif_predicted_28W)*100) %>%
  mutate(TifPOSTBD_PCT_32W= (TifPOSTBD_32W/Tif_predicted_32W)*100) %>%
  mutate(TifPOSTBD_PCT_36W= (TifPOSTBD_36W/Tif_predicted_36W)*100) %>%
  mutate(TifPOSTBD_PCT_40W= (TifPOSTBD_40W/Tif_predicted_40W)*100) %>%
  mutate(TifPOSTBD_PCT_44W= (TifPOSTBD_44W/Tif_predicted_44W)*100) %>%
  mutate(TifPOSTBD_PCT_48W= (TifPOSTBD_48W/Tif_predicted_48W)*100) %>%
  mutate(TifPOSTBD_PCT_52W= (TifPOSTBD_52W/Tif_predicted_52W)*100) %>%
  
  
  mutate(TifPOSTBD_1W= FEV1POSTBD_L_1W/FVCPOSTBD_L_1W) %>%
  mutate(TifPOSTBD_2W= FEV1POSTBD_L_2W/FVCPOSTBD_L_2W) %>%
  mutate(TifPOSTBD_4W= FEV1POSTBD_L_4W/FVCPOSTBD_L_4W) %>%
  mutate(TifPOSTBD_8W= FEV1POSTBD_L_8W/FVCPOSTBD_L_8W) %>%
  mutate(TifPOSTBD_12W= FEV1POSTBD_L_12W/FVCPOSTBD_L_12W) %>%
  mutate(TifPOSTBD_16W= FEV1POSTBD_L_16W/FVCPOSTBD_L_16W) %>%
  mutate(TifPOSTBD_20W= FEV1POSTBD_L_20W/FVCPOSTBD_L_20W) %>%
  mutate(TifPOSTBD_21W= FEV1POSTBD_L_21W/FVCPOSTBD_L_21W) %>%
  mutate(TifPOSTBD_24W= FEV1POSTBD_L_24W/FVCPOSTBD_L_24W) %>%
  mutate(TifPOSTBD_26W= FEV1POSTBD_L_26W/FVCPOSTBD_L_26W) %>%
  mutate(TifPOSTBD_28W= FEV1POSTBD_L_28W/FVCPOSTBD_L_28W) %>%
  mutate(TifPOSTBD_32W= FEV1POSTBD_L_32W/FVCPOSTBD_L_32W) %>%
  mutate(TifPOSTBD_36W= FEV1POSTBD_L_36W/FVCPOSTBD_L_36W) %>%
  mutate(TifPOSTBD_40W= FEV1POSTBD_L_40W/FVCPOSTBD_L_40W) %>%
  mutate(TifPOSTBD_44W= FEV1POSTBD_L_44W/FVCPOSTBD_L_44W) %>%
  mutate(TifPOSTBD_48W= FEV1POSTBD_L_48W/FVCPOSTBD_L_48W) %>%
  mutate(TifPOSTBD_52W= FEV1POSTBD_L_52W/FVCPOSTBD_L_52W)

#======================================================================================================================================
#Preparation of different dataset - ACQ
Raw_ACQ_FEV1_without_Sequential_and_trial <- Raw_ACQ_FEV1 %>%
  dplyr::select(
    -Sequential_number,
    -Enrolled_Trial_name,
    -FEV1PREBD_L_52W,
    -FEV1PREBD_PCT_52W,
    -FVCPREBD_L_52W,
    -FVCPREBD_PCT_52W,
    -FEV1POSTBD_L_52W,
    -FEV1POSTBD_PCT_52W,
    -FVCPOSTBD_L_52W,
    -FVCPOSTBD_PCT_52W
  )

colnames(Raw_ACQ_FEV1)
Data_Oracle$Subject_ID

Merge_ACQ_with_ORACLE <- Data_Oracle %>%
  left_join(Raw_ACQ_FEV1_without_Sequential_and_trial, by = "Subject_ID")
summary(Merge_ACQ_with_ORACLE$ACQ5_score_mean_24W)

#ACQ : Study selection
Studies_included_ACQ<- unique(Raw_ACQ_FEV1$Enrolled_Trial_name)
Data_Oracle_ACQ<-Merge_ACQ_with_ORACLE %>%
  filter(Enrolled_Trial_name%in%Studies_included_ACQ)

unique(Data_Oracle_ACQ$Enrolled_Trial_name)
nrow(Data_Oracle_ACQ)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Preparation for analysis of ACQ placebo change

Data_Oracle_ACQ$ACQ_baseline_score_mean

# Calculate ACQ change using mutate with case_when
Data_Oracle_ACQ <- Data_Oracle_ACQ %>%
  mutate(
    # Determine which follow-up measurement to use with case_when
    ACQ5_followup = case_when(
      !is.na(ACQ5_score_mean_24W) ~ ACQ5_score_mean_24W,
      !is.na(ACQ5_score_mean_26W) ~ ACQ5_score_mean_26W,
      !is.na(ACQ5_score_mean_26W) ~ ACQ5_score_mean_22W,
      !is.na(ACQ5_score_mean_26W) ~ ACQ5_score_mean_23W,
      !is.na(ACQ5_score_mean_26W) ~ ACQ5_score_mean_25W,
      TRUE ~ NA_real_  # If neither is available
    ),
    
    # Calculate change from baseline
    delta_ACQ = case_when(
      !is.na(ACQ5_followup) & !is.na(ACQ_baseline_score_mean) ~  ACQ_baseline_score_mean-ACQ5_followup,
      TRUE ~ NA_real_
    )
  )

cols <- c(
  "ACQ5_score_mean_22W", 
  "ACQ5_score_mean_23W", 
  "ACQ5_score_mean_24W", 
  "ACQ5_score_mean_25W", 
  "ACQ5_score_mean_26W"
)

# Count NAs per column
sum(rowSums(is.na(Data_Oracle_ACQ[ , cols])) == 5)
#Data_Oracle_ACQ_filterNA <- Data_Oracle_ACQ %>%
#filter(rowSums(is.na(across(all_of(cols)))) < 5)

nrow(Data_Oracle_ACQ)

#Data_Oracle_ACQ_filterNA<-Data_Oracle_ACQ %>%
#filter(!is.na(delta_ACQ))

#Put Treatment step in ordinal
Data_Oracle_ACQ$Treatment_step_ordinal <- as.numeric(gsub("Step ", "", as.character(Data_Oracle_ACQ$Treatment_step)))
summary(Data_Oracle_ACQ$Treatment_step_ordinal)

#======================================================================================================================================
#======================================================================================================================================
##  FOR ASTHMA ATTACKS
#======================================================================================================================================
# Figure S19.A - Simple model
#======================================================================================================================================


Univariate_deltaattack_model <- MI_estimates(
  data = Data_Oracle_attack,
  outcome_var = "Delta_attack",
  predictor_vars = c("Age_per_10_imputated", "Gender", "BMI_per_5_imputated", 
                     "Smoking_Statut_imputated", "Atopy_history_imputated",
                     "Airborne_allergen_sensibilisation_imputated", "Allergic_rhinitis_imputated", 
                     "Eczema_imputated", "CRSsNP_imputated", "CRSwNP_imputated", 
                     "Treatment_step_ordinal", "ACQ_baseline_score_mean_imputated", 
                     "Attack_12mo_Nb_imputated", "FEV1_preBD_per10_Baseline_imputated", 
                     "FEV1_per10_reversibilityBD_imputated", "Eosinophils_Log_imputated",
                     "FeNO_Log_imputated", "IgE_Log_imputated"),
  covariables = c("Attack_12mo_Nb_imputated"),
  imp_col = ".imp",
  followup_offset = "Yes",
  followup_col = "Follow_up_duration_days",
  random_intercept = "Yes",
  random_intercept_var = "Enrolled_Trial_name",
  model_type = "lm"
)

Univariate_deltaattack_predictors <- attr(Univariate_deltaattack_model, "combined_results")
Univariate_deltaattack_predictors
row.names(Univariate_deltaattack_predictors) <- c(
  "Age (per 10 years)", "Sex (male)", "BMI (per 5 kg/m2)", "Smoking history (Ex vs no)",
  "Atopy history", "Airborne allergen sensitization", "Allergic rhinitis", "Eczema",
  "CRSsNP", "CRSwNP","Treatment step (per step increase)", "ACQ-5",
  "Attack history", "FEV1 pre-BD (per 10% decrease)", "FEV1 reversibility (per 10%)",
  "Eosinophils (Log)", "FeNO (Log)", "IgE (Log)"
)


forplo(
  as.data.frame(Univariate_deltaattack_predictors[, c("estimate", "2.5 %", "97.5 %")]),
  row.labels = row.names(Univariate_deltaattack_predictors),
  left.align = FALSE,
  shade.every = 1,
  shade.col = 'gray',
  linreg= TRUE,
  left.bar = FALSE,
  margin.left = 14,
  margin.right = 12,
  groups = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 3)),
  grouplabs = c("Demographic", "Comorbidities", "Asthma history", "Baseline lung function", "Inflammatory markers"),
  add.arrow.right=TRUE,
  arrow.right.length=3
)

#======================================================================================================================================
# Figure S19.B - Complete model (with covariates)
#======================================================================================================================================
Multivariate_deltaattack_model <- MI_estimates(
  data = Data_Oracle_attack,
  outcome_var = "Delta_attack",
  predictor_vars = c("Age_per_10_imputated", "Gender", "BMI_per_5_imputated", 
                     "Smoking_Statut_imputated", "Atopy_history_imputated",
                     "Airborne_allergen_sensibilisation_imputated", "Allergic_rhinitis_imputated", 
                     "Eczema_imputated", "CRSsNP_imputated", "CRSwNP_imputated", 
                     "Treatment_step_ordinal", "ACQ_baseline_score_mean_imputated", 
                     "Attack_12mo_Nb_imputated", "FEV1_preBD_per10_Baseline_imputated", 
                     "FEV1_per10_reversibilityBD_imputated", "Eosinophils_Log_imputated",
                     "FeNO_Log_imputated", "IgE_Log_imputated"),
  covariables = c(  "Attack_12mo_Nb_imputated",
                    "Gender",
                    "Smoking_Statut_imputated",
                    "Allergic_rhinitis_imputated",
                    "CRSsNP_imputated",
                    "Treatment_step_ordinal",
                    "ACQ_baseline_score_mean_imputated",
                    "FEV1_preBD_per10_Baseline_imputated",
                    "Eosinophils_Log_imputated",
                    "FEV1_per10_reversibilityBD_imputated",
                    "FeNO_Log_imputated"),
  imp_col = ".imp",
  followup_offset = "Yes",
  followup_col = "Follow_up_duration_days",
  random_intercept = "Yes",
  random_intercept_var = "Enrolled_Trial_name",
  model_type = "lm"
)

#Old covariables : "Gender","Atopy_history_imputated","CRSsNP_imputated","CRSwNP_imputated","Treatment_step_ordinal", "ACQ_baseline_score_mean_imputated","Attack_12mo_Nb_imputated","FEV1_preBD_per10_Baseline_imputated","Eosinophils_Log_imputated","FeNO_Log_imputated"

Multivariate_deltaattack_predictors <- attr(Multivariate_deltaattack_model, "combined_results")
Multivariate_deltaattack_predictors
row.names(Multivariate_deltaattack_predictors) <- c(
  "Age (per 10 years)", "Sex (male)", "BMI (per 5 kg/m2)", "Smoking history (Ex vs no)",
  "Atopy history", "Airborne allergen sensitization", "Allergic rhinitis", "Eczema",
  "CRSsNP", "CRSwNP","Treatment step (per step increase)", "ACQ-5",
  "Attack history", "FEV1 pre-BD (per 10% decrease)", "FEV1 reversibility (per 10%)",
  "Eosinophils (Log)", "FeNO (Log)", "IgE (Log)"
)


forplo(
  as.data.frame(Multivariate_deltaattack_predictors[, c("estimate", "2.5 %", "97.5 %")]),
  row.labels = row.names(Multivariate_deltaattack_predictors),
  left.align = FALSE,
  shade.every = 1,
  shade.col = 'gray',
  linreg= TRUE,
  left.bar = FALSE,
  margin.left = 14,
  margin.right = 12,
  groups = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 3)),
  grouplabs = c("Demographic", "Comorbidities", "Asthma history", "Baseline lung function", "Inflammatory markers"),
  add.arrow.right=TRUE,
  arrow.right.length=3
)

#======================================================================================================================================
###FINAL DESIGN - (Model multivariable) ASTHMA ATTACKS


# Save the plot using base R functions (not ggsave)
png("/Users/joel/ORACLE - Placebo - Exacerbations/Project/Forrest_plot_Predictors_Attacks_multivariable_final.png",
    width = 15, height = 12, units = "in", res = 300)

pdf("/Users/joel/ORACLE - Placebo - Exacerbations/Project/Forrest_plot_Predictors_Attacks_multivariable_final.pdf",
    width = 15, height = 12)

Forrest_plot_Predictors_Attacks_multivariable_final<-forplo(as.data.frame(Multivariate_deltaattack_predictors[, c("estimate", "2.5 %", "97.5 %")]),
                                                            
                                                            #Define the units
                                                            em="aRC",
                                                            
                                                            #Define the labels
                                                            row.labels = row.names(Multivariate_deltaattack_predictors),
                                                            left.align=FALSE,
                                                            
                                                            #Define the limits of the X axis
                                                            #xlim= c(0.5,4),
                                                            
                                                            #Add centered title (will be added after plot creation)
                                                            # title functionality not available in forplo - add with title() function
                                                            
                                                            #Define the arrows
                                                            add.arrow.right=FALSE,
                                                            arrow.right.length=20,  # Reduced length to move arrow more to center
                                                            add.arrow.left=FALSE,
                                                            arrow.left.length=20,   # Reduced length to move arrow more to center
                                                            
                                                            #Remove the left bar
                                                            left.bar = FALSE,
                                                            
                                                            
                                                            #Define the shading
                                                            shade.every = 1,
                                                            shade.col = 'grey',
                                                            shade.alpha = 0.2,
                                                            
                                                            
                                                            #Define the margin (increased top margin for title)
                                                            margin.left=15,
                                                            margin.right=12,
                                                            
                                                            
                                                            
                                                            groups = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 3)),
                                                            grouplabs = c("Demographic", "Comorbidities", "Asthma history", "Baseline lung function", "Inflammatory biomarkers"),
                                                            
                                                            
                                                            
                                                            # Characteristics of the points (char is equivalent to pch)
                                                            char = 20,
                                                            size = 1.5,
                                                            col=c(rep("darkgreen",6),
                                                                  rep("purple",8),
                                                                  rep("darkred",5),
                                                                  rep("darkblue",4),
                                                                  rep("darkorange",4)
                                                            ),
                                                            
                                                            #Adding the p-value
                                                            #pval=p_round(Multivariable_deltaFEV1_mL_compilation[,"p.value"], digits = 2),
                                                            
                                                            
                                                            #Adding AICc
                                                            #add.columns=round(Models_preBD_Z[,"AICC_mean"], 0),
                                                            #add.colnames=c('AICc'),
                                                            
                                                            
)

dev.off()


#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S19. LUNG FUNCTION - Comparison with Results from Original imputation (Meulmeester et al., Lancet Res, 2025)
# Old imputation vs. NEW imputation (by SML)
#======================================================================================================================================

#======================================================================================================================================
# Figure S19.A - Simple model
#======================================================================================================================================

Univariate_deltaFEV1_mL_model <- MI_estimates(
  data = Data_Oracle_lung,
  outcome_var = "delta_FEV1_preBD_mL",
  predictor_vars = c("Age_per_10_imputated", "Gender", "BMI_per_5_imputated", 
                     "Smoking_Statut_imputated", "Atopy_history_imputated",
                     "Airborne_allergen_sensibilisation_imputated", "Allergic_rhinitis_imputated", 
                     "Eczema_imputated", "CRSsNP_imputated", "CRSwNP_imputated", 
                     "Treatment_step_ordinal", "ACQ_baseline_score_mean_imputated", 
                     "Attack_12mo_Nb_imputated", "FEV1_preBD_per10_Baseline_imputated", 
                     "FEV1_per10_reversibilityBD_imputated", "Eosinophils_Log_imputated",
                     "FeNO_Log_imputated", "IgE_Log_imputated"),
  covariables = c("FEV1_preBD_per10_Baseline_imputated"),
  imp_col = ".imp",
  followup_offset = "No",
  random_intercept = "Yes",
  random_intercept_var = "Enrolled_Trial_name",
  model_type = "lm"
)
Univariate_deltaFEV1_mL_model
Univariate_deltaFEV1_mL_predictors <- attr(Univariate_deltaFEV1_mL_model, "combined_results")
Univariate_deltaFEV1_mL_predictors
row.names(Univariate_deltaFEV1_mL_predictors) <- c(
  "Age (per 10 years)", "Sex (male)", "BMI (per 5 kg/m2)", "Smoking history (Ex vs no)",
  "Atopy history", "Airborne allergen sensitization", "Allergic rhinitis", "Eczema",
  "CRSsNP", "CRSwNP","Treatment step (per step increase)", "ACQ-5",
  "Attack history", "FEV1 pre-BD (per 10% decrease)", "FEV1 reversibility (per 10%)",
  "Eosinophils (Log)", "FeNO (Log)", "IgE (Log)"
)


forplo(
  as.data.frame(Univariate_deltaFEV1_mL_predictors[, c("estimate", "2.5 %", "97.5 %")]),
  row.labels = row.names(Univariate_deltaFEV1_mL_predictors),
  left.align = FALSE,
  shade.every = 1,
  shade.col = 'gray',
  linreg= TRUE,
  left.bar = FALSE,
  margin.left = 14,
  margin.right = 12,
  groups = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 2)),
  grouplabs = c("Demographic", "Comorbidities", "Asthma history", "Baseline lung function", "Inflammatory markers"),
  add.arrow.right=TRUE,
  arrow.right.length=3
)

#======================================================================================================================================
# Figure S19.B - Complete model (with covariates)
#======================================================================================================================================

Multivariate_deltaFEV1_mL_model <- MI_estimates(
  data = Data_Oracle_lung,
  outcome_var = "delta_FEV1_preBD_mL",
  predictor_vars = c("Age_per_10_imputated", "Gender", "BMI_per_5_imputated", 
                     "Smoking_Statut_imputated", "Atopy_history_imputated",
                     "Airborne_allergen_sensibilisation_imputated", "Allergic_rhinitis_imputated", 
                     "Eczema_imputated", "CRSsNP_imputated", "CRSwNP_imputated", 
                     "Treatment_step_ordinal", "ACQ_baseline_score_mean_imputated", 
                     "Attack_12mo_Nb_imputated", "FEV1_preBD_per10_Baseline_imputated", 
                     "FEV1_per10_reversibilityBD_imputated", "Eosinophils_Log_imputated",
                     "FeNO_Log_imputated", "IgE_Log_imputated"),
  covariables = c(    "Age_per_10_imputated",
                      "Allergic_rhinitis_imputated",
                      "Treatment_step_ordinal",
                      "FEV1_preBD_per10_Baseline_imputated",
                      "FEV1_per10_reversibilityBD_imputated",
                      "FeNO_Log_imputated"),
  imp_col = ".imp",
  followup_offset = "No",
  random_intercept = "Yes",
  random_intercept_var = "Enrolled_Trial_name",
  model_type = "lm"
)

#Old imputation covariates : "Age_per_10_imputated","Gender", "BMI_per_5_imputated","Smoking_Statut_imputated","Atopy_history_imputated","Allergic_rhinitis_imputated","Treatment_step_ordinal", "FEV1_preBD_per10_Baseline_imputated","FEV1_per10_reversibilityBD_imputated","Eosinophils_Log_imputated","FeNO_Log_imputated"

Multivariate_deltaFEV1_mL_model_predictors <- attr(Multivariate_deltaFEV1_mL_model, "combined_results")
Multivariate_deltaFEV1_mL_model_predictors
row.names(Multivariate_deltaFEV1_mL_model_predictors) <- c(
  "Age (per 10 years)", "Sex (male)", "BMI (per 5 kg/m2)", "Smoking history (Ex vs no)",
  "Atopy history", "Airborne allergen sensitization", "Allergic rhinitis", "Eczema",
  "CRSsNP", "CRSwNP","Treatment step (per step increase)", "ACQ-5",
  "Attack history", "FEV1 pre-BD (per 10% decrease)", "FEV1 reversibility (per 10%)",
  "Eosinophils (Log)", "FeNO (Log)", "IgE (Log)"
)


forplo(
  as.data.frame(Multivariate_deltaFEV1_mL_model_predictors[, c("estimate", "2.5 %", "97.5 %")]),
  row.labels = row.names(Multivariate_deltaFEV1_mL_model_predictors),
  left.align = FALSE,
  shade.every = 1,
  shade.col = 'gray',
  linreg= TRUE,
  left.bar = FALSE,
  margin.left = 14,
  margin.right = 12,
  groups = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 3)),
  grouplabs = c("Demographic", "Comorbidities", "Asthma history", "Baseline lung function", "Inflammatory markers"),
  add.arrow.right=TRUE,
  arrow.right.length=3
)

#======================================================================================================================================
#======================================================================================================================================
# FINAL DESIGN - Model multivariable (deltaFEV1)

# Save the plot using base R functions (not ggsave)
png("/Users/joel/ORACLE - Placebo - Exacerbations/Project/Forrest_plot_Predictors_Lung_multivariable_final.png",
    width = 16, height = 12, units = "in", res = 300)

pdf("/Users/joel/ORACLE - Placebo - Exacerbations/Project/Forrest_plot_Predictors_Lung_multivariable_final.pdf",
    width = 16, height = 12)

Forrest_plot_Predictors_Lung_multivariable_final<-forplo(as.data.frame(Multivariate_deltaFEV1_mL_model_predictors[, c("estimate", "2.5 %", "97.5 %")]),
                                                         
                                                         #Define the units
                                                         em="aRC",
                                                         
                                                         #Define the labels
                                                         row.labels = row.names(Multivariate_deltaFEV1_mL_model_predictors),
                                                         left.align=FALSE,
                                                         
                                                         #Define the limits of the X axis
                                                         #xlim= c(0.5,4),
                                                         
                                                         #Add centered title (will be added after plot creation)
                                                         # title functionality not available in forplo - add with title() function
                                                         
                                                         #Define the arrows
                                                         add.arrow.right=FALSE,
                                                         arrow.right.length=20,  # Reduced length to move arrow more to center
                                                         add.arrow.left=FALSE,
                                                         arrow.left.length=FALSE,   # Reduced length to move arrow more to center
                                                         
                                                         #Remove the left bar
                                                         left.bar = FALSE,
                                                         
                                                         
                                                         #Define the shading
                                                         shade.every = 1,
                                                         shade.col = 'grey',
                                                         shade.alpha = 0.2,
                                                         
                                                         
                                                         #Define the margin (increased top margin for title)
                                                         margin.left=16,
                                                         margin.right=13,
                                                         
                                                         
                                                         
                                                         groups = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 3)),
                                                         grouplabs = c("Demographic", "Comorbidities", "Asthma history", "Baseline lung function", "Inflammatory biomarkers"),
                                                         
                                                         
                                                         
                                                         # Characteristics of the points (char is equivalent to pch)
                                                         char = 20,
                                                         size = 1.5,
                                                         col=c(rep("darkgreen",6),
                                                               rep("purple",8),
                                                               rep("darkred",5),
                                                               rep("darkblue",4),
                                                               rep("darkorange",4)
                                                         ),
                                                         
                                                         #Adding the p-value
                                                         #pval=p_round(Multivariable_deltaFEV1_mL_compilation[,"p.value"], digits = 2),
                                                         
                                                         
                                                         #Adding AICc
                                                         #add.columns=round(Models_preBD_Z[,"AICC_mean"], 0),
                                                         #add.colnames=c('AICc'),
                                                         
                                                         
)


dev.off()



#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# Figure S20. ACQ-5 - Comparison with Results from Original imputation (Meulmeester et al., Lancet Res, 2025)
# Old imputation vs. NEW imputation (by SML)
#======================================================================================================================================

#======================================================================================================================================
# Figure S20.A - Simple model
#======================================================================================================================================

Univariate_deltaACQ_model <- MI_estimates(
  data = Data_Oracle_ACQ,
  outcome_var = "delta_ACQ",
  predictor_vars = c("Age_per_10_imputated", "Gender", "BMI_per_5_imputated", 
                     "Smoking_Statut_imputated", "Atopy_history_imputated",
                     "Airborne_allergen_sensibilisation_imputated", "Allergic_rhinitis_imputated", 
                     "Eczema_imputated", "CRSsNP_imputated", "CRSwNP_imputated", 
                     "Treatment_step_ordinal", "ACQ_baseline_score_mean_imputated", 
                     "Attack_12mo_Nb_imputated", "FEV1_preBD_per10_Baseline_imputated", 
                     "FEV1_per10_reversibilityBD_imputated", "Eosinophils_Log_imputated",
                     "FeNO_Log_imputated", "IgE_Log_imputated"),
  covariables = c("ACQ_baseline_score_mean_imputated"),
  imp_col = ".imp",
  followup_offset = "No",
  random_intercept = "Yes",
  random_intercept_var = "Enrolled_Trial_name",
  model_type = "lm"
)
Univariate_deltaACQ_predictors <- attr(Univariate_deltaACQ_model, "combined_results")
Univariate_deltaACQ_predictors
row.names(Univariate_deltaACQ_predictors) <- c(
  "Age (per 10 years)", "Sex (male)", "BMI (per 5 kg/m2)", "Smoking history (Ex vs no)",
  "Atopy history", "Airborne allergen sensitization", "Allergic rhinitis", "Eczema",
  "CRSsNP", "CRSwNP","Treatment step (per step increase)", "ACQ-5",
  "Attack history", "FEV1 pre-BD (per 10% decrease)", "FEV1 reversibility (per 10%)",
  "Eosinophils (Log)", "FeNO (Log)", "IgE (Log)"
)


forplo(
  as.data.frame(Univariate_deltaACQ_predictors[, c("estimate", "2.5 %", "97.5 %")]),
  row.labels = row.names(Univariate_deltaACQ_predictors),
  left.align = FALSE,
  shade.every = 1,
  shade.col = 'gray',
  linreg= TRUE,
  left.bar = FALSE,
  margin.left = 14,
  margin.right = 12,
  groups = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 2)),
  grouplabs = c("Demographic", "Comorbidities", "Asthma history", "Baseline lung function", "Inflammatory markers"),
  add.arrow.right=TRUE,
  arrow.right.length=3
)

#======================================================================================================================================
# Figure S20.B - Complete model (with covariates)
#======================================================================================================================================
Data_Oracle_ACQ<-Data_Oracle_ACQ %>%
  filter(!is.na(delta_ACQ))
sum(is.na(Data_Oracle_ACQ$delta_ACQ))
summary(Data_Oracle_ACQ$delta_ACQ)

# MULTIVARIABLE - Delta ACQ-5 (using MI_estimates)

Multivariate_deltaACQ_model <- MI_estimates(
  data = Data_Oracle_ACQ,
  outcome_var = "delta_ACQ",
  predictor_vars = c("Age_per_10_imputated", "Gender", "BMI_per_5_imputated", 
                     "Smoking_Statut_imputated", "Atopy_history_imputated",
                     "Airborne_allergen_sensibilisation_imputated", "Allergic_rhinitis_imputated", 
                     "Eczema_imputated", "CRSsNP_imputated", "CRSwNP_imputated", 
                     "Treatment_step_ordinal", "ACQ_baseline_score_mean_imputated", 
                     "Attack_12mo_Nb_imputated", "FEV1_preBD_per10_Baseline_imputated", 
                     "FEV1_per10_reversibilityBD_imputated", "Eosinophils_Log_imputated",
                     "FeNO_Log_imputated", "IgE_Log_imputated"),
  covariables = c(    "Age_per_10_imputated",
                      "Treatment_step_ordinal",
                      "ACQ_baseline_score_mean_imputated",
                      "FEV1_preBD_per10_Baseline_imputated",
                      "FeNO_Log_imputated"),
  imp_col = ".imp",
  followup_offset = "No",
  random_intercept = "Yes",
  random_intercept_var = "Enrolled_Trial_name",
  model_type = "lm"
)

#Old imputation covariates : "Age_per_10_imputated","BMI_per_5_imputated","Treatment_step_ordinal","ACQ_baseline_score_mean_imputated","Attack_12mo_Nb_imputated","FEV1_preBD_per10_Baseline_imputated","FeNO_Log_imputated"

Multivariate_deltaACQ_model_predictors <- attr(Multivariate_deltaACQ_model, "combined_results")
Multivariate_deltaACQ_model_predictors
row.names(Multivariate_deltaACQ_model_predictors) <- c(
  "Age (per 10 years)", "Sex (male)", "BMI (per 5 kg/m2)", "Smoking history (Ex vs no)",
  "Atopy history", "Airborne allergen sensitization", "Allergic rhinitis", "Eczema",
  "CRSsNP", "CRSwNP","Treatment step (per step increase)", "ACQ-5",
  "Attack history", "FEV1 pre-BD (per 10% decrease)", "FEV1 reversibility (per 10%)",
  "Eosinophils (Log)", "FeNO (Log)", "IgE (Log)"
)


forplo(
  as.data.frame(Multivariate_deltaACQ_model_predictors[, c("estimate", "2.5 %", "97.5 %")]),
  row.labels = row.names(Multivariate_deltaACQ_model_predictors),
  left.align = FALSE,
  shade.every = 1,
  shade.col = 'gray',
  linreg= TRUE,
  left.bar = FALSE,
  margin.left = 14,
  margin.right = 12,
  groups = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 3)),
  grouplabs = c("Demographic", "Comorbidities", "Asthma history", "Baseline lung function", "Inflammatory markers"),
  add.arrow.right=TRUE,
  arrow.right.length=3
)

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
# FINAL DESIGN - Model multivariable (ACQ)

# Save the plot using base R functions (not ggsave)
png("/Users/joel/ORACLE - Placebo - Exacerbations/Project/Forrest_plot_Predictors_ACQ_multivariable_final.png",
    width = 16, height = 12, units = "in", res = 300)

pdf("/Users/joel/ORACLE - Placebo - Exacerbations/Project/Forrest_plot_Predictors_ACQ_multivariable_final.pdf",
    width = 16, height = 12)

Forrest_plot_Predictors_ACQ_multivariable_final<-forplo(as.data.frame(Multivariate_deltaACQ_model_predictors[, c("estimate", "2.5 %", "97.5 %")]),
                                                        
                                                        #Define the units
                                                        em="aRC",
                                                        
                                                        #Define the labels
                                                        row.labels = row.names(Multivariate_deltaACQ_model_predictors),
                                                        left.align=FALSE,
                                                        
                                                        #Define the limits of the X axis
                                                        xlim= c(-0.4,0.8),
                                                        
                                                        #Add centered title (will be added after plot creation)
                                                        # title functionality not available in forplo - add with title() function
                                                        
                                                        #Define the arrows
                                                        add.arrow.right=FALSE,
                                                        arrow.right.length=20,  # Reduced length to move arrow more to center
                                                        add.arrow.left=FALSE,
                                                        arrow.left.length=10,   # Reduced length to move arrow more to center
                                                        
                                                        #Remove the left bar
                                                        left.bar = FALSE,
                                                        
                                                        
                                                        #Define the shading
                                                        shade.every = 1,
                                                        shade.col = 'grey',
                                                        shade.alpha = 0.2,
                                                        
                                                        
                                                        #Define the margin (increased top margin for title)
                                                        margin.left=16,
                                                        margin.right=13,
                                                        
                                                        
                                                        
                                                        groups = c(rep(1, 4), rep(2, 6), rep(3, 3), rep(4, 2), rep(5, 3)),
                                                        grouplabs = c("Demographic", "Comorbidities", "Asthma history", "Baseline lung function", "Inflammatory biomarkers"),
                                                        
                                                        
                                                        
                                                        # Characteristics of the points (char is equivalent to pch)
                                                        char = 20,
                                                        size = 1.5,
                                                        col=c(rep("darkgreen",6),
                                                              rep("purple",8),
                                                              rep("darkred",5),
                                                              rep("darkblue",4),
                                                              rep("darkorange",4)
                                                        ),
                                                        
                                                        #Adding the p-value
                                                        #pval=p_round(Multivariable_deltaFEV1_mL_compilation[,"p.value"], digits = 2),
                                                        
                                                        
                                                        #Adding AICc
                                                        #add.columns=round(Models_preBD_Z[,"AICC_mean"], 0),
                                                        #add.colnames=c('AICc'),
                                                        
                                                        
)


dev.off()









