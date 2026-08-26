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
# === VERIFY DATASETS EXIST BEFORE LOADING === FROM SCRIPT : PlaceboAsthma_MAIN_Script_clean.R 

required_datasets <- c(
  "ORACLE_after_imputation",
  "ORACLE_Placebo_Attack_imputed",
  "ORACLE_Placebo_Lung_imputed",
  "ORACLE_Placebo_ACQ_imputed"
)

# Check which exist in current environment
existing    <- required_datasets[required_datasets %in% ls()]
missing     <- required_datasets[!required_datasets %in% ls()]

cat("✓ Found:\n")
print(existing)

cat("\n✗ Missing:\n")
print(missing)

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
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================

# SUPPLEMENTARY FIGURES

#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#======================================================================================================================================
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#======================================================================================================================================
# FIGURE S1 - FLOWCHART
#======================================================================================================================================

# RUN Main script before

# Combined 3-panel flowchart (A / B / C)

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


#======================================================================================================================================
#======================================================================================================================================
# Figure S2 - Fraction of NA by trial by variable
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
# Figure S3 - Distribution of placebo change - Bar graph
#======================================================================================================================================
# For Delta Asthma Attacks

#FIGURE S3.A - Bar graph

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

#FIGURE S3.B - Bar graph

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

#FIGURE S3.C - Bar graph

subset_3 <- Data_Oracle_ACQ %>%
  filter(.imp == 0, !is.na(delta_ACQ))

pdf("Figure_S3C_Histogram_DeltaACQ.pdf",
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
# Figure S4 - Q-Q PLOTS WITH STATISTICAL ANNOTATIONS
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
# FIGURE S4.A — Delta Asthma Attack Rate
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
  filename = paste0(base_path, "Figure_S4A_QQplot_DeltaAsthmaAttacks.pdf"),
  plot     = Fig_E2A,
  width    = 7,
  height   = 7,
  dpi      = 300
)

ggsave(
  filename = paste0(base_path, "Figure_S4A_QQplot_DeltaAsthmaAttacks.png"),
  plot     = Fig_E2A,
  width    = 7,
  height   = 7,
  dpi      = 300
)

cat("✓ Figure S4A saved.\n")

#======================================================================================================================================
# FIGURE S4.B — Delta FEV1
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

cat("✓ Figure S4B saved.\n")

#======================================================================================================================================
# FIGURE S4.C — Delta ACQ-5
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

cat("✓ Figure S4C saved.\n")

#======================================================================================================================================
# COMBINED FIGURE S4 — All 3 panels side by side
#======================================================================================================================================

library(patchwork)

Fig_E2_combined <- (Fig_E2A | Fig_E2B | Fig_E2C) +
  plot_annotation(
    title   = "Figure S4 — Normality Assessment of Placebo Change Distributions",
    caption = "Red dashed line = theoretical normal distribution. 
               Shapiro-Wilk test performed on imputation 0 (.imp == 0).
               For n > 5000, Shapiro-Wilk performed on a random sample of 5000 observations.",
    theme   = theme(
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.caption = element_text(size =  9, color = "grey50", hjust = 0.5)
    )
  )

ggsave(
  filename = paste0(base_path, "Figure_S4_QQplots_Combined.pdf"),
  plot     = Fig_E2_combined,
  width    = 21,
  height   = 8,
  dpi      = 300
)

ggsave(
  filename = paste0(base_path, "Figure_S4_QQplots_Combined.png"),
  plot     = Fig_E2_combined,
  width    = 21,
  height   = 8,
  dpi      = 300
)

cat("✓ Figure S4 combined saved.\n")

#======================================================================================================================================
#======================================================================================================================================
# FIGURE S5 : MODEL 1.A - LINEAR - ATTACKS - FOREST PLOT
#======================================================================================================================================
# SEE MAIN SCRIPT

#======================================================================================================================================
#======================================================================================================================================
# Figure S6 - Asthma attack - Forest plot’s results per trial + I2 (for FeNO and BEC)
#======================================================================================================================================

#======================================================================================================================================
# Figure S6.A - For FeNO (across trial)
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
# Figure S6.B - For BEC (across trial)
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
# Figure S7 - Asthma attack - Sensitivity analysis : Negative binomial regression
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
# Figure S7.A - Simple model
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
# Figure S7.B - Complete model (with covariates)
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
# Figure S8. Asthma attack - Sensitivity analysis – CAPTAIN + PACT 
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
#======================================================================================================================================
# Figure S10. Asthma attacks - Distribution (histogram) of follow-up duration by participants
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
# Figure S11. Asthma attacks - Sensitivity analysis for short follow-up duration 
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
# Figure S12. Asthma attacks - Sensitivity analysis for zero-inflation pre-trial exacerbations 
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
# FIGURE S13 : MODEL 1.A - LINEAR - ATTACKS - FOREST PLOT
#======================================================================================================================================
# SEE MAIN SCRIPT

#======================================================================================================================================
#======================================================================================================================================
# Figure S14 - Lung function - Forest plot’s results per trial + I2 (for FeNO)
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
# Figure S15. LUNG FUNCTION - Sensitivity analysis using ∆FEV₁ as % of predicted
#======================================================================================================================================
#======================================================================================================================================
# Figure S15.A - Simple model
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
# Figure S16. LUNG FUNCTION - Sensitivity analysis – CAPTAIN
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
# FIGURE S17 : MODEL 1.A - LINEAR - ATTACKS - FOREST PLOT
#======================================================================================================================================
# SEE MAIN SCRIPT

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







