# =========================================================================
# Author: Jeremi Chabros
# Affiliations:
# 1) University of Cambridge School of Clinical Medicine
# 2) Computational Neuroscience Outcomes Center, Department of Neurosurgery,
#    Brigham & Women’s Hospital, Harvard Medical School
# Date: 09/09/2023
# Description: This script reads and preprocesses matched patient data, generates
# and saves survival curves, and performs a Cox Proportional Hazards analysis.
# =========================================================================

# Load required libraries
library(survival)
library(lubridate)
library(ggplot2)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(rms)
library(survminer)
library(adjustedCurves)

# Set relative directory for the project
# setwd("path/to/your/directory")
setwd("/Users/jjc/Research/GBM-FUS/")

# =========================================================================
# Functions
# =========================================================================

read_and_process_data <- function(file_path) {
  df <- read.csv(file_path)
  df <- df[!(df$DeathCensorDate == "Unknown"), ]
  df$SurgeryDate <- dmy(df$SurgeryDate)
  df$DeathCensorDate <- dmy(df$DeathCensorDate)
  df$DiagnosisDate <- dmy(df$DiagnosisDate)
  df$ProgressionDate <- dmy(df$ProgressionDate)
  df$Survival <- as.duration(df$DiagnosisDate %--% df$DeathCensorDate) / dmonths(1)
  df$PFS <- as.duration(df$DiagnosisDate %--% df$ProgressionDate) / dmonths(1)
  df$LeadTime <- as.duration(df$SurgeryDate %--% df$DiagnosisDate) / dweeks(1)
  # df$Dead <- as.factor(df$Dead)
  # df$Progression <- as.factor(df$Progression)
  return(df)
}

generate_survival_curve_plot <- function(df, fit_surv) {
  plot <- ggsurvplot(
    fit_surv,
    data = df,
    palette = c("orange", "darkblue"),
    xlab = "Time (Months)",
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.height = 0.35,
    legend.labs = c("BWH", "BT008")
  )
  return(plot)
}

# Workaround to enable saving the survival table
save_plot <- function(plot, file_name) {
  ggsave_workaround <- function(g) {
    survminer:::.build_ggsurvplot(
      x = g,
      surv.plot.height = NULL,
      risk.table.height = NULL,
      ncensor.plot.height = NULL
    )
  }
  g_to_save <- ggsave_workaround(plot)
  ggsave(filename = file_name, plot = g_to_save, width = 18, height = 15, dpi = 300, unit = "cm")
}

run_cox_analysis <- function(df, surv_model) {
  res.cox <- coxph(surv_model, data = df)
  list(res = res.cox, summary = summary(res.cox), zph = cox.zph(res.cox))
}

save_cox_plot <- function(fit, vname) {
  adjsurv <- adjustedsurv(
    data = df,
    variable = "FUS",
    ev_time = "Survival",
    event = "Dead",
    method = "direct",
    outcome_model = fit,
    conf_int = TRUE
  )

  ptsv <- plot(adjsurv,
    conf_int = TRUE,
    labels = c("BWH", "BT008"),
    xlab = "Time (Months)",
  ) + labs(x = "Time (Months") +
    scale_color_manual(values = c("orange", "darkblue"), labels = c("BWH", "BT008")) +
    scale_fill_manual(values = c("orange", "darkblue"), labels = c("BWH", "BT008")) +
    guides(color = guide_legend(override.aes = list(fill = c("orange", "darkblue")))) +
    theme(legend.title = element_blank()) # Removes the legend title
  ggsave(filename = paste0("results/", vname, " COX survival curve.png"), plot = ptsv, width = 18, height = 15, dpi = 300, unit = "cm")
}

# =========================================================================
# Main code starts here
# =========================================================================

# Read and preprocess matched data
data_path <- "data/MatchedData.csv"
df <- read_and_process_data(data_path)

# Generate and save survival curve plot
fit_surv <- survfit2(Surv(Survival, Dead) ~ FUS, data = df)
surv_plot <- generate_survival_curve_plot(df, fit_surv)
save_plot(surv_plot, "results/OS_curves.png")

fit_surv <- survfit2(Surv(df$PFS, df$Progression) ~ FUS, data = df)
surv_plot <- generate_survival_curve_plot(df, fit_surv)
save_plot(surv_plot, "results/PFS_curves.png")

# Cox Proportional Hazards Model Analysis
# surv_model <- as.formula(Surv(df$Survival, df$Dead) ~ FUS + frailty(subclass)) # Shared frailty model
surv_model <- as.formula(Surv(df$Survival, df$Dead) ~ FUS)
survdiff(surv_model, df)
survfit2(Surv(df$Survival, df$Dead) ~ FUS, data = df)
cox_results <- run_cox_analysis(df, surv_model)

# Save Cox Analysis Results
sink("results/cox_results.txt")
print(cox_results)
sink()

# Progression Free Survival
surv_model <- as.formula(Surv(df$PFS, df$Progression) ~ FUS)
survdiff(surv_model, df)
survfit2(surv_model, data = df)
cox_results <- run_cox_analysis(df, surv_model)


df$FUS <- as.factor(df$FUS)
fit <- coxph(Surv(Survival, Dead) ~ FUS, data = df, x = TRUE)
save_cox_plot(fit, "OS")

# df$Progression <- as.factor(df$Progression)
fit <- coxph(Surv(PFS, Progression) ~ FUS, data = df, x = TRUE)
save_cox_plot(fit, "PFS")

# Save PFS Analysis Results
sink("results/cox_pfs_results.txt")
print(cox_results)
sink()

# =========================================================================
# Model comparison for sensitivity analysis

cox_fit1 <- run_cox_analysis(df, as.formula(Surv(df$Survival, df$Dead) ~ FUS))
cox_fit2 <- run_cox_analysis(df, as.formula(Surv(df$Survival, df$Dead) ~ FUS + IDH))

# sink("results/cox_results_IDH_covariate.txt") # Save the result with IDH as a co-variate
# print(cox_fit2)
# sink()

# Extract AIC and BIC from Model 1
aic1 <- AIC(cox_fit1$res)
bic1 <- BIC(cox_fit1$res)

# Extract AIC and BIC from Model 2
aic2 <- AIC(cox_fit2$res)
bic2 <- BIC(cox_fit2$res)

# Compare
print(paste("Model 1 - AIC: ", aic1, " BIC: ", bic1))
print(paste("Model 2 - AIC: ", aic2, " BIC: ", bic2))


control <- subset(df, FUS == 0)
treatment <- subset(df, FUS == 1)
summary(control)
summary(treatment)
