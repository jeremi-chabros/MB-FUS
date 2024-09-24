# =========================================================================
# Author: Jeremi Chabros
# Affiliations:
# 1) Department of Clinical Neurosciences, University of Cambridge
# 2) Computational Neuroscience Outcomes Center, Department of Neurosurgery,
#    Brigham & Women’s Hospital, Harvard Medical School
# 3) Department of Biostatistics; Department of Epidemiology,
#    Harvard T. H. Chan School of Public Health
# Date: 09/09/2023
# Last revised: 09/24/2024
# Description: This script reads and preprocesses matched patient data, generates
# and saves K-M survival curves, and performs Cox Proportional Hazards analysis.
# =========================================================================

# Load required libraries
{
  library(survival)
  library(lubridate)
  library(ggplot2)
  library(ggsurvfit)
  library(gtsummary)
  library(tidycmprsk)
  library(rms)
  library(survminer)
  library(adjustedCurves)
}

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
  return(df)
}

generate_survival_curve_plot <- function(df, fit_surv, vname) {
  if (vname == "PFS") {
    ylabel <- "PFS Probability"
  }
  if (vname == "OS") {
    ylabel <- "OS Probability"
  }
  plot <- ggsurvplot(
    fit_surv,
    data = df,
    palette = c("orange", "darkblue"),
    xlab = "Time (Months)",
    ylab = ylabel,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.height = 0.35,
    # surv.median.line = "hv",
    legend.labs = c("BWH", "BT008"),
    pval = TRUE,
    xlim = c(0, 70),
    ggtheme = theme(
      text = element_text(size = 20), # Base text size; affects titles and labels
      axis.title = element_text(size = 20), # X and Y axis titles
      axis.text = element_text(size = 16), # X and Y axis text (ticks)
      legend.title = element_text(size = 20), # Legend title
      legend.text = element_text(size = 20),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(colour = "black")
    )
    # xlim = c(0, 54.6)
  )
  x <- which(plot$plot$scales$find("x"))
  plot$plot$scales$scales[[x]] <- scale_x_continuous(breaks = seq(0, 70, by = 10))
  x <- which(plot$table$scales$find("x"))
  plot$table$scales$scales[[x]] <- scale_x_continuous(breaks = seq(0, 70, by = 10))

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

save_cox_plot <- function(fit, vname, supp_ext) {
  p_value <- summary(fit)$coefficients[, "Pr(>|z|)"][1]
  if (vname == "PFS") {
    ylabel <- "Adjusted PFS Probability"
    adjsurv <- adjustedsurv(
    data = df,
    variable = "FUS",
    ev_time = "PFS",
    event = "Progression",
    method = "direct",
    outcome_model = fit,
    conf_int = TRUE,
    mi_extrapolation = FALSE
  )
  }
  if (vname == "OS") {
    ylabel <- "Adjusted OS Probability"
    adjsurv <- adjustedsurv(
    data = df,
    variable = "FUS",
    ev_time = "Survival",
    event = "Dead",
    method = "direct",
    outcome_model = fit,
    conf_int = TRUE,
    mi_extrapolation = FALSE
  )
  }

  ptsv <- plot(adjsurv,
    conf_int = TRUE,
    labels = c("BWH", "BT008"),
    xlab = "Time (Months)",
    ylab = ylabel,
  ) + labs(x = "Time (Months)") +
    scale_color_manual(values = c("orange", "darkblue"), labels = c("BWH", "BT008")) +
    scale_fill_manual(values = c("orange", "darkblue"), labels = c("BWH", "BT008")) +
    guides(color = guide_legend(override.aes = list(fill = c("orange", "darkblue")))) +
    theme(legend.title = element_blank(), legend.position = "none", text = element_text(size = 20), plot.margin = margin(10, 10, 10, 10)) +
    annotate("text",
      x = 25,
      y = 0.1,
      label = paste0("p = ", round(p_value, 4)),
      hjust = 1, vjust = 1,
      size = 5,
      color = "black"
    ) + scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by = 10))
  ggsave(filename = paste0("results/", vname, " COX survival curve.svg"), plot = ptsv, width = 18, height = 9.75, dpi = 300, unit = "cm")
}


# =========================================================================
# Main code starts here
# =========================================================================

# Read and preprocess matched data
data_path <- "data/MatchedData.csv"
df <- read_and_process_data(data_path)

#  --------------------------------------------------------------------------------------------
# Generate and save K-M survival curve plots

# Overall Survival
fit_surv <- survfit2(Surv(Survival, Dead) ~ FUS, data = df)
fit_surv # print KM median OS
survdiff(Surv(Survival, Dead) ~ FUS, df)$pval # print p-val for log-rank test
surv_plot <- generate_survival_curve_plot(df, fit_surv, "OS")
save_plot(surv_plot, "results/OS_curves.svg")

# Progression-Free Survival
fit_surv <- survfit2(Surv(df$PFS, df$Progression) ~ FUS, data = df)
fit_surv # print KM median PFS
survdiff(Surv(df$PFS, df$Progression) ~ FUS, df)$pval # print p-val for log-rank test
surv_plot <- generate_survival_curve_plot(df, fit_surv, "PFS")
save_plot(surv_plot, "results/PFS_curves.svg")


#  --------------------------------------------------------------------------------------------
# Cox Proportional Hazards Model OS Analysis

surv_model <- as.formula(Surv(df$Survival, df$Dead) ~ FUS + TumorSize)
survdiff(surv_model, df)
survfit2(surv_model, data = df)
cox_results <- run_cox_analysis(df, surv_model)

# res.cox <- coxph(surv_model, data = df, weights = weights)
# list(res = res.cox, summary = summary(res.cox), zph = cox.zph(res.cox))

# Save Cox Analysis OS Results
sink("results/cox_results.txt")
print(cox_results)
sink()

# Cox Proportional Hazards Model PFS Analysis
surv_model <- as.formula(Surv(df$PFS, df$Progression) ~ FUS + TumorSize)
survdiff(surv_model, df)
survfit2(surv_model, data = df)
cox_results <- run_cox_analysis(df, surv_model)
AIC(cox_results$res)


# Save Cox Analysis PFS Results
sink("results/cox_pfs_results.txt")
print(cox_results)
sink()

# Save Cox-adjusted survival curves
df$FUS <- as.factor(df$FUS)
fit <- coxph(Surv(Survival, Dead) ~ FUS + TumorSize, data = df, x = TRUE)
save_cox_plot(fit, "OS", max(subset(df, FUS == 1)$Survival))

# df$Progression <- as.factor(df$Progression)
fit <- coxph(Surv(PFS, Progression) ~ FUS + TumorSize, data = df, x = TRUE)
save_cox_plot(fit, "PFS", max(subset(df, FUS == 1)$PFS))


# =========================================================================
