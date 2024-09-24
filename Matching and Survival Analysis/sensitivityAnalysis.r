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
# Description: This script performs and saves results of sensitivity analyses.
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

run_cox_analysis <- function(df, surv_model) {
  res.cox <- coxph(surv_model, data = df)
  list(res = res.cox, summary = summary(res.cox), zph = cox.zph(res.cox))
}

analyze_cox_models <- function(cox_models) {
  # Initialize list to store results for each model
  table_rows <- vector("list", length = length(cox_models))
  names(table_rows) <- names(cox_models)
  
  # Iterate through each Cox model
  for (n in names(cox_models)) {
    model <- cox_models[[n]]
    
    # Extract proportional hazards test p-value
    zph <- model$zph[1]$table[, "p"]
    summary_model <- model$summary
    
    # Extract p-value and z-value, handling potential differences in model output
    p_value <- tryCatch(
      summary_model$coefficients[1, "Pr(>|z|)"],
      error = function(e) summary_model$coefficients[1, "p"]
    )
    z_value <- tryCatch(
      summary_model$coefficients[1, "z"],
      error = function(e) summary_model$coefficients[1, "coef"][1] / cox_models$all_frailty$summary$coefficients[1, "se(coef)"][1]
    )
    
    # Compile model statistics
    vals <- list(
      HR = round(summary_model$conf.int[1, "exp(coef)"], 2),
      SE_coef = round(summary_model$coefficients[1, "se(coef)"], 2),
      z_value = round(z_value, 2),
      p_value = p_value,
      CI = paste(round(summary_model$conf.int[1, "lower .95"], 2), "-", round(summary_model$conf.int[1, "upper .95"], 2)),
      zph = round(zph, 2),
      AIC = round(AIC(model$res), 0),
      BIC = round(BIC(model$res), 0)
    )
    
    # Combine statistics into a single string
    table_rows[[n]] <- paste(vals, collapse = ", ")
  }
  
  # Save sensitivity analysis results to file
  sink("results/sensitivity_analysis.txt")
  print(table_rows)
  sink()
  
  return(table_rows)
}

# =========================================================================
# Main code starts here
# =========================================================================

data_path <- "data/MatchedData.csv"
df <- read_and_process_data(data_path)

# Define models
cox_models_OS <- list(
  native = run_cox_analysis(df, as.formula(Surv(df$Survival, df$Dead) ~ FUS)),
  TumorSize = run_cox_analysis(df, as.formula(Surv(df$Survival, df$Dead) ~ FUS + TumorSize)),
  frailty = run_cox_analysis(df, as.formula(Surv(df$Survival, df$Dead) ~ FUS + frailty(subclass))),
  TumorSize_frailty = run_cox_analysis(df, as.formula(Surv(df$Survival, df$Dead) ~ FUS + TumorSize + frailty(subclass))),
  all = run_cox_analysis(df, as.formula(Surv(df$Survival, df$Dead) ~ FUS + TumorSize + Gender + Race)),
  all_frailty = run_cox_analysis(df, as.formula(Surv(df$Survival, df$Dead) ~ FUS + TumorSize + Gender + Race + frailty(subclass)))
)

cox_models_PFS <- list(
  native = run_cox_analysis(df, as.formula(Surv(PFS, Progression) ~ FUS)),
  TumorSize = run_cox_analysis(df, as.formula(Surv(PFS, Progression) ~ FUS + TumorSize)),
  frailty = run_cox_analysis(df, as.formula(Surv(PFS, Progression) ~ FUS + frailty(subclass))),
  TumorSize_frailty = run_cox_analysis(df, as.formula(Surv(PFS, Progression) ~ FUS + TumorSize + frailty(subclass))),
  all = run_cox_analysis(df, as.formula(Surv(PFS, Progression) ~ FUS + TumorSize + Gender + Race)),
  all_frailty = run_cox_analysis(df, as.formula(Surv(PFS, Progression) ~ FUS + TumorSize + Gender + Race + frailty(subclass)))
)

analyze_cox_models(cox_models_OS)
analyze_cox_models(cox_models_PFS)
