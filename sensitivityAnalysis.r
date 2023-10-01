# =========================================================================
# Author: Jeremi Chabros
# Affiliations:
# 1) University of Cambridge School of Clinical Medicine
# 2) Computational Neuroscience Outcomes Center, Department of Neurosurgery,
#    Brigham & Women’s Hospital, Harvard Medical School
# Date: 09/09/2023
# Description: This script performs bootstrapping to run multiple Cox Proportional Hazards
# models and calculates the proportion of significant p-values.
# =========================================================================

# Load required libraries
library(survival)
library(lubridate)
library(ggplot2)
library(tidyverse)

# Set relative directory for the project
# setwd("path/to/your/directory")
setwd("/Users/jjc/Research/GBM-FUS/")

# Function to convert dates and calculate survival time
read_and_process_data <- function(file_path) {
  df <- read.csv(file_path)
  df <- df[!(df$DeathCensorDate == "Unknown"), ]
  df$SurgeryDate <- dmy(df$SurgeryDate)
  df$DeathCensorDate <- dmy(df$DeathCensorDate)
  df$DiagnosisDate <- dmy(df$DiagnosisDate)
  df$Survival <- as.duration(df$DiagnosisDate %--% df$DeathCensorDate) / dmonths(1)
  return(df)
}

# Function to bootstrap and fit Cox models
bootstrap_cox_model <- function(control, treatment, n_iter = 1000) {
  results <- vector("list", length = n_iter)

  for (i in seq_len(n_iter)) {
    control_sample <- control[sample(nrow(control), nrow(treatment), replace = TRUE), ]
    temp_data <- rbind(control_sample, treatment)
    cox_model <- coxph(Surv(Survival, Dead) ~ FUS, data = temp_data)
    results[[i]] <- summary(cox_model)$coefficients
  }

  results
}

# Function to calculate p-values
calculate_p_values <- function(bootstrap_results) {
  sapply(bootstrap_results, function(x) x["FUS", "Pr(>|z|)"])
}

# =========================================================================
# Main code starts here
# =========================================================================

# Prepare data
df_prepared <- read_and_process_data("data/MatchedData.csv")

control <- subset(df_prepared, FUS == 0)
treatment <- subset(df_prepared, FUS == 1)

# Bootstrap Cox models
bootstrap_results <- bootstrap_cox_model(control, treatment)

# Calculate p-values
bootstrap_p_values <- calculate_p_values(bootstrap_results)

# Display results
print(paste("Median of p-values: ", median(bootstrap_p_values)))

# Plot histogram of p-values
png(file = "results/bootstrap_p_values.png", units = "in", width = 5, height = 5, res = 300)
hist(bootstrap_p_values,
  breaks = 40, xaxt = "n",
  xlab = "Bootstrap p-values", ylab = "Frequency", main = ""
)
axis(1, at = seq(0, 1, by = 0.05), las = 2)
abline(v = 0.05, col = "red")
dev.off()
