# =========================================================================
# Author: Jeremi Chabros
# Affiliations:
# 1) University of Cambridge School of Clinical Medicine
# 2) Computational Neuroscience Outcomes Center, Department of Neurosurgery,
#    Brigham & Women’s Hospital, Harvard Medical School
# Date: 09/09/2023
# Description: This script prepares data by matching cases using
# Coarsened Exact Matching in the MatchIt package, then saves the matched
# data, summary statistics, and density plots.
# To avoid errors, it is recommended to run line-by-line or run twice as whole.
# =========================================================================

# Load required libraries
library(MatchIt)
library("marginaleffects")
library(survival)
library(lubridate)
library(ggplot2)

# Set relative directory for the project
# setwd("path/to/your/directory")
setwd("/Users/jjc/Research/GBM-FUS/")

# =========================================================================
# Functions
# =========================================================================

# Define functions to read and preprocess data
read_and_process_data <- function(file_path) {
    df <- read.csv(file_path)
    # df <- na.omit(df)
    df <- subset(df, !(is.na(DeathCensorDate) | DeathCensorDate == "Unknown"))
    df <- subset(df, X1p19q != 1)
    df <- subset(df, Chemotherapy == 1)
    df <- subset(df, Radiotherapy == 1)
    # df <- subset(df, TumorSize != "Missing")
    # Dichotomize race
    df$Race <- ifelse(df$Race == "White", "White", "Non-white")
    return(df)
}

# Define function for matching
perform_matching <- function(df, formula, method, distance, cutpoints, ratio) {
    m.out <- matchit(formula,
        data = df,
        method = method,
        distance = distance,
        cutpoints = cutpoints,
        ratio = ratio
    )
    return(m.out)
}

# Function to save density plots
save_density_plot <- function(m.out, which_vars, file_path) {
    png(file = file_path, units = "in", width = 5, height = 5, res = 300)
    plot(m.out, type = "density", interactive = FALSE, which.xs = which_vars, main="")
    title(main = "")
    dev.off()
}

# =========================================================================
# Main code starts here
# =========================================================================

# Read and preprocess data
file_path <- "data/FinalData.csv"
df <- read_and_process_data(file_path)

trt <- subset(df, TumorSize != "Unknown")
trt$TumorSize <- as.numeric(trt$TumorSize)
df$TumorSize[df$TumorSize == "Unknown"] <- median(trt$TumorSize)
df$TumorSize[is.na(df$TumorSize)] <- as.numeric(median(trt$TumorSize))
df$TumorSize <- as.numeric(df$TumorSize)

# Define model formula
mymodel <- as.formula(FUS ~ Age + Gender + Race + TumorSize)
# mymodel <- as.formula(FUS ~ Age + Gender + Race + TumorSize + IDH)

m.out0 <- matchit(mymodel,
    data = df,
    method = NULL, distance = "glm"
)
m.out0

# Perform matching
cutpoints <- list(Age = c(0, 20, 40, 60, 80, 100), TumorSize = "q3")
ratio <- 1
m.out2 <- perform_matching(df, mymodel, "cem", "glm", cutpoints, ratio)

summary(m.out2)

# Save summary to a text file
sink("results/CEM_summary.txt")
summary(m.out2, un = FALSE)
sink()

control <- subset(match.data(m.out2), FUS == 0)
summary(control)

table(control$Gender)
table(control$Race)
table(control$Location)
table(control$MGMT)
table(control$IDH)

# Save density plots to PNG files
# save_density_plot(m.out2, ~ Age + Gender + Race, "results/matching_check1.png")
# save_density_plot(m.out2, ~TumorSize, "results/matching_check2.png")

# Save matched data to a CSV file
matched_data <- match.data(m.out2)
write.csv(matched_data, "data/MatchedData.csv", row.names = FALSE)

# =========================================================================
















