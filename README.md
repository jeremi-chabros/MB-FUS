# Microbubble-enhanced focused ultrasound and temozolomide for high-grade gliomas

This repository contains the code for genetic analysis, matching and survival analysis for 

*Microbubble-enhanced focused ultrasound and temozolomide for high-grade gliomas: safety, feasibility, efficacy, and sono-liquid biomarker analyses*.

The project is divided into two main branches:

1. Trajectory Figures (genetic analysis)
2. Matching and Survival Analysis

## 1. DNA Analysis Trajectories

### Step 1: C-Score Analysis
Run the R file named `c-score` to obtain the two groups of outcome of interest for cfDNA (EP and LTS).

### Step 2: F-Score Analysis
Run the R file named `f-score` to obtain the two groups of outcome of interest for FLRatio (EP and LTS).

### Step 3: Combined Plot Generation
Run the R file named `all in one plot` to generate the trajectory plot for all four groups.

### Step 4: cfDNA Last Observed Box Plot
Run the R file named `peak vs last and time trend for cfdna` to generate the last observed box plot for cfDNA, including *active tumors*, *early progressors*, and *long-term survivors*.

### Step 5: FLratio Last Observed Box Plot
Run the R file named `peak vs last and time trend for fragment` to generate the last observed box plot for FLratio, including *active tumors*, *early progressors*, and *long-term survivors*.

## 2. Matching and Survival Analysis

This section contains data and code necessary to perform coarsened exact matching and survival analysis.

### Data
- Patient datasheet: `data/FinalData.csv`

### Analysis Steps

1. **Matching Analysis**
   Run `code/matchingCEM.r` to generate:
   - `data/MatchedData.csv`: Final matched data sheet
   - `results/CEM_baseline.txt`: Summary of balance for unmatched data
   - `results/CEM_summary.txt`: Summary of balance for matched data
   - `results/matching_check.png`: Comparison density plot of variables pre- and post-matching

2. **Survival Analysis**
   Run `code/survivalAnalysis.r` to generate:
   - `results/PFS_curves.svg`: Plot of Kaplan-Meier PFS curves
   - `results/OS_curves.svg`: Plot of Kaplan-Meier OS curves
   - `results/PFS COX survival curve.svg`: Plot of Cox-adjusted PFS curves
   - `results/OS COX survival curve.svg`: Plot of Cox-adjusted OS curves
   - `results/cox_results.txt`: Cox proportional hazards regression results, including proportional hazards assumption test results

3. **Sensitivity Analysis**
   Run `code/sensitivityAnalysis.r` to generate:
   - `results/sensitivity_analysis_OS.txt`: Results of sensitivity analyses for OS
   - `results/sensitivity_analysis_PFS.txt`: Results of sensitivity analyses for PFS

## Usage

1. Clone this repository
2. Install the required R packages
3. Run the R scripts in the order specified above

## Citation

If you use this code or data in your research, please cite our paper:

[Citation details to be added upon publication]

## License

[Specify the license under which this project is released]

## Contact

For any questions or issues, please open an issue in this repository or contact Jeremi Chabros [jchabros@hsph.harvard.edu].
