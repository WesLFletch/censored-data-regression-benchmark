# Censored Data Regression Benchmark Materials

Supplementry data and code materials for performing all parts of the benchmark study, including:

- Acquiring, manipulating, and analyzing the public bladder cancer dataset TCGA-BLCA
- Performing both settings of our simulation studies

## Quick Start Guide

The benchmark is split into three separate studies, two settings of simulation studies and an analysis of a real bladder cancer dataset (TCGA-BLCA). Setting-II of simulations and the real data analysis first require 1) acquiring the bladder cancer dataset and 2) performing our preliminary feature selection method on it to reduce data dimensionality and multicoinearity. Setting-I of simulations does not require these steps.

### Acquiring the Bladder Cancer Dataset (TCGA-BLCA)

1. All required code files to acquire and manipulate the bladder cancer dataset for use in the benchmark are found in `TCGA-BLCA`. Place all files in a common directory.

2. In `TCGA_BLCA.R`, modify the string on line 3 with the directory containing all required files.

3. Install all required R packages, the names of which are found on lines 6-8 of `TCGA_BLCA.R`. All packages are either hosted/archived on CRAN, or in GitHub repositories.

4. Run the entirety of `TCGA_BLCA.R` to download the required data for the study from the Xena browser into the working directory and manipulate its structure. The raw data will be stored in `blca_mrna.txt`, `blca_surv.txt`, and `blca_clin.txt`, and the joined/restructured data will be stored in `BLCA.RDS`. The basic structure of this object is described on lines 78-83 of `TCGA_BLCA.R`.

### Preliminary Feature Selection

1. All required code files to perform the preliminary feature selection are found in `preliminary-feature-selection`. Additionally, the data file `BLCA.RDS` from the above process is required. Place all files in a common directory.

2. In `PreliminaryFS.R` and `InterpretPFSResults.R`, modify the strings on line 3 with the directory containing all required files.

3. Install all required R packages, the names of which are found on lines 6-10 of `PreliminaryFS.R`. All packages are either hosted/archived on CRAN, or in GitHub repositories.

4. Run the entirety of `PreliminaryFS.R` to perform the preliminary feature selection, and record the results in `PrelimFSResults.RDS` in the working directory. **This process may take extremely long, even on high-end computing hardware.**

5. Run the entirety of `InterpretPFSResults.R` to acquire the names of the mRNA features that passed the preliminary feature selection process. This object will be stored in `UncorMRNA.RDS` in the working directory.

### Setting-I of Simulations

1. All required code files to run setting-i of simulations are found in `simulation-setting-i`. Place all required files in a common directory.

2. In `RunSettingI.R`, modify the string on line 3 with the directory containing all required files.

3. Install all required R packages, the names of which are found on lines 6-16 in `RunSettingI.R`. All packages are either hosted/archived on CRAN, or in GitHub repositories.

4. Set the simulation parameters. Lines 36-48 contain the data characteristics we examine in our studies, as well as their descriptions. When changing these values, be mindful that $\texttt{num\_trials}\cdot\texttt{length(s\_vals)}\cdot\texttt{length(a\_vals)}\cdot\texttt{length(s\_vals)}\cdot\texttt{length(sig\_str\_vals)}$ synthetic datasets are analyzed by all methods, each with $n$ observations and $p$ data features. Adding more elements to these vectors will greatly increase the duration of the simulations.

5. Run the entirety of `RunSettingI.R` to perform the simulation study. The results of the study will be stored in `setting-i-results.RDS` within the working directory once completed. The resulting object's structure is described on lines 55-59 in `RunSettingI.R`. **The simulations may take extremely long, even on high-end computing hardware.**

### Setting-II of Simulations

1. All required code files to run setting-ii of simulations are found in `simulation-setting-ii`. Additionally, the data files `BLCA.RDS` and `UncorMRNA.RDS` from above processes are required. Place all files in a common directory.

2. In `RunSettingII.R`, modify the string on line 3 with the directory containing all required files.

3. Install all required R packages, the names of which are found on lines 6-18 in `RunSettingII.R`. All packages are either hosted/archived on CRAN, or in GitHub repositories.

4. Set the simulation parameters. Lines 33-49 of `RunSettingII.R` contain parameters detailing the number of times to simulate survival times, the approximate censorship rates to examine, and how strong the signals of the $\beta$ vector are. When changing these values, be mindful that $\texttt{num\_reps}\cdot\texttt{length(quantiles)}$ simulations are performed.

5. Run the entirety of `RunSettingII.R` to perform the simulation study. The results of the study will be stored in `setting-ii-results.RDS` within the working directory once competed. The resulting object's structure is described on lines 128-131 in `RunSettingII.R`. **The simulations may take extremely long, even on high-end computing hardware.**

### Real Data Analysis (TCGA-BLCA)

1. All required code files to run the real data analysis are found in `real-data-analysis`. Additionally, the data files `BLCA.RDS` and `UncorMRNA.RDS` from above processes are required. Place all files in a common directory.

2. In `RunRealDataAnalysis.R`, modify the string on line 3 with the directory containing all required filed.

3. Install all required R packages, the names of which are found on lines 6-18 of `RunRealDataAnalysis.R`. All packages are either hosted/archived on CRAN, or in GitHub repositories.

4. Choose which methods to perform. lines 23-24 of `RunRealDataAnalysis.R` contains a string vector of all methods, remove a string to remove the corresponding method from the study.

5. Run the entirety of `RunRealDataAnalysis.R` to perform the real data analysis. The results of the study will be stored in `real-data-analysis-results.RDS` within the working directory once completed. The resulting object's structure is described on lines 30-34 of `RunRealDataAnalysis.R`. **The analysis may take some time, even on high-end computing hardware.**
