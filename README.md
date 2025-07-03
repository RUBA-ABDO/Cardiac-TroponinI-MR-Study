# Cardiac-TroponinI-MR-Study
This repository contains R scripts and data files for the master's thesis: Cardiac Troponin I and Risk of Coronary Artery Disease and Dilated Cardiomyopathy – A Mendelian Randomization Study. It includes all steps for preparing summary statistics,harmonizing datasets, and performing the MR analysis.
# Data 
Exposure data: Cardiac troponin I (cTnI) GWAS.
Outcome data: Coronary artery disease (CAD) and Dilated cardiomyopathy (DCM).
# Sources:
The CAD summary data were obtained from the CARDIoGRAMplusC4D consortium (Nikpay et al., 2015).The cTnI and DCM datasets were provided directly by the corresponding authors upon request.
# Notes:
Please ensure that the exposure and outcome GWAS datasets do not contain overlapping samples.
Use a consistent genome build (hg19/GRCh37 or liftOver to GRCh38) for all datasets.
# How to Reuse This Script
1) Replace file paths in fread() or read.csv() with your own exposure and outcome summary files.
2) Adapt column names in read_exposure_data() and read_outcome_data() for MR analysis format: SNP ID, effect size (beta), standard error, p-value, effect allele, other allele, EAF, sample size.
3) Confirm that SNP IDs and positions match the chosen genome build and are formatted consistently (rsID or chr:pos).
4) If SNPs are missing after harmonization, use LDlinkR to identify suitable proxies (LD r² > 0.8 recommended).
5) Review output files:
   write.csv(): saves harmonized or cleaned data.
   saveRDS(): saves R objects for later use.
   Plots: scatter, funnel, forest, and leave-one-out plots saved as PNG.
 6) Final MR results are printed in the console and can be saved to PDF or HTML using mr_report().
# Inputs
Exposure file: see df <- fread(...)
Outcome files: see CAD37 <- fread(...) or DCM37 <- fread(...)
Harmonization: see harmonise_data()
# Outputs
Clean .csv files: exposure/outcome ready for MR.
Harmonized datasets: datIC (CAD) and datID (DCM).
MR results: summary tables, ORs, CIs, and plots.
Final report: optional PDF/HTML from mr_report().
# RUBA ABDO #
