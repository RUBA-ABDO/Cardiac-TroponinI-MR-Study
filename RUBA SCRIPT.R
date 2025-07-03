# Title: Cardiac Troponin I and Risk of Coronary Artery Disease and Dilated Cardiomyopathy - a Mendelian Randomization Study
# Master's Thesis Script – RUBA ABDO

# INSTALL REQUIRED PACKAGES 

install.packages("stringi")          # For string manipulation (e.g., text parsing/formatting)
install.packages("devtools")         # To install packages from GitHub (e.g., TwoSampleMR)
install.packages("plyr")             # For data manipulation (older but still used in some functions)
install.packages("ggplot2")          # For data visualization and MR plot outputs
install.packages("LDlinkR")          # For finding proxy SNPs based on linkage disequilibrium
install.packages("tidyr")            # For reshaping and tidying data
install.packages("tinytex")          # Required to create PDF output using RMarkdown
tinytex::install_tinytex()           # Installs TinyTeX LaTeX engine for rendering PDFs
install.packages("rmarkdown")        # For writing reports in RMarkdown format
install.packages("knitr")            # Converts RMarkdown to final documents (e.g., PDF, HTML)
install.packages("markdown")         # Markdown support for rendering reports
install.packages("Cairo")            # For high-quality graphics, useful in saving plots for thesis

# LOAD LIBRARIES 

library(TwoSampleMR)                 # Core package for Mendelian Randomization analysis
library(plyr)                        # Used for data summarization and manipulation
library(ggplot2)                     # Plotting MR results (scatter, forest, funnel)
library(data.table)                 # Fast and efficient data reading/manipulation
library(rtracklayer)                # For genomic annotations if needed
library(GenomicRanges)              # For handling genomic interval data
library(stringr)                    # Modern string manipulation (more robust than base R functions)
library(dplyr)                      # For modern data manipulation (filter, select, mutate, etc.)
library(devtools)                   # Needed for installing packages from GitHub
library(LDlinkR)                    # Interface with LDlink API for LD proxy search
library(tidyr)                      # For reshaping MR results into tidy formats
library(stringr)                    # (Redundant but safe) for string operations
library(ggrepel)                    # Improves ggplot label readability (e.g., in MR scatter plots)

# INSTALL DEVELOPMENT VERSION OF TwoSampleMR FROM GITHUB 
install_github("MRCIEU/TwoSampleMR")  # Ensures latest development version is used for advanced MR features

# Load full summary statistics for cardiac troponin I GWAS
df <- fread("C:\\Users\\hp\\Downloads\\troponin_gwas_summary_statistics (1).txt")  # Read GWAS summary for troponin I
Troponin <- df  # Assign it to a named variable

# Load GWAS summary statistics for CAD outcome
file_path <- "C:/Users/hp/Desktop/CAD/cad.add.160614.website.txt"  # Set path to CAD data
CAD37 <- fread(file_path)  # Read CAD summary statistics

# Define list of selected SNPs for exposure (from Mokesnes.et.al article)
markers <- c(
  "3:12800305:A:G", "4:186218785:A:G", "5:177415473:T:C", "6:32715903:G:T",
  "10:74103992:C:T", "10:113588287:G:A", "11:22250324:A:T", "14:103137158:G:A",
  "15:84825188:A:T", "17:66212167:C:G", "18:36681510:G:A", "18:59353974:T:C"
)
# Filter exposure data to keep only the selected SNPs
filtered_df <- df[df$MarkerName %in% markers, ]  # Keep rows with matching SNPs
print(filtered_df)

# Rename columns to match MR package expectations
colnames(filtered_df)[colnames(filtered_df) == "MarkerName"] <- "SNP"
colnames(filtered_df)[colnames(filtered_df) == "Effect"] <- "beta"
colnames(filtered_df)[colnames(filtered_df) == "StdErr"] <- "SE"
colnames(filtered_df)[colnames(filtered_df) == "P-value"] <- "p.value"
colnames(filtered_df)[colnames(filtered_df) == "Freq1"] <- "FREQ"
colnames(filtered_df)[colnames(filtered_df) == "Allele1"] <- "effect_allele"
colnames(filtered_df)[colnames(filtered_df) == "Allele2"] <- "noneffect_allele"
filtered_df$samplesize <- 48115  # Add sample size from original GWAS

# Save prepared exposure data to working directory
setwd("C:/Users/hp/Desktop/last last last analysis")  # Set directory for saving output
write.csv(filtered_df, "ctni.csv", row.names = FALSE)  # Export filtered exposure data

# Read and format exposure data into TwoSampleMR structure
Ctni_exp_data <- read_exposure_data(
  filename = "ctni.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  pval_col = "p.value",
  eaf_col = "FREQ",
  effect_allele_col = "effect_allele",
  other_allele_col = "noneffect_allele",
  samplesize_col = "samplesize"
)
Ctni_exp_data$units.exposure <- "SD(inverse rank-normalized)" # Add unit type
Ctni_exp_data$exposure <- "cardiac troponin I"  # Label exposure for output

# Calculate r from beta and SE (for Rsq and F-statistics)
Ctni_exp_data$r.exposure <- get_r_from_bsen(
  b = Ctni_exp_data$beta.exposure,
  se = Ctni_exp_data$se.exposure,
  n = Ctni_exp_data$samplesize.exposure
)
Ctni_exp_data$rsq.exposure <- Ctni_exp_data$r.exposure^2  # Calculate R²
Ctni_exp_data$F_statistic <- with(Ctni_exp_data, (rsq.exposure * (samplesize.exposure - 2)) / (1 - rsq.exposure))  # F-statistic

# Review Rsq values
summary(Ctni_exp_data$rsq.exposure)

# Clean SNP names to remove allele info, keep only chr:position
Ctni_exp_data$SNP <- sub("(:[0-9]+):[a-zA-Z]+:[a-zA-Z]+", "\\1", Ctni_exp_data$SNP)

#Ensure "chr" prefix in chromosome column for CAD summary
if (!grepl("^chr", CAD37$chr[1])) CAD37[, chr := paste0("chr", chr)]

# Prepare BED file for liftOver from hg19 to hg38 (UCSC browser)
BED <- CAD37[, .(chr, start = bp_hg19 - 1L, end = bp_hg19, name = markername)]
fwrite(BED, "C:/Users/hp/Desktop/CAD37_hg19.bed", sep = "\t", col.names = FALSE)

# Read liftOver output (converted coordinates)
bed_data <- fread("C:/Users/hp/Desktop/hglft_genome_3838ea_3a5910.bed", header=FALSE)
setnames(bed_data, c("chr", "start", "end", "markername", "extra"))

# Convert both datasets to data.table for merging
setDT(CAD37)
setDT(bed_data)

# Merge liftOver coordinates into CAD GWAS data
bed_sub <- bed_data[, .(markername, chr_hg38 = chr, pos_hg38 = start)]
merged_data <- merge(CAD37, bed_sub, by = "markername", all.x = TRUE, sort = FALSE)

# Preview merged data to ensure correct mapping
head(merged_data[, .(effect_allele, noneffect_allele, chr_hg38, pos_hg38)])

# Remove redundant columns if present
merged_data[, c("POS_hg38", "CHR_hg38") := NULL]

# Adjust position +1 (UCSC coordinates start from 0)
merged_data[, pos_hg38 := pos_hg38 + 1]

# Remove "chr" prefix
merged_data[, chr := sub("^chr", "", chr)]

# Create harmonized SNP column as chr:pos
merged_data[, SNP := paste0(chr, ":", pos_hg38)]

# Save cleaned and updated CAD outcome data
fwrite(merged_data, "C:/Users/hp/Desktop/merged_data_cleaned.txt", sep = "\t")

# Rename columns for compatibility with TwoSampleMR
CAD_out_data<- merged_data
setnames(CAD_out_data, old = c("noneffect_allele", "se_dgc", "p_dgc", "effect_allele_freq","SNP","beta"),
         new = c("other_allele", "se", "pval", "eaf","SNP","beta"), skip_absent=TRUE)

# Select and format outcome-specific columns
CAD_out_data <- CAD_out_data[, .(SNP, effect_allele, other_allele, beta, se, pval, eaf)]
CAD_out_data$outcome <- "Coronary artery disease"  # Label for MR output
head(CAD_out_data)
dim(CAD_out_data)

# Define unit of measurement for outcome (log odds)
CAD_out_data$units.outcome <- "Logodd"

# Remove rows with missing SNPs
CAD38_subset_cleaned <- CAD_out_data %>% filter(!is.na(SNP))

# Remove duplicated SNPs
CAD38_subset_cleaned <- CAD38_subset_cleaned %>% distinct(SNP, .keep_all = TRUE)

# Save cleaned CAD outcome data
write.csv(CAD38_subset_cleaned, "CAD38_cleaned.csv", row.names = FALSE)
CAD_out_data <- CAD38_subset_cleaned

# Final rename of columns for outcome data as required by TwoSampleMR
data.table::setnames(CAD_out_data, 
                     old = c("effect_allele", "other_allele", "beta", "se", "pval", "eaf"),
                     new = c("effect_allele.outcome", "other_allele.outcome", "beta.outcome", 
                             "se.outcome", "pval.outcome", "eaf.outcome"))

# Final check of columns
names(CAD_out_data)
CAD_out_data$id.outcome <- "Coronary artery disease"  # For MR harmonization and output

write.csv(Ctni_exp_data, "Ctni_exp_data.csv", row.names = FALSE)
write.csv(CAD37, "CAD37.csv", row.names = FALSE)
write.csv(CAD_out_data, "CAD_out_data.csv", row.names = FALSE)

# save the files as rDS
saveRDS(Ctni_exp_data, "Ctni_exp_data.rds")
saveRDS(CAD37, "CAD37.rds")
saveRDS(CAD_out_data, "CAD_out_data.rds")

#Harmonization of exposure and outcome data before adding proxy SNPs
datIC <- harmonise_data(
  exposure_dat = Ctni_exp_data,   # cardiac troponin I exposure data
  outcome_dat = CAD_out_data,     # coronary artery disease outcome data
  action = 1                      # automatically resolve strand ambiguities
)

# Identify missing SNPs that were not harmonized
original_snps <- Ctni_exp_data$SNP
harmonised_snps <- datIC$SNP
missing_snps <- setdiff(original_snps, harmonised_snps)
missing_snps  # SNPs that are unmatched and need proxies
#missing_snps"18:36681510(rs151313792)" "15:84825188(rs8024538)"

#Prepare for proxy SNP insertion
Troponin[, MarkerName := sub("(:[^:]*):.*", "\\1", MarkerName)]  # Clean SNP name formatting
setDT(Troponin)

#find proxy SNP
CAD_out_data[CAD_out_data$SNP == "15:84661669", ] # chosen proxy for missing SNP rs8024538
setDT(Troponin)
Troponin[MarkerName == "15:84661669"]

# Remove unmatched SNP from exposure data
Ctni_exp_data <- Ctni_exp_data[Ctni_exp_data$SNP != "15:84825188", ]  # rs8024538 removed
Ctni_exp_data

# Add proxy SNP to exposure data
new_proxy <- data.frame(
  effect_allele.exposure = "A",           
  other_allele.exposure = "G",            
  eaf.exposure = 0.5293,                  
  beta.exposure = 0.0324,                
  se.exposure = 0.0059,                   
  pval.exposure = 3.007e-08,               
  samplesize.exposure = 48115,            
  SNP = "15:84661669", # proxy SNP ID                
  exposure = "cardiac troponin I",
  mr_keep.exposure = TRUE,                
  pval_origin.exposure = "reported",      
  id.exposure = "lsvyno",        
  data_source.exposure = "textfile"   
)

# Add missing columns to match exposure data
new_proxy$r.exposure <- NA
new_proxy$rsq.exposure <- NA
new_proxy$F_statistic <- NA
new_proxy$units.exposure <- NA  

# Match column order and append
newSNP <- new_proxy[, names(Ctni_exp_data)]
Ctni_exp_data <- rbind(Ctni_exp_data, newSNP)

#Recalculate r, Rsq, F-statistics after proxy insertion
Ctni_exp_data$r.exposure <- get_r_from_bsen(
  b = Ctni_exp_data$beta.exposure,
  se = Ctni_exp_data$se.exposure,
  n = Ctni_exp_data$samplesize.exposure
)
Ctni_exp_data$rsq.exposure <- Ctni_exp_data$r.exposure^2
Ctni_exp_data$F_statistic <- with(Ctni_exp_data, (rsq.exposure * (samplesize.exposure - 2)) / (1 - rsq.exposure))

#Combined rsq and combined F
combined_rsq <- sum(Ctni_exp_data$rsq.exposure, na.rm = TRUE)
k <- nrow(Ctni_exp_data)
N <- mean(Ctni_exp_data$samplesize.exposure, na.rm = TRUE)
combined_F <- (combined_rsq / (1 - combined_rsq)) * ((N - k - 1) / k)

#Harmonization after adding proxy SNPs
datIC <- harmonise_data(
  exposure_dat = Ctni_exp_data,
  outcome_dat = CAD_out_data,
  action = 1
)

# Check for palindromic and ambiguous SNPs
palindromic_snps <- datIC %>% filter(palindromic == TRUE)
cat("Palindromic SNPs:\n")
print(palindromic_snps$SNP)

ambiguous_snps <- datIC %>% filter(ambiguous == TRUE)
cat("Ambiguous SNPs:\n")
print(ambiguous_snps$SNP)

# Main MR Analysis
results <- mr(datIC)
print(results)

# Convert to Odds Ratios with 95% Confidence Interval
or_results <- generate_odds_ratios(results)
exp(results$b[1])
exp(results$b[1] - 1.96 * results$se[1])
exp(results$b[1] + 1.96 * results$se[1])

# Penalized Weighted Median
result_pwm <- mr_penalised_weighted_median(
  b_exp = datIC$beta.exposure,
  b_out = datIC$beta.outcome,
  se_exp = datIC$se.exposure,
  se_out = datIC$se.outcome,
  parameters = list(penk = 20, nboot = 1000)
)
print(result_pwm)

# MR-Egger intercept test
egg.int <- mr_pleiotropy_test(datIC)
print(egg.int)

# Read the PheWAS data from a CSV file
phewas_data <- fread("C:/Users/hp/Desktop/phewas_data.csv")

# Adding the -log10(P-value) column
phewas_data <- phewas_data %>% mutate(logp = -log10(`P-value`))

# Arrange the data by increasing logp and convert Phenotype to a factor with levels ordered accordingly
phewas_data <- phewas_data %>%
  arrange(logp) %>%
  mutate(Phenotype = factor(Phenotype, levels = unique(Phenotype)))

# Create a horizontal Manhattan-style bar plot:
# - x axis: Phenotype names (after coordinate flip, will be vertical axis)
# - y axis: -log10(p-value)
# - Bars colored by SNP for easy distinction
# Add a red dashed horizontal line to mark the significance threshold (p = 2.8e-5)
p <- ggplot(phewas_data, aes(x = Phenotype, y = logp, fill = SNP)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  geom_hline(yintercept = -log10(2.8e-5), color = "red", linetype = "dashed") +
  labs(title = "PheWAS Plot for 5 SNPs",
       x = "Phenotype",
       y = "-log10(P-value)",
       fill = "SNP") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# Display the plot in the R graphics window
print(p)

# Save the plot as a PNG file with specified dimensions
ggsave("phewas_plot_5snps_combined.png", plot = p, width = 10, height = 8)

# Heterogeneity test (Cochran’s Q)
het <- mr_heterogeneity(datIC, method_list = "mr_ivw")
print(het)

# Sensitivity analysis
res_single <- mr_singlesnp(datIC, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
print(res_single)

# Rsq for Outcome + Steiger Directionality
CAD_out_data$samplesize.outcome <- 184305
CAD_out_data$r.outcome <- get_r_from_bsen(
  b = CAD_out_data$beta.outcome,
  se = CAD_out_data$se.outcome,
  n = CAD_out_data$samplesize.outcome
)
CAD_out_data$rsq.outcome <- CAD_out_data$r.outcome^2

# Harmonize again with outcome r and Rsq
datIC <- harmonise_data(
  exposure_dat = Ctni_exp_data,
  outcome_dat = CAD_out_data,
  action = 1
)

# Steiger test: directionality (Ctni to CAD)
directionality_test(datIC)

# Manual advanced Steiger test
steiger_results <- mr_steiger2(
  r_exp = datIC$r.exposure,
  r_out = datIC$r.outcome,
  n_exp = datIC$samplesize.exposure[1],
  n_out = datIC$samplesize.outcome[1]
)
print(steiger_results$correct_causal_direction)
print(steiger_results$steiger_test)
print(steiger_results$sensitivity_plot)           

# Steiger filtering: keep only SNPs where exposure explains more variance
dat_steiger_manual <- datIC[datIC$rsq.exposure > datIC$rsq.outcome, ]

# P-values for Steiger filtering
calc_steiger_pval <- function(r_exp, r_out, n_exp) {
  z <- (r_exp^2 - r_out^2) / sqrt(2 / (n_exp - 3))
  p <- 2 * pnorm(-abs(z))
  return(p)
}
dat_steiger_manual$steiger_pval <- mapply(
  calc_steiger_pval,
  r_exp = sqrt(dat_steiger_manual$rsq.exposure),
  r_out = sqrt(dat_steiger_manual$rsq.outcome),
  n_exp = dat_steiger_manual$samplesize.exposure
)
dat_steiger_manual$steiger_dir <- TRUE

# Visualization
# Scatter plot
scatter <- mr_scatter_plot(results, datIC)
print(scatter)
png("scatter.png")
dev.off()
#Scatter plot with labeling of SNPs
scatter <- mr_scatter_plot(results, datIC)
p <- scatter[[1]]  
p <- p + geom_text(
  aes(label = SNP),
  hjust = -0.1,
  vjust = 0.5,
  size = 3,
  data = p$data  
)
print(p)
dev.off()

# Funnel plot
p2 <- mr_funnel_plot(res_single)
print(p2)
png("funnel.png")
print(p2)
dev.off()

# Forest plot (effect size per SNP)
p3 <- mr_forest_plot(res_single)
print(p3)

# Leave-one-out sensitivity analysis
leaveoneout_results <- mr_leaveoneout(datIC)
mr_leaveoneout_plot(leaveoneout_results)

#MR Report Export
output_dir <- "."  # Use the current working directory 

mr_report(
  dat = datIC,
  output_path = output_dir,
  output_type = "pdf",  # or "html"
  author = "RUBA ABDO",
  study = "Cardiac Troponin I ON CAD MR study"
)

#=============================================================================================================#
# DCM Outcome #2

# Load CtnI
CtnI <- read.csv("C:/Users/hp/Desktop/last last last analysis/الجداول النهائية/Ctni_exp_data.csv", 
                 fileEncoding = "UTF-8")  

# Load raw DCM data (with hg19)
file_path <- "C:/Users/hp/Desktop/Files/MY THESIS/letter and terms and papers/DATAAAA NEW/اخر نسخة تلزمني/Allchr DCM.txt"
DCM37 <- fread(file_path)

# Load DCM data (After change build from liftover to genomic build 38)
file_path <- "C:/Users/hp/Desktop/اخر ملفات للرسالة/DCM_saved.csv"
dcm_data <- fread(file_path)

# Separate the frequency range from the FREQ column
dcm_data <- dcm_data %>%
  mutate(
    FREQ_RANGE = str_extract(FREQ, "\\(.*?\\)"),               # << Extract range in parentheses
    freq = str_remove(str_extract(FREQ, "^\\.?\\d+"), "\\.")   # << Extract main frequency value
  )

# Convert extracted frequency column to numeric
dcm_data$freq <- as.numeric(dcm_data$freq)  # << Convert freq to numeric for analysis

# Rename original FREQ column to FREQR for clarity
dcm_data <- dcm_data %>%
  rename(
    FREQR = FREQ  # << Keep original frequency data under new name
  )

# Save modified DCM dataset
file_path <- "C:/Users/hp/Desktop/last last last analysis/dcm_data_MR.csv"
write.csv(dcm_data, file = file_path, row.names = FALSE)  # << Export cleaned data

# Set working directory to project folder
setwd("C:/Users/hp/Desktop/last last last analysis")

# Load the saved file for further checks
data_check <- read.csv("dcm_data_MR.csv")  # << Confirm correct saving
colnames(data_check)  # << Display column names
names(data_check)[names(data_check) == "freq"] <- "FREQ"  # << Rename freq to FREQ for MR package compatibility


# Convert frequency to decimal format (if given as large integer)
data_check$FREQ <- data_check$FREQ / 10000  # << Normalize allele frequency values

#Save final cleaned version for MR input
write.csv(data_check, "dcm_data_fixed MR.csv", row.names = FALSE)

# Load DCM outcome dataset using TwoSampleMR function
DCM_out_data <- read_outcome_data(
  filename = "dcm_data_fixed MR.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  pval_col = "pval",
  eaf_col = "FREQ",
  effect_allele_col = "effect_allele",
  other_allele_col = "noneffect_allele"
)
DCM_out_data$outcome <- "Dilated cardiomyopathy" # << Label outcome for MR plot and results

# Clean DCM outcome data
# Remove rows with missing SNP IDs (important for harmonization)
DCM_out_data_cleaned <- DCM_out_data %>% filter(!is.na(SNP)) 

# Remove duplicate SNPs
DCM_out_data_cleaned <- DCM_out_data_cleaned %>% distinct(SNP, .keep_all = TRUE)

# Save cleaned dataset
write.csv(DCM_out_data_cleaned, "DCM_cleaned.csv", row.names = FALSE)

#Rename data
DCM_out_data <- DCM_out_data_cleaned

# Assign sample size (n) of the DCM GWAS
DCM_out_data$samplesize.outcome <- 8706
Ctni_exp_data <- CtnI
saveRDS(DCM_out_data, "DCM_out_data.rds")

# Harmonise the ctni (exposure) and DCM (outcome) datasets for MR
datID <- harmonise_data(
  exposure_dat = Ctni_exp_data,   # Exposure data: cardiac troponin I
  outcome_dat = DCM_out_data,     # Outcome data: dilated cardiomyopathy
  action = 1                      # Automatically align alleles and remove ambiguous/palindromic SNPs
)
head(datID)

# Identify missing SNPs after harmonisation
original_snps <- Ctni_exp_data$SNP
harmonised_snps <- datID$SNP
missing_snps <- setdiff(original_snps, harmonised_snps)  # SNPs removed in harmonisation
head(missing_snps)  #11:22250324/rs7481951

# Edit MarkerName column to retain only SNP ID (before colon)
Troponin[, MarkerName := sub("(:[^:]*):.*", "\\1", MarkerName)]

# Search for proxy SNPs in European populations for a removed SNP
proxies <- LDproxy(
  snp = "rs7481951",
  pop = c("CEU", "TSI", "FIN", "GBR", "IBS"),
  r2d = "r2",
  token = "52950e216937",
  genome_build = "grch38",
  win_size = 500000
)

# Filter proxy SNPs with Rsq ≥ 0.8 (strong LD)
r2_threshold <- 0.8
proxies_filtered <- subset(proxies, as.numeric(R2) >= r2_threshold)
head(proxies_filtered)

# Remove the SNP that was excluded (rs7481951 /11:22250324)
Ctni_exp_data <- Ctni_exp_data[Ctni_exp_data$SNP != "11:22250324", ]
Ctni_exp_data

#instead this is proxy rs10833712		
#find proxy SNP (RSQ = 0.8382 )
# Check in outcome and Troponin data 
DCM_out_data[DCM_out_data$SNP == "11:22179461", ]
setDT(Troponin)
Troponin[MarkerName == "11:22179461"]

#Add proxy SNPs
new_proxy <- data.frame(
  effect_allele.exposure = "T",
  other_allele.exposure = "G",
  eaf.exposure = 0.4244,
  beta.exposure = -0.0503,
  se.exposure = 0.0059,
  pval.exposure = 1.232e-17,
  samplesize.exposure = 48115,
  SNP = "11:22179461",
  exposure = "cardiac troponin I",
  mr_keep.exposure = TRUE,
  pval_origin.exposure = "reported",
  id.exposure = "ilFJIx",
  data_source.exposure = "textfile"
)

# Add missing columns
new_proxy$units.exposure <- NA
new_proxy$r.exposure <- NA
new_proxy$rsq.exposure <- NA
new_proxy$F_statistic <- NA

# Match structure to original dataset
new_proxy <- new_proxy[, names(Ctni_exp_data)]
Ctni_exp_data <- rbind(Ctni_exp_data, new_proxy)
tail(Ctni_exp_data)

# Identify palindromic_snps and ambiguous
palindromic_snps <- datID %>% filter(palindromic == TRUE)
cat("Palindromic SNPs:\n")
print(palindromic_snps$SNP)

ambiguous_snps <- datID %>% filter(ambiguous == TRUE)
cat("Ambiguous SNPs:\n")
print(ambiguous_snps$SNP)

# Save updated exposure data (with proxies added)
write.csv(Ctni_exp_data, 
          file = "C:/Users/hp/Desktop/last last last analysis/Ctni_exp_dataproxy.csv", 
          row.names = FALSE)
#Recalculate r, Rsq, F-statistics after proxy insertion
Ctni_exp_data$r.exposure <- get_r_from_bsen(
  b = Ctni_exp_data$beta.exposure,
  se = Ctni_exp_data$se.exposure,
  n = Ctni_exp_data$samplesize.exposure
)
Ctni_exp_data$rsq.exposure <- Ctni_exp_data$r.exposure^2
Ctni_exp_data$F_statistic <- with(Ctni_exp_data, (rsq.exposure * (samplesize.exposure - 2)) / (1 - rsq.exposure))

#Combined rsq and combined F statistic
combined_rsq <- sum(Ctni_exp_data$rsq.exposure, na.rm = TRUE)
k <- nrow(Ctni_exp_data)  # number of SNPs (IVs)
N <- mean(Ctni_exp_data$samplesize.exposure, na.rm = TRUE)  # or use consistent N if all same
combined_F <- (combined_rsq / (1 - combined_rsq)) * ((N - k - 1) / k)

#Harmonization after adding proxy SNPs
datID <- harmonise_data(
  exposure_dat = Ctni_exp_data,
  outcome_dat = DCM_out_data,
  action = 1
)


# Perform main MR analysis
results <- mr(datID)
print(results)

# Calculate OR and 95% CI
results$OR <- exp(results$b)
results$CI_lower <- exp(results$b - 1.96 * results$se)
results$CI_upper <- exp(results$b + 1.96 * results$se)

# Penalised Weighted Median method
result_pwm <- mr_penalised_weighted_median(
  b_exp = datID$beta.exposure,
  b_out = datID$beta.outcome,
  se_exp = datID$se.exposure,
  se_out = datID$se.outcome,
  parameters = list(penk = 20, nboot = 1000)
)
print(result_pwm)

# Pleiotropy (Egger intercept test)
egg.int <- mr_pleiotropy_test(datID)
print(egg.int)

# Heterogeneity (Cochran’s Q test)
het <- mr_heterogeneity(datID, method_list = "mr_ivw")
print(het)

# Single SNP effect analysis
res_single <- mr_singlesnp(datID, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
print(res_single)

# Calculate R and Rsq for outcome
DCM_out_data$r.outcome <- get_r_from_bsen(
  b = DCM_out_data$beta.outcome,
  se = DCM_out_data$se.outcome,
  n = DCM_out_data$samplesize.outcome
)
DCM_out_data$rsq.outcome <- DCM_out_data$r.outcome^2
head(DCM_out_data[, c("SNP", "beta.outcome", "se.outcome", "samplesize.outcome", "r.outcome", "rsq.outcome")])

datID <- harmonise_data(
  exposure_dat = Ctni_exp_data,   # Exposure data: cardiac troponin I
  outcome_dat = DCM_out_data,     # Outcome data: dilated cardiomyopathy
  action = 1                      # Automatically align alleles and remove ambiguous/palindromic SNPs
)
head(datID)

# Steiger directionality test (confirm causal direction)
directionality_test(datID)

# Advanced Steiger test using MR-Steiger2
r_exp <- datID$r.exposure
r_out <- datID$r.outcome
n_exp <- datID$samplesize.exposure[1]
n_out <- datID$samplesize.outcome[1]

steiger_results <- mr_steiger2(
  r_exp = r_exp,
  r_out = r_out,
  n_exp = n_exp,
  n_out = n_out
)
print(steiger_results$correct_causal_direction)
print(steiger_results$steiger_test)
print(steiger_results$sensitivity_plot)

# Manual Steiger filtering
# Filter SNPs where exposure Rsq > outcome Rsq
dat_steiger_manual <- datID[datID$rsq.exposure > datID$rsq.outcome, ]

# Calculate Steiger p-values for each SNP
calc_steiger_pval <- function(r_exp, r_out, n_exp) {
  z <- (r_exp^2 - r_out^2) / sqrt(2 / (n_exp - 3))
  p <- 2 * pnorm(-abs(z))
  return(p)
}

dat_steiger_manual$steiger_pval <- mapply(
  calc_steiger_pval,
  r_exp = sqrt(dat_steiger_manual$rsq.exposure),
  r_out = sqrt(dat_steiger_manual$rsq.outcome),
  n_exp = dat_steiger_manual$samplesize.exposure
)

dat_steiger_manual$steiger_dir <- TRUE
head(dat_steiger_manual[, c("rsq.exposure", "rsq.outcome", "steiger_dir", "steiger_pval")])

#Visualization
# Scatter plot 
scatter <- mr_scatter_plot(results, datID)
print(scatter)
png("scatter.png")
dev.off()

# Label SNPs on scatter plot
scatter <- mr_scatter_plot(results, datID)
p <- scatter[[1]]  

p <- p + 
  geom_text_repel(
    aes(label = SNP),
    size = 3,
    max.overlaps = Inf, 
    segment.color = "grey50"
  )

print(p)

# Funnel plot
p2 <- mr_funnel_plot(res_single)
print(p2)
png("funnel.png")
print(p2)
dev.off()

# Forest plot
p3 <- mr_forest_plot(res_single)
print(p3)

# Leave-one-out analysis to assess single SNP influence
leaveoneout_results <- mr_leaveoneout(datID, parameters = default_parameters(), method = mr_ivw)
mr_leaveoneout_plot(leaveoneout_results)

#MR Report Export
output_dir <- "."  # Use the current working directory 

mr_report(
  dat = datID,
  output_path = output_dir,
  output_type = "pdf",  # or "html"
  author = "RUBA ABDO",
  study = "Cardiac Troponin I ON DCM MR study"
)
