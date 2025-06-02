# Title:Cardiac Troponin I and Risk of Coronary Artery Disease and Dilated Cardiomyopathy - a Mendelian Randomization Study
# Master's Thesis Script – RUBA ABDO

# Install packages
install.packages('stringi')
install.packages("devtools")
install.packages("plyr")
install.packages("ggplot2")
install.packages("LDlinkR") 
install.packages("tidyr")
install.packages('mr.raps', repos = c('https://mrcieu.r-universe.dev', 'https://cloud.r-project.org'))

# Load library
library(TwoSampleMR)
library(plyr)
library(ggplot2)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(stringr)
library(dplyr)
library(devtools)
library(LDlinkR)
library(RadialMR)
library(tidyr)
library(stringr)
library(ggrepel)
install_github("MRCIEU/TwoSampleMR")

# Set working directory
setwd("C:/Users/hp/Desktop/اخر ملفات للرسالة")

# Load ctni (exposure) data
ctni_exp_data <- read_exposure_data(
  filename = "ctni_saved.csv",
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
ctni_exp_data$exposure <- "cardiac troponin i"

# units.exposure "SD (inverse rank-normalized)"
ctni_exp_data$units.exposure <- "SD (inverse rank-normalized)"

# Get r
ctni_exp_data$r.exposure <- get_r_from_bsen(
  b = ctni_exp_data$beta.exposure,
  se = ctni_exp_data$se.exposure,
  n = ctni_exp_data$samplesize.exposure
)

# Convert to R-squared
ctni_exp_data$rsq.exposure <- ctni_exp_data$r.exposure^2

#  F for each SNP
ctni_exp_data$F_statistic <- with(ctni_exp_data, (rsq.exposure * (samplesize.exposure - 2)) / (1 - rsq.exposure))
summary(ctni_exp_data$rsq.exposure)
hist(ctni_exp_data$rsq.exposure)
hist(ctni_exp_data$F_statistic)

#Load Cleaned Outcome Data (CAD)
file_path <- "C:/Users/hp/Desktop/CAD/cad.add.160614.website.txt"
CAD37 <- fread(file_path)
CAD37[, chr_chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]

#change the genomic build from hg19 to hg38
gr_hg19 <- GRanges(seqnames = CAD37$chr_chr, ranges = IRanges(start = CAD37$bp_hg19, width = 1))
chain <- import.chain("C:/Users/hp/Desktop/hg19ToHg38.over.chain")
# liftOver
lifted <- liftOver(gr_hg19, chain)
# Get converted positions
chr_hg38 <- sapply(lifted, function(x) if (length(x) > 0) as.character(seqnames(x)) else NA)
pos_hg38 <- sapply(lifted, function(x) if (length(x) > 0) start(x) else NA)
CAD37[, chr_hg38 := gsub("chr", "", chr_hg38)]
CAD37[, bp_hg38 := pos_hg38]
fwrite(CAD37, "C:/Users/hp/Desktop/CAD37_lifted_hg38.txt", sep = "\t")
head(CAD37[, .(markername, chr, bp_hg19, chr_hg38, bp_hg38)])

# Load CHD (outcome) data
CAD_out_data <- read_outcome_data(filename="CAD38_saved.csv", sep=",", snp_col="SNP", beta_col="beta", 
                                  se_col="SE", pval_col="p.value", eaf_col="FREQ", effect_allele_col="effect_allele", other_allele_col="noneffect_allele")    
CAD_out_data$outcome <- "Coronary artery disease"

# Look at the data and check the dimensions to ensure all SNPs are present
head(CAD_out_data)
dim(CAD_out_data)

#unit of outcome
CAD_out_data$units.outcome <- "SD (inverse rank-normalized)"

# Remove raw with SNPs is NA
CAD38_subset_cleaned <- CAD_out_data %>% filter(!is.na(SNP))

# remove duplicated SNPs
CAD38_subset_cleaned <- CAD38_subset_cleaned %>% distinct(SNP, .keep_all = TRUE)
write.csv(CAD38_subset_cleaned, "CAD38_cleaned.csv", row.names = FALSE)
CAD_out_data <-CAD38_subset_cleaned

# Harmonise the ctni and CAD datasets so that the effect alleles are the same.
datCC <- harmonise_data(
  exposure_dat = ctni_exp_data,
  outcome_dat = CAD_out_data,
  action = 1
)

head(datCC)

original_snps <- ctni_exp_data$SNP
harmonised_snps <- datCC$SNP
missing_snps <- setdiff(original_snps, harmonised_snps)

#Add Proxy SNPs
r2_threshold <- 0.8
proxies <- LDproxy(
  snp = "rs8024538",
  pop = c("CEU", "TSI", "FIN", "GBR", "IBS"),
  r2d = "r2",
  token = "52950e216937",
  genome_build = "grch38",
  win_size = 500000
)
proxies_filtered <- subset(proxies, as.numeric(R2) >= r2_threshold)
head(proxies_filtered)

CAD_out_data[CAD_out_data$SNP == "15:84661669:A:G", ]
setDT(troponin)
troponin[MarkerName == "15:84661669:A:G"]

# removed SNP to add proxies ones 
ctni_exp_data <- ctni_exp_data[ctni_exp_data$SNP != "11:22205338", ]
ctni_exp_data

# Insert proxy snp in exposure data
new_proxy_2 <- data.frame(
  effect_allele.exposure = "T",           
  other_allele.exposure = "C",            
  eaf.exposure = 0.4592,                  
  beta.exposure = -0.0507,                
  se.exposure = 0.0059,                   
  pval.exposure = 1.23e-17,               
  samplesize.exposure = 48115,            
  SNP = "11:22205338",                    
  exposure = "cardiac troponin i",
  mr_keep.exposure = TRUE,                
  pval_origin.exposure = "reported",      
  id.exposure = "vqwpb2",        
  data_source.exposure = "textfile"   
)
newSNP <- new_proxy_2[, names(ctni_exp_data)]
ctni_exp_data <- rbind(ctni_exp_data, newSNP)

# check palindromic and Ambiguous SNPs
palindromic_snps <- datCC %>% filter(palindromic == TRUE)
cat("Palindromic SNPs:\n")
print(palindromic_snps$SNP)

ambiguous_snps <- datCC %>% filter(ambiguous == TRUE)
cat("Ambiguous SNPs:\n")
print(ambiguous_snps$SNP)

# MR analysis
results <- mr(datCC)
res <- mr(datCC, method_list = c("mr_ivw", "mr_simple_median", "mr_weighted_median", "mr_egger_regression"))
print(results)

#OR and CI
generate_odds_ratios(results)
or_results <- generate_odds_ratios(results)
exp(res$b[1])
exp(res$b[1] - 1.96 * res$se[1])
exp(res$b[1] + 1.96 * res$se[1])

# Sensitivity analysis
# Single SNP analysis
res_single <- mr_singlesnp(datCC, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
print(res_single)

# mr_heterogeneity-IVW
het <- mr_heterogeneity(datCC, method_list = "mr_ivw")
print(het)

# MR-Egger
egg.int <- mr_pleiotropy_test(datCC)
print(egg.int)

# Run MR Rucker framework on harmonized data
result_rucker <- mr_rucker(datCC)
print(result_rucker)

#RadialMR
radial_dat <- dat_to_RadialMR(datCC)

#RadialMR_IVW
result_ivw <- ivw_radial(radial_dat[[1]], alpha = 0.05, weights = 3)
print(result_ivw)

# Run MR-PRESSO on your harmonised data
mr_presso_results <- run_mr_presso(datCC, NbDistribution = 1000, SignifThreshold = 0.05)
print(mr_presso_results)

# Penalised Weighted Median MR
b_exp <- datCC$beta.exposure
b_out <- datCC$beta.outcome
se_exp <- datCC$se.exposure
se_out <- datCC$se.outcome
result_pwm <- mr_penalised_weighted_median(
  b_exp = b_exp,
  b_out = b_out,
  se_exp = se_exp,
  se_out = se_out,
  parameters = list(penk = 20, nboot = 1000)  # optional, adjust as needed
)
print(result_pwm)

#Maximum likelihood MR method
result_ml <- mr_two_sample_ml(
  b_exp = datCC$beta.exposure,
  b_out = datCC$beta.outcome,
  se_exp = datCC$se.exposure,
  se_out = datCC$se.outcome,
  parameters = list()  # optional: can use default settings
)
print(result_ml)
#mr_raps robust test
# Extract vectors from merged dataset
b_exp <- datCC$beta.exposure
se_exp <- datCC$se.exposure
b_out <- datCC$beta.outcome
se_out <- datCC$se.outcome

# Run MR-RAPS with default parameters (robust loss, overdispersion = TRUE)
result_raps <- mr_raps(
  b_exp = b_exp,
  b_out = b_out,
  se_exp = se_exp,
  se_out = se_out
)
print(result_raps)

#Estimate r vaule for outcome
CAD_out_data$samplesize.outcome <- 184305
CAD_out_data$r.outcome <- get_r_from_bsen(
  b = CAD_out_data$beta.outcome,
  se = CAD_out_data$se.outcome,
  n = CAD_out_data$samplesize.outcome
)
# Convert to R-squared
CAD_out_data$rsq.outcome <- CAD_out_data$r.outcome^2
head(CAD_out_data[, c("SNP", "beta.outcome", "se.outcome", "samplesize.outcome", "r.outcome", "rsq.outcome")])

# Harmonise again after adding r valuse of outcome and rsquare
datCC <- harmonise_data(
  exposure_dat = ctni_exp_data,
  outcome_dat = CAD_out_data,
  action = 1
)
head(dat)

#MR Steiger directionality test
directionality_test(datCC)

#steiger2 Advanced test
r_exp <- datCC$r.exposure
r_out <- datCC$r.outcome

# Assuming sample sizes are consistent across SNPs, take the first
n_exp <- datCC$samplesize.exposure[1]
n_out <- datCC$samplesize.outcome[1]
steiger_results <- mr_steiger2(
  r_exp = r_exp,
  r_out = r_out,
  n_exp = n_exp,
  n_out = n_out
)
print(steiger_results$correct_causal_direction)    
print(steiger_results$steiger_test)                 
print(steiger_results$sensitivity_plot)             
# steiger filtering
dat_steiger_manual <- datCC[datCC$rsq.exposure > datCC$rsq.outcome, ]

# calculate pvalue for each SNP in steiger
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

# scatter plot 
scatter <- mr_scatter_plot(results, datCC)
print(scatter)
png("scatter.png")
dev.off()
scatter <- mr_scatter_plot(results, datCC)
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

# Forest plot
p3 <- mr_forest_plot(res_single)
print(p3)
p3 <- mr_forest_plot(or_results)
or_results
png("forest.png")
print(p3)
dev.off()
p3

# Forest Plot
ggplot(plot_data, aes(x = reorder(SNP, or), y = or)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  coord_flip() +
  labs(x = "SNP", y = "Odds Ratio (OR)", title = "Forest Plot for OR") +
  theme_minimal()

# leaveoneout test
leaveoneout_results <- mr_leaveoneout(dat, parameters = default_parameters(), method = mr_ivw)
mr_leaveoneout_plot(leaveoneout_results)

#________________________________________________________________________________________#
#DCM Outcome #2

file_path <- "C:/Users/hp/Desktop/Files/MY THESIS/letter and terms and papers/DATAAAA NEW/اخر نسخة تلزمني/Allchr DCM.txt"
DCM37 <- fread(file_path)
file_path <- "C:/Users/hp/Desktop/اخر ملفات للرسالة/DCM_saved.csv"
dcm_data <- fread(file_path)
#to separate FREQ column
dcm_data <- dcm_data %>%
  mutate(
    FREQ_RANGE = str_extract(FREQ, "\\(.*?\\)"),              
    freq = str_remove(str_extract(FREQ, "^\\.?\\d+"), "\\.")  
  )
dcm_data$freq <- as.numeric(dcm_data$freq)
# rename columns
dcm_data <- dcm_data %>%
  rename(
    FREQR = FREQ
  )
dcm_data <- dcm_data %>%
  rename(
    FREQ = freq
  )
desktop_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop")

file_path <- file.path(desktop_path, "dcm_data_MR.csv")

write.csv(dcm_data, file = file_path, row.names = FALSE)

names(dcm_data)
data_check <- read.csv("dcm_data_MR.csv")
colnames(data_check)
head(data_check$FREQ)
str(data_check$FREQ)
summary(data_check$FREQ)
data_check$FREQ <- data_check$FREQ / 10000
summary(data_check$FREQ)
write.csv(data_check, "dcm_data_fixed MR.csv", row.names = FALSE)

# load lifted DCM outcome data
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
DCM_out_data$outcome <- "Dilated cardiomyopathy"


# Removed NA SNPs
DCM_out_data_cleaned <- DCM_out_data %>% filter(!is.na(SNP))

# removed duplicated SNP
DCM_out_data_cleaned <- DCM_out_data_cleaned %>% distinct(SNP, .keep_all = TRUE)
write.csv(DCM_out_data_cleaned, "DCM_cleaned.csv", row.names = FALSE)
DCM_out_data <- DCM_out_data_cleaned
DCM_out_data$samplesize.outcome <- 8706

# Harmonise the ctni and DCM datasets
datCD <- harmonise_data(
  exposure_dat = ctni_exp_data,
  outcome_dat = DCM_out_data,
  action = 1
)
head(dat)
original_snps <- ctni_exp_data$SNP
harmonised_snps <- dat$SNP
missing_snps <- setdiff(original_snps, harmonised_snps)
head(missing_snps)
11:22250324:A:T
#full ctni summary statisitics 
file_path <- "C:/Users/hp/Desktop/اخر ملفات للرسالة/troponin_gwas_summary_statistics.txt"
troponin_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(troponin_data)

# structure
str(troponin_data)
# استخدام LDproxy 
library(ldlinkr)
r2_threshold <- 0.8
proxies <- LDproxy(
  snp = "rs8024538",
  pop = c("CEU", "TSI", "FIN", "GBR", "IBS"),
  r2d = "r2",
  token = "52950e216937",
  genome_build = "grch38",
  win_size = 500000
)
proxies_filtered <- subset(proxies, as.numeric(R2) >= r2_threshold)
head(proxies_filtered)
#investigate if the proxy SNP is in the data
DCM_out_data[DCM_out_data$SNP == "11:22205338", ]
setDT(troponin_data)
troponin_data[MarkerName == "11:22205338:T:C"]

str(troponin_data$MarkerName)
"11:22205338:T:C" %in% troponin_data$MarkerName

# delete raw related to rs7481951 in ctni_exp_data to add proxy snp
ctni_exp_data <- ctni_exp_data[ctni_exp_data$SNP != "11:22250324", ]

#Add proxy SNPs
new_proxy <- data.frame(
  effect_allele.exposure = "T",
  other_allele.exposure = "C",
  eaf.exposure = 0.4592,
  beta.exposure = -0.0507,
  se.exposure = 0.0059,
  pval.exposure = 1.23e-17,
  samplesize.exposure = 48115,
  SNP = "11:22205338",
  exposure = "cardiac troponin i",
  mr_keep.exposure = TRUE,
  pval_origin.exposure = "reported",
  id.exposure = "SgRqRQ",
  data_source.exposure = "textfile"
)
newSNP <- new_proxy[, names(ctni_exp_data)]
ctni_exp_data <- rbind(ctni_exp_data, newSNP)
names(ctni_exp_data)
names(new_proxy)
# make some editings
new_proxy$units.exposure <- NA
new_proxy$r.exposure <- NA
new_proxy$rsq.exposure <- NA
new_proxy$F_statistic <- NA

new_proxy <- new_proxy[, names(ctni_exp_data)]
ctni_exp_data <- rbind(ctni_exp_data, new_proxy)
tail(ctni_exp_data)

#palindromic_snps and ambiguous
palindromic_snps <- dat %>% filter(palindromic == TRUE)
cat("Palindromic SNPs:\n")
print(palindromic_snps$SNP)
ambiguous_snps <- dat %>% filter(ambiguous == TRUE)
cat("Ambiguous SNPs:\n")
print(ambiguous_snps$SNP)
write.csv(ctni_exp_data, 
          file = "C:/Users/hp/Desktop/tables otf thesis/ctni_exp_data.csv", 
          row.names = FALSE)

# Harmonise the ctni and DCM datasets
datCD <- harmonise_data(
  exposure_dat = ctni_exp_data,
  outcome_dat = DCM_out_data,
  action = 1
)
head(datCD)
# MR analysis
results <- mr(datCD)
res <- mr(datCD, method_list = c("mr_ivw", "mr_simple_median", "mr_weighted_median", "mr_egger_regression"))
print(results)
write.csv(results, 
          file = "C:/Users/hp/Desktop/tables otf thesis/results.csv", 
          row.names = FALSE)

# OR and CI
or_results <- generate_odds_ratios(results)
results$OR <- exp(results$b)
results$CI_lower <- exp(results$b - 1.96 * results$se)
results$CI_upper <- exp(results$b + 1.96 * results$se)
write.csv(or_results, 
          file = "C:/Users/hp/Desktop/tables otf thesis/results.csv", 
          row.names = FALSE)
# Single SNP analysis
res_single <- mr_singlesnp(datCD, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
print(res_single)
write.csv(res_single, 
          file = "C:/Users/hp/Desktop/tables otf thesis/res_single.csv", 
          row.names = FALSE)

# Heterogeneity & Pleiotropy
het <- mr_heterogeneity(datCD, method_list = "mr_ivw")
print(het)

egg.int <- mr_pleiotropy_test(datCD)
print(egg.int)

# MR Rucker framework
result_rucker <- mr_rucker(datCD)
print(result_rucker)

# Radial MR
radial_dat <- dat_to_RadialMR(datCD)
result_ivw <- ivw_radial(radial_dat[[1]], alpha = 0.05, weights = 3)
print(result_ivw)

# MR-PRESSO
mr_presso_results <- run_mr_presso(datCD, NbDistribution = 1000, SignifThreshold = 0.05)
print(mr_presso_results)

# Penalised Weighted Median
result_pwm <- mr_penalised_weighted_median(
  b_exp = datCD$beta.exposure,
  b_out = datCD$beta.outcome,
  se_exp = datCD$se.exposure,
  se_out = datCD$se.outcome,
  parameters = list(penk = 20, nboot = 1000)
)
print(result_pwm)

# Maximum likelihood method
result_ml <- mr_two_sample_ml(
  b_exp = datCD$beta.exposure,
  b_out = datCD$beta.outcome,
  se_exp = datCD$se.exposure,
  se_out = datCD$se.outcome,
  parameters = list()
)
print(result_ml)

# MR-RAPS
result_raps <- mr_raps(
  b_exp = datCD$beta.exposure,
  b_out = datCD$beta.outcome,
  se_exp = datCD$se.exposure,
  se_out = datCD$se.outcome
)
print(result_raps)

# R and R-squared of outcome
DCM_out_data$r.outcome <- get_r_from_bsen(
  b = DCM_out_data$beta.outcome,
  se = DCM_out_data$se.outcome,
  n = DCM_out_data$samplesize.outcome
)
DCM_out_data$rsq.outcome <- DCM_out_data$r.outcome^2
head(DCM_out_data[, c("SNP", "beta.outcome", "se.outcome", "samplesize.outcome", "r.outcome", "rsq.outcome")])

# Steiger directionality test
directionality_test(datCD)

# steiger2 analysis Advanced analysis
r_exp <- datCD$r.exposure
r_out <- datCD$r.outcome
n_exp <- datCD$samplesize.exposure[1]
n_out <- datCD$samplesize.outcome[1]

steiger_results <- mr_steiger2(
  r_exp = r_exp,
  r_out = r_out,
  n_exp = n_exp,
  n_out = n_out
)
print(steiger_results$correct_causal_direction)
print(steiger_results$steiger_test)
print(steiger_results$sensitivity_plot)

# steiger filtering

summary(datCD$rsq.exposure)
summary(datCD$rsq.outcome)
table(datCD$rsq.exposure > datCD$rsq.outcome)
install.packages("remotes")  # لو ما عندك
remotes::install_github("MRCIEU/TwoSampleMR")
#filtering manually 
dat_steiger_manual <- datCD[datCD$rsq.exposure > datCD$rsq.outcome, ]

#  p-value of Steiger to each SNP
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
# scatter plot 
scatter <- mr_scatter_plot(results, datCD)
print(scatter)
png("scatter.png")
dev.off()

#Label SNPs
scatter <- mr_scatter_plot(results, datCD)
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

#OR as numeric
or_results$or <- as.numeric(or_results$or)
or_results$or_lci95 <- as.numeric(or_results$or_lci95)
or_results$or_uci95 <- as.numeric(or_results$or_uci95)

# Forest Plot 
ggplot(or_results, aes(x = method, y = or)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = or_lci95, ymax = or_uci95), width = 0.2) +
  coord_flip() +
  labs(x = "Method", y = "Odds Ratio (OR)", title = "Forest Plot by Method") +
  theme_minimal()

# leaveoneout test
leaveoneout_results <- mr_leaveoneout(datCD, parameters = default_parameters(), method = mr_ivw)
mr_leaveoneout_plot(leaveoneout_results)

#_________________________________________________________________________
#Bidirectional MR

#load exposure data
CADEX <- read.csv("C:/Users/hp/Desktop/اخر ملفات للرسالة/CAD SNPS.csv", fileEncoding = "UTF-8")
#filter SNPs based on pvalue
cadexsignificant_snps <- subset(CADEX, pvalue < 5e-8)
#editings on data
cadexsignificant_snps <- cadexsignificant_snps %>%
  separate(Alleles, into = c("effect_allele", "noneffect_allele"), sep = "/")
names(cadexsignificant_snps)[names(cadexsignificant_snps) == "lead.OR..95..CI."] <- "OR AND CI"
cadexsignificant_snps$OR_num <- as.numeric(str_extract(cadexsignificant_snps$OR, "^[0-9.]+"))
cadexsignificant_snps$`95%CI` <- str_extract(cadexsignificant_snps$OR, "\\(.*\\)")
#change the name of columns
names(cadexsignificant_snps)[names(cadexsignificant_snps) == "OR"] <- "OR_CI"
names(cadexsignificant_snps)[names(cadexsignificant_snps) == "OR_num"] <- "OR"
cadex <- cadexsignificant_snps
#save data
write.csv(cadex, "C:/Users/hp/Desktop/cadex.csv", row.names = FALSE)
#load saved data
cade38 <- read.csv("C:/Users/hp/Desktop/cadex.csv", fileEncoding = "UTF-8")
cade38$beta <- log(as.numeric(cade38$OR))
names(cade38)
cade38 <- as.data.table(cade38)
setnames(cade38, old = c("beta", "effect_allele", "noneffect_allele", "FREQ", "pvalue"),
         new = c("beta.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "pval.exposure"))
cade38[, exposure := "CAD"]
names(cade38)
cade38$samplesize.exposure <-184305
for (i in seq_along(cade38)) {
  if (any(grepl("^[0-9XY]+:[0-9]+$", cade38[[i]]))) {
    names(cade38)[i] <- "CHR:POS"
  }
}
if (all(c("CHR:POS", "effect_allele.exposure", "other_allele.exposure") %in% names(cade38))) {
  
  cade38$SNP <- paste(
    cade38$`CHR:POS`,
    cade38$effect_allele.exposure,
    cade38$other_allele.exposure,
    sep = ":"
  )
}

CAD_exp_datavf <- read.csv("C:/Users/hp/Desktop/CAD_exp_datavf.csv")
CAD_exp_datavf$r.exposure <- NULL
CAD_exp_datavf$F_statistic <- NULL
CAD_exp_datavf$CHR.POS <- NULL
CAD_exp_datavf$SE <- NULL
CAD_exp_datavf$MarkerName <- NULL
write.csv(CAD_exp_data, "C:/Users/hp/Desktop/CAD_exp_datavf.csv", row.names = FALSE)


# Set working directory
setwd("C:/Users/hp/Desktop/اخر ملفات للرسالة")
CAD_out_data <- read_outcome_data(filename="CAD38_saved.csv", sep=",", snp_col="SNP", beta_col="beta", 
                                  se_col="SE", pval_col="p.value", eaf_col="FREQ", effect_allele_col="effect_allele", other_allele_col="noneffect_allele")    
CAD_out_data$outcome <- "Coronary artery disease"
names(CAD_out_data)[names(CAD_out_data) == "SNP"] <- "CHR:POS"
CAD_out_data$SNP <- paste0(CAD_out_data$`CHR:POS`, ":", CAD_out_data$effect_allele.outcome, ":", CAD_out_data$other_allele.outcome)
CAD_exp_datavf$se.exposure <- CAD_out_data$se.outcome[match(CAD_exp_datavf$SNP, CAD_out_data$SNP)]
cade38$se.exposure <- CAD_out_data$se.outcome[match(cade38$SNP, CAD_out_data$SNP)]
# Get r
CAD_exp_datavf$r.exposure <- get_r_from_bsen(
  b = CAD_exp_datavf$beta.exposure,
  se = CAD_exp_datavf$se.exposure,
  n = CAD_exp_datavf$samplesize.exposure
)
cade38$r.exposure <- get_r_from_bsen(
  b = cade38$beta.exposure,
  se = cade38$se.exposure,
  n = cade38$samplesize.exposure
)
# Convert to R-squared
CAD_exp_datavf$rsq.exposure <- CAD_exp_datavf$r.exposure^2
cade38$rsq.exposure <- cade38$r.exposure^2

#  F of each SNP
CAD_exp_datavf$F_statistic <- with(CAD_exp_datavf, (rsq.exposure * (samplesize.exposure - 2)) / (1 - rsq.exposure))
summary(CAD_exp_datavf$rsq.exposure)
cade38$se.exposure <- 0.01
cade38$F_statistic <- with(cade38, (rsq.exposure * (samplesize.exposure - 2)) / (1 - rsq.exposure))

#units.exposure "SD (inverse rank-normalized)"
CAD_exp_datavf$units.exposure <- "SD (inverse rank-normalized)"

#load outcome data
troponin_data <- read.table("C:/Users/hp/Desktop/اخر ملفات للرسالة/troponin_gwas_summary_statistics.txt", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE, fileEncoding = "UTF-8")
names(troponin_data)
troponin_outcome <- as.data.table(troponin_data)
# Rename columns to match TwoSampleMR outcome expectations
setnames(troponin_outcome,
         old = c("MarkerName","Allele1", "Allele2", "Effect", "StdErr", "P.value","Freq1"),
         new = c("SNP","effect_allele.outcome", "other_allele.outcome", "beta.outcome", "se.outcome", "pval.outcome","eaf.outcome"))
# Add outcome name
troponin_outcome[, outcome := "Cardiac Troponin I"]
troponin_outcome[, id.outcome := "troponin_i"]
# Keep only relevant columns
troponin_outcome <- troponin_outcome[, .(SNP, effect_allele.outcome, other_allele.outcome,
                                         beta.outcome, se.outcome, pval.outcome, outcome, id.outcome ,eaf.outcome)]
head(troponin_outcome)
Ctni_outcome <- troponin_outcome

# harmonisation 
dat <- harmonise_data(
  exposure_dat = CAD_exp_datavf,
  outcome_dat = Ctni_outcome,
  action = 1
)
original_snps <- CAD_exp_datavf$SNP
harmonised_snps <- dat$SNP
missing_snps <- setdiff(original_snps, harmonised_snps)
head(missing_snps)
missing_snps
# palindromic_snps and ambiguous
palindromic_snps <- dat %>% filter(palindromic == TRUE)
cat("Palindromic SNPs:\n")
print(palindromic_snps$SNP)
ambiguous_snps <- dat %>% filter(ambiguous == TRUE)
cat("Ambiguous SNPs:\n")
print(ambiguous_snps$SNP)

#get proxy SNPS instead of removed SNPs 
install.packages("LDlinkR")
library(LDlinkR)
#Set population and r2 threshold
r2_threshold <- 0.8
population <- c("CEU", "TSI", "FIN", "GBR", "IBS")
# Set your API token
token <- "52950e216937" #this is mine
#  list of SNPs
snp_list <- c(
  "rs11206510", "rs6689306", "rs67180937", "rs16986953", "rs515135",
  "rs6544713", "rs7568458", "rs17678683", "rs6725887", "rs9818870",
  "rs4593108", "rs72689147", "rs9349379", "rs56336142", "rs12202017",
  "rs55730499", "rs4252185", "rs2107595", "rs11556924", "rs2891168",
  "rs2519093", "rs2487928", "rs1870634", "rs11838776", "rs1412444",
  "rs11191416", "rs2128739", "rs2681472", "rs3184504", "rs11838776",
  "rs10139550", "rs4468572", "rs56289821", "rs4420638", "rs28451064",
  "rs17087335", "rs3918226", "rs10840293", "rs56062135", "rs8042271",
  "rs7212798", "rs663129", "rs180803"
)

# Remove duplicates from the list
snp_list <- unique(snp_list)
# Create an empty list to store results
ld_results <- list()
# Loop through each SNP
for (snp in snp_list) {
  message(paste("Processing SNP:", snp))
  tryCatch({
    res <- LDproxy(
      snp = snp,
      pop = population,
      r2d = "r2",
      token = token
    )
    
    # Filter results by r² threshold
    res_filtered <- subset(res, as.numeric(R2) >= r2_threshold)
    
    # Save to list
    ld_results[[snp]] <- res_filtered
  }, error = function(e) {
    message(paste("Error with SNP:", snp, "->", e$message))
    ld_results[[snp]] <- NULL
  })
}

# To investigate if they are in the outcome data and to get their information from CAD FULL SUMMARY 
Ctni_outcome[Ctni_outcome$SNP == "22:24283080:T:G", ]
setDT(CAD_out_data)
CAD_out_data[SNP == "22:24283080:T:G"]

# delete SNP that removed after harmonization to add proxy instead
CAD_exp_datavf <- CAD_exp_datavf[CAD_exp_datavf$rsID != "rs142130958", ]
# تحويل السطر الجديد إلى data.table
new_cad_proxy <- data.table(
  Known.locus = "LDLR",
  CHR = 19,
  rsID = "rs142130958",
  effect_allele.exposure = "G",
  other_allele.exposure = "A",
  eaf.exposure = 0.900107,
  pval.exposure = 7.149703e-14,
  SNP = "19:11079976:G:A",
  POS = "11079976",
  `CHR:POS` = "19:11079976",
  beta.exposure = 0.12564700,
  id.exposure = "CAD",
  exposure = "CAD",
  samplesize.exposure = 184305
)

# Convert to data.tables
dt1 <- as.data.table(CAD_exp_datavf)
dt2 <- as.data.table(new_cad_proxy)

# Combine with fill = TRUE to allow different columns
combined_data <- rbindlist(list(dt1, dt2), fill = TRUE)

# Convert back to data.frame 
CAD_exp_datavf <- as.data.frame(combined_data)
write.csv(CAD_exp_data, "C:/Users/hp/Desktop/CAD_exp_datavf.csv", row.names = FALSE)

#repeat harmonization after adding proxies SNPS 
dat <- harmonise_data(
  exposure_dat = CAD_exp_datavf,
  outcome_dat = Ctni_outcome,
  action = 1
)
head(dat)
Ctni_outcome$samplesize.outcome <- 48115

# Mendelian randomization analysis
results <- mr(dat)
res <- mr(dat, method_list = c("mr_ivw", "mr_simple_median", "mr_weighted_median", "mr_egger_regression"))
print(results)
#OR and CI 
or_results <- generate_odds_ratios(results)
res$OR <- exp(res$b)
res$CI_lower <- exp(res$b - 1.96 * res$se)
res$CI_upper <- exp(res$b + 1.96 * res$se)

# Single SNP Analysis
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
print(res_single)

# Heterogeneity & Pleiotropy
het <- mr_heterogeneity(dat, method_list = "mr_ivw")
print(het)

egg.int <- mr_pleiotropy_test(dat)
print(egg.int)

# MR Rucker framework
result_rucker <- mr_rucker(dat)
print(result_rucker)

# Radial MR
radial_dat <- dat_to_RadialMR(dat)
result_ivw <- ivw_radial(radial_dat[[1]], alpha = 0.05, weights = 3)
print(result_ivw)

# MR-PRESSO
mr_presso_results <- run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)
print(mr_presso_results)

# Penalised Weighted Median
result_pwm <- mr_penalised_weighted_median(
  b_exp = dat$beta.exposure,
  b_out = dat$beta.outcome,
  se_exp = dat$se.exposure,
  se_out = dat$se.outcome,
  parameters = list(penk = 20, nboot = 1000)
)
print(result_pwm)

# Maximum likelihood method
result_ml <- mr_two_sample_ml(
  b_exp = dat$beta.exposure,
  b_out = dat$beta.outcome,
  se_exp = dat$se.exposure,
  se_out = dat$se.outcome,
  parameters = list()
)
print(result_ml)

# MR-RAPS
result_raps <- mr_raps(
  b_exp = dat$beta.exposure,
  b_out = dat$beta.outcome,
  se_exp = dat$se.exposure,
  se_out = dat$se.outcome
)
print(result_raps)
Ctni_outcome$samplesize.outcome <- 48115

# R and R-squared of outcome
Ctni_outcome$r.outcome <- get_r_from_bsen(
  b = Ctni_outcome$beta.outcome,
  se = Ctni_outcome$se.outcome,
  n = Ctni_outcome$samplesize.outcome
)
Ctni_outcome$rsq.outcome <- Ctni_outcome$r.outcome^2
head(Ctni_outcome[, c("SNP", "beta.outcome", "se.outcome", "samplesize.outcome", "r.outcome", "rsq.outcome")])

# Steiger directionality test
directionality_test(dat)

# steiger2 analysis addvanced one 
r_exp <- dat$r.exposure
r_out <- dat$r.outcome
n_exp <- dat$samplesize.exposure[1]
n_out <- dat$samplesize.outcome[1]

steiger_results <- mr_steiger2(
  r_exp = r_exp,
  r_out = r_out,
  n_exp = n_exp,
  n_out = n_out
)
print(steiger_results$correct_causal_direction)
print(steiger_results$steiger_test)
print(steiger_results$sensitivity_plot)

# Make sure these exist
summary(dat$rsq.exposure)
summary(dat$rsq.outcome)
table(dat$rsq.exposure > dat$rsq.outcome)
install.packages("remotes")  
remotes::install_github("MRCIEU/TwoSampleMR")

# steiger filtering
dat_steiger_manual <- dat[dat$rsq.exposure > dat$rsq.outcome, ]

# calculate pvalue for each SNP in steiger
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

# visualization
p1 <- mr_scatter_plot(results, dat)[[1]]
png("scatter.png"); print(p1); dev.off()
p1
p2 <- mr_funnel_plot(res_single)[[1]]
png("funnel.png"); print(p2); dev.off()
p2
p3 <- mr_forest_plot(res_single)[[1]]
png("forest.png"); print(p3); dev.off()
p3
# Leave-one-out analysis
leaveoneout_results <- mr_leaveoneout(dat, parameters = default_parameters(), method = mr_ivw)
p4 <- mr_leaveoneout_plot(leaveoneout_results)[[1]]
png("leaveoneout.png"); print(p4); dev.off()
p4
scatter <- mr_scatter_plot(results, datCD)
print(scatter)
png("scatter.png")
dev.off()

#Label SNPs
scatter <- mr_scatter_plot(results, dat)
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

#OR as numeric
or_results$or <- as.numeric(or_results$or)
or_results$or_lci95 <- as.numeric(or_results$or_lci95)
or_results$or_uci95 <- as.numeric(or_results$or_uci95)

# Forest Plot 
ggplot(or_results, aes(x = method, y = or)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = or_lci95, ymax = or_uci95), width = 0.2) +
  coord_flip() +
  labs(x = "Method", y = "Odds Ratio (OR)", title = "Forest Plot by Method") +
  theme_minimal()

# leaveoneout test
leaveoneout_results <- mr_leaveoneout(datCD, parameters = default_parameters(), method = mr_ivw)
mr_leaveoneout_plot(leaveoneout_results)
