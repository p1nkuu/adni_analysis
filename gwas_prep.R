# QUALITY CONTROL DATA

# MAF
frq<-read.table('adni_freq.frq', header=TRUE)
frq_df<-data.frame(frq)
num_variants<-nrow(frq)
maf_0.05<-frq_df[frq_df$MAF < 0.05,]
pct_maf_0.05<-nrow(maf_0.05) / num_variants
# 54% have MAF>0.05 from a total number of variants of 2379855

# MISSING DATA
miss<-data.frame(read.table('missing.lmiss', header=TRUE))
miss_filtered<-miss[miss$F_MISS>0.05,]
pct_miss_variant_0.05<-nrow(miss_filtered)/nrow(miss)
#0.876% of SNPs have missingness

miss_i<-data.frame(read.table('missing.imiss', header=TRUE))
miss_i_filt<-miss_i[miss_i$F_MISS>0.05,]
pct_miss_i_0.05<-nrow(miss_i_filt)/nrow(miss_i)
# None of the individuals have >0.05 missingness

# HARDY WEINBERG
hwe<-data.frame(read.table('all_hwe.hwe',header=TRUE))
hwe_filtered<-hwe[hwe$P<0.000001,]
pct_hwe_0.05<-nrow(hwe_filtered)/nrow(hwe)
# 0.336% have hwe<0.000001

"
1 = Am Indian/Alaskan
     2 = Asian 
     3 = Hawaiian/Other PI
     4 = Black 
     5 = White 
     6 = More than one 
     7 = Unknown
"
# create dummy variables for the PTRACCAT (race) column for plink analysis
merged$am_indian_alaskan<-ifelse(merged$PTRACCAT == 1, 1, 0)
merged$asian<-ifelse(merged$PTRACCAT == 2, 1, 0)
merged$pacific_islander<-ifelse(merged$PTRACCAT == 3, 1, 0)
merged$black<-ifelse(merged$PTRACCAT == 4, 1, 0)
merged$white<-ifelse(merged$PTRACCAT == 5, 1, 0)
merged$more_than_one<-ifelse(merged$PTRACCAT == 6, 1, 0)
merged$unknown<-ifelse(merged$PTRACCAT == 7, 1, 0)

# check distribution of the races
native<-sum(merged$am_indian_alaskan == 1)
asian<-sum(merged$asian == 1)
pacific_islander<-sum(merged$pacific_islander == 1)
black<-sum(merged$black == 1)
white<-sum(merged$white == 1)
more_than_one<-sum(merged$more_than_one == 1)
unknown<-sum(merged$unknown == 1)
race_data<-data.frame(
  Race = c("AI", "Asian", "PI", "Black", "White", "MTO", "Unknown"),
  Number_of_Subjects = c(native, asian, pacific_islander, black, white, more_than_one, unknown)
)
barplot(race_data$Number_of_Subjects, names=race_data$Race, ylab = "Number of Participants", main ="Race Distribution of ADNI dataset")
print(summary(merged, select=c(PTRACCAT)))

# check distribution of females and males
female_percent<-sum(merged$female == 1) / 812

# check distribution of APOE4 allele count
no_apoe4<-sum(merged$APOE4 == 0) / 812
one_apoe4<-sum(merged$APOE4 == 1) / 812
two_apoe4<-sum(merged$APOE4 == 2) / 812

# PHENOTYPE ANALYSIS
library("tidyverse")
library("haven")
pheno<-data.frame(read_sas("adni_pheno_2023july28.sas7bdat"))
fam<-data.frame(read.table('WGS_Omni25_BIN_wo_ConsentsIssues.fam', header=FALSE))
# new data with those with genotype data
merged<-data.frame(merge(pheno, fam, by.fam="V2")) # oops didn't do this right
merged<-subset(merged, select = -c(V1, V2, V3, V4, V5, V6))
merged<-merged[!duplicated(merged), ]
# Plots of the MRI Scan data
mri_col<-subset(merged, select = c(Hippocampus_bl, WholeBrain_bl, Entorhinal_bl))
mri_clean<-mri_col[complete.cases(mri_col[, 1:3]), ]
hist(mri_clean$Hippocampus_bl)
hist(mri_clean$WholeBrain_bl)
hist(mri_clean$Entorhinal_bl)

merged_summary<-summary(mri_clean)
print(merged_summary)

# preparing phenotype dataframe for PLINK
pheno_data<-subset(merged, select=c(FID, IID, Entorhinal_bl, Hippocampus_bl, WholeBrain_bl))
pheno_data_wholebrain<-subset(merged, select=c(FID, IID, WholeBrain_bl))
pheno_data_hippocampus<-subset(merged, select=c(FID, IID, Hippocampus_bl))
pheno_data_entorhinal<-subset(merged, select=c(FID, IID, Entorhinal_bl))
pheno_data_no_na<-pheno_data[!is.na(pheno_data$Hippocampus_bl & pheno_data$WholeBrain_bl),]
write.table(pheno_data, file="phenotype.txt", sep=" ", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(pheno_data_no_na, file="phenotypes_cleaned.txt", sep=" ", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(pheno_data_wholebrain, file="phenotypes_wholebrain.txt", sep=" ", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(pheno_data_hippocampus, file="phenotypes_hippocampus.txt", sep=" ", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(pheno_data_entorhinal, file="phenotypes_entorhinal.txt", sep=" ", quote=FALSE, row.names=FALSE, col.names=TRUE)

adni_qcd<-read.table('adni_qcd.bim', header=FALSE)