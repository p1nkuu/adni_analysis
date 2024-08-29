# plotting PC1 vs PC2
pca_eigenvec<-read.table("pca_qcd_2.eigenvec", header=F)
pca_eigenvec<-subset(pca_eigenvec, select = c(V3, V4))
pheno<-read.table("phenotype.txt", header=T)
covar<-read.table("covar_race.txt", header=T)

merged<-cbind(pheno, pca_eigenvec)
merged2<-cbind(merged, covar$PTRACCAT)
colnames(merged2)[6] <- 'PC1'
colnames(merged2)[7] <- 'PC2'

"
1 = Am Indian/Alaskan
     2 = Asian 
     3 = Hawaiian/Other PI
     4 = Black 
     5 = White 
     6 = More than one 
     7 = Unknown
"

library(ggplot2)
ggplot(merged2, aes(x=PC1, y=PC2, color=factor(covar$PTRACCAT))) +
  geom_point() + labs(title= "PC1 vs PC2 of ADNI", x="PC1", y="PC2", color="Race")

# create new covar file with sex, age, dummy races PC1, PC2, ICV, APOE4
library("haven")
all_subject_data<-data.frame(read_sas("adni_pheno_2023july28.sas7bdat"))
covar<-subset(covar, select = -c(PTRACCAT))
merged2<-subset(merged2, select=c(PC1, PC2, am_indian_alaskan, asian, pacific_islander, black, white, more_than_one, unknown))
ICV_APOE4<-subset(all_subject_data, select=c(APOE4, ICV_bl))
covar2<-cbind(covar, merged2)
new_covar<-cbind(covar2, ICV_APOE4)
write.table(new_covar, file="new_covar.txt", sep=" ", quote=FALSE, row.names=FALSE, col.names=TRUE)

# create plots for PC adjusted gwas
install.packages("qqman")
library(qqman)

wholebrain_pc<-read.table("pc_wholebrain.assoc.linear", header=TRUE)
wholebrain_apoe4<-read.table("pc_apoe4_wholebrain.assoc.linear", header=TRUE)
wholebrain_pc_man<-subset(wholebrain_pc[!is.na(wholebrain_pc$P),], select=c(CHR, SNP, BP, P))
wholebrain_apoe4_man<-subset(wholebrain_apoe4[!is.na(wholebrain_apoe4$P),], select=c(CHR, SNP, BP, P))
manhattan(wholebrain_pc_man, chr="CHR", bp="BP", snp="SNP", p="P", main="Whole Brain (adj. sex, age, PC1, PC2)", annotatePval = 0.00001)
manhattan(wholebrain_apoe4, chr="CHR", bp="BP", snp="SNP", p="P", main="Whole Brain (adj. sex, age, PC1, PC2, APOE4)", annotatePval = 0.00001)

entorhinal_pc<-read.table("pc_entorhinal.assoc.linear", header=TRUE)
entorhinal_apoe4<-read.table("pc_apoe4_entorhinal.assoc.linear", header=TRUE)
entorhinal_pc_man<-subset(entorhinal_pc[!is.na(entorhinal_pc$P),], select=c(CHR, SNP, BP, P))
entorhinal_apoe4_man<-subset(entorhinal_apoe4[!is.na(entorhinal_apoe4$P),], select=c(CHR, SNP, BP, P))
manhattan(entorhinal_pc_man, chr="CHR", bp="BP", snp="SNP", p="P", main="Entorhinal (adj. sex, age, PC1, PC2)", annotatePval = 0.00001)
manhattan(entorhinal_apoe4_man, chr="CHR", bp="BP", snp="SNP", p="P", main="Entorhinal (adj. sex, age, PC1, PC2, APOE4)", annotatePval = 0.00001)

hippo_pc<-read.table("pc_hippo.assoc.linear", header=TRUE)
hippo_apoe4<-read.table("pc_apoe4_hippo.assoc.linear", header=TRUE)
hippo_pc_man<-subset(hippo_pc[!is.na(entorhinal_pc$P),], select=c(CHR, SNP, BP, P))
hippo_apoe4_man<-subset(hippo_apoe4[!is.na(entorhinal_apoe4$P),], select=c(CHR, SNP, BP, P))
manhattan(hippo_pc_man, chr="CHR", bp="BP", snp="SNP", p="P", main="Hippocampus (adj. sex, age, PC1, PC2)", annotatePval = 0.00001)
manhattan(hippo_apoe4_man, chr="CHR", bp="BP", snp="SNP", p="P", main="Hippocampus (adj. sex, age, PC1, PC2, APOE4)", annotatePval = 0.00001)