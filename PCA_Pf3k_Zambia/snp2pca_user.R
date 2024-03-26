####################################################
## R script to generate PCA plot using VCF file.   ##
####################################################

##  Please install the R packages first
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("gdsfmt", force = TRUE)
BiocManager::install("SNPRelate", force = TRUE)

## This is the version info of the R packages used 
#SNPRelate_1.16.1 gdsfmt_1.38

library(gdsfmt)
library(SNPRelate)

## Pf3k + Zambia 
#!1) Please change the name of vcf file to the one you want to use
setwd('/Users/giovannacarpi/Manuscripts/Manuscript_Zambia_MIS2018_Pf_genomics/Pf3k_and_Zambia')
input='Pf_capture1001_final_MAFfilter_002.recode.vcf'


#!2) Please change the name of the prefix of output files if needed
name = "pf3k_Zambia_maf0.2"


#!3) Please change the file name if needed
# File format: tab delimited
# 1st column is sample name
# 2nd column is country name 
poplist <-read.table("meta_info.sample_id.country.tsv", sep="\t")
colnames(poplist) <- c("sample", "country")

#!4) Please change the colors if needed
col.list <- c("blue", "forestgreen", "darkslategray3", "red",
              "orchid1", "tan3","wheat4","purple","black") 


vcf.fn <-input

ccm = paste('ccm_', name, '.gds', sep = "")
ccm

## Read the VCF file and save it as GDS format file
snpgdsVCF2GDS(vcf.fn, ccm,  method="biallelic.only")

## Open the SNP GDS file
genofile <- snpgdsOpen(ccm)

#### PCA analysis
## To calculate the eigenvectors and eigenvalues for principal component analysis.
ccm_pca<-snpgdsPCA(genofile, autosome.only=FALSE)

summary(ccm_pca)

## Get data "sample.id" from a GDS node
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

EV1 = ccm_pca$eigenvect[,1]    # the first eigenvector
EV2 = ccm_pca$eigenvect[,2]    # the second eigenvector

## Get a list of countries following the order of the sample IDs
pop = factor(poplist$country)[match(ccm_pca$sample.id, poplist$sample)]

tab <- data.frame(ccm_pca$sample.id,pop, EV1, EV2, stringsAsFactors = FALSE)
write.table(tab, file=paste(name,'-PCA.tsv', sep=""), 
                  quote = F, row.names = F, sep="\t")
head(tab)

list(ccm_pca$eigenval[1:10])

## Variance explained by 32 PCs
ccm_pca$varprop[1:32]*100



figpdf = paste(name,'-PCA.pdf', sep="")
pdf(file = figpdf)
#plot(EV1, EV2,xlab="PC1", ylab="PC2", col=col.list[as.integer(pop)], pch=19,cex=0.3)
#legend("bottomright", legend=levels(pop), bg="transparent",pch=19, cex=1,col=col.list,text.col=col.list)

plot(EV1, EV2, xlab="PC1", ylab="PC2",xlim= c(-0.06, 0.06),ylim= c(-0.05, 0.05), col=col.list[as.integer(pop)], pch=19,cex=0.3)
legend("bottomleft", legend=levels(pop), bg="transparent",pch=19, cex=0.7,col=col.list,text.col=col.list)
abline(v=0, h= 0, col="grey", cex= 0.3, lty = "dashed")
#text(EV1, EV2, pos = 1, labels = sample.id, cex =0.7 )


dev.off()
sessionInfo()



## East Africa: Pf3k (DRC, Malawi, Tanzania) + Zambia

#!1) Please change the name of vcf file to the one you want to use
setwd('/Users/giovannacarpi/Manuscripts/Manuscript_Zambia_MIS2018_Pf_genomics/Pf3k_and_Zambia/Pf_capture.EastAfrica')
input='Pf_capture.EastAfrica.missness.MAFfilter_0.02.recode.vcf'


#!2) Please change the name of the prefix of output files if needed
name = "EastAfrica_maf0.2"


#!3) Please change the file name if needed
# File format: tab delimited
# 1st column is sample name
# 2nd column is country name 
poplist_EA <-read.table("EastAfrica_meta_info.sample_id.country.tsv", sep="\t")
colnames(poplist_EA) <- c("sample", "country")

#!4) Please change the colors if needed
col.list_EA <- c("blue", "red",
                 "orchid1", "tan3","wheat4","purple","black") 

vcf.fn <-input

ccm_EA = paste('ccm_EA_', name, '.gds', sep = "")
ccm_EA

## Read the VCF file and save it as GDS format file
snpgdsVCF2GDS(vcf.fn, ccm_EA,  method="biallelic.only")

## Open the SNP GDS file
genofile <- snpgdsOpen(ccm_EA)

#### PCA analysis
## To calculate the eigenvectors and eigenvalues for principal component analysis.
ccm_EA_pca<-snpgdsPCA(genofile, autosome.only=FALSE)

summary(ccm_EA_pca)

## Get data "sample.id" from a GDS node
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

EV1 = ccm_EA_pca$eigenvect[,1]    # the first eigenvector
EV2 = ccm_EA_pca$eigenvect[,2]    # the second eigenvector

## Get a list of countries following the order of the sample IDs
popEA = factor(poplist_EA$country)[match(ccm_EA_pca$sample.id, poplist_EA$sample)]

tab <- data.frame(ccm_EA_pca$sample.id, popEA, EV1, EV2, stringsAsFactors = FALSE)
write.table(tab, file=paste(name,'-PCA.tsv', sep=""), 
            quote = F, row.names = F, sep="\t")
head(tab)

list(ccm_EA_pca$eigenval[1:10])
list(ccm_EA_pca$eigenval)
tab

figpdf = paste(name,'-PCA.pdf', sep="")
pdf(file = figpdf)
plot(EV1, EV2, xlab="PC1", ylab="PC2", col=col.list_EA[as.integer(popEA)], pch=19,cex=0.3)
legend("bottomright", legend=levels(popEA), bg="transparent",pch=19, cex=1,col=col.list_EA,text.col=col.list_EA)
abline(v=0, h= 0, col="grey", cex= 0.3, lty = "dashed")


#variance explained by 32 PCs
ccm_EA_pca$varprop[1:32]*100

dev.off()
sessionInfo()