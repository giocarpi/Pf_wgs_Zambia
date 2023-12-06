# IBD analyses----------

library(tidyverse)
library(isoRelate)
#Uses big data format and load the following library() 
library(doMC)
registerDoMC(cores = 8)

## Generate ped and map files----------

## input files are store in the input folder

Pf241_samples<- "input/Pf_capture241_final_MAFfilter_002.recode.vcf" 

### Step 1. convert vcf to gds----
library(SeqArray)
vcf_header <- seqVCF_Header(Pf241_samples)
# info columns to retain
info.import <- c("AC", "AF", "AN", "BaseQRankSum", "DP", "DS",
                 "ExcessHet", "FS", "InbreedingCoeff", "MQ",
                 "MQRankSum", "QD", "ReadPosRankSum", "SOR", "EFF")
# format columns to retain
format.import <- c("AD", "DP", "GQ", "GT", "PL", "RGQ", "SB")

# convert VCF to GDS

Pf241_samples_gds<-seqVCF2GDS(Pf241_samples,
                              "Pf241_samples.gds",genotype.var.name="GT",
                              header=vcf_header, info.import=info.import,
                              fmt.import=format.import)
# Use isolates file for downstream analysis 
isolates241samples_coregenome <- seqOpen(Pf241_samples_gds) 

## fws_metadata.csv is the meta data that contains population, fws and moi info
fws_metadata <- read.csv("input/fws_metadata.csv")

### Step 2. extract ped from gds----
library(moimix)
# create map and ped file for isoRelate 
# moi.estimates once set through fws caluclation
# use.hets is redundant
# for moi=1 isolates, het markers will be coded as missing
# moi=2 isolates, het markers are allowed
# sample.id <- seqGetData(isolates241samples_coregenome, "sample.id")
zambia_ped_map_file_241samples<-extractPED(isolates241samples_coregenome, 
                                           moi.estimates = fws_metadata$moi, 
                                           outfile = NULL) 
# Change the recombination rate according 2016 paper
zambia_ped_map_file_241samples$map$genetic_distance<-zambia_ped_map_file_241samples$map$pos/13500
# Change famID to population info
zambia_ped_map_file_241samples$ped$famID <- fws_metadata$fid

saveRDS(zambia_ped_map_file_241samples, file = "input/zambia_ped_map_file_241samples.rds")

zambia_ped_map_file_241samples <-readRDS("input/zambia_ped_map_file_241samples.rds")
# get separate EFF info and save as a separate meta data file

EFFinfo <- seqGetData(isolates241samples_coregenome, 
                      c("variant.id", "position", "chromosome","annotation/info/EFF"))
EFFinfo<-as.data.frame(EFFinfo)
colnames(EFFinfo)<-c("variant.id","pos_bp","chr","EFF")
EFFinfo<-EFFinfo %>% tidyr::separate(EFF,into=c("Effect","Effect_Impact",
   "Functional_Class","Codon_Change","Amino_Acid_Change",
   "Gene_Name","Transcript_BioType","Gene_Coding","Transcript_ID","Exon_Rank ",
   "Genotype_Number"), sep="[(|)]",remove=F,convert=F)

write.csv(EFFinfo, "EffInfo.csv", row.names=F)
EFFinfo<- read.csv("EffInfo.csv")
# close the gds connection
seqClose(isolates241samples_coregenome) 

## Genome-level significance--------------
# remove MOI=1 samples, to ensure diploid scenario

xmoi2<- samples241_final_genotypes$genotypes[,6:ncol(samples241_final_genotypes$genotypes)]
xmoi2 <- t(as.matrix(xmoi2))
xmoi2 <- xmoi2[fws_metadata$moi==2,]
# calculate pairwise correlation for all the SNP sites
cld_moi2 <- cor(xmoi2,use="pairwise.complete.obs")

# get eigen values for every 1000 sites, 
# and calculate how many independent SNPs can be recovered
eigen995 <- c()
for (i in seq(1,27230,100)){
  eigen_result <- eigen(cld_moi2[i:(i+999),i:(i+999)])
  tempEigen <- cumsum(eigen_result$values)/1000
  eigen995 <- c(eigen995,min(which(tempEigen>0.995)))
}


## isoRelate-------

##Next we estimate the model parameters. These parameters can be useful as an initial measure of relatedness.
# Estimates the number of meiosis and the probabilities of sharing 0, 1 and 2 alleles IBD between all pairwise combinations of isolates

samples241_final_genotypes<-
  getGenotypes(ped.map = zambia_ped_map_file_241samples,
  reference.ped.map = NULL,
  maf = 0.01, # you can modify the threshold based on you data type
  isolate.max.missing = 0.2,# you can modify the threshold based on you data type
  snp.max.missing = 0.1,# you can modify the threshold based on you data type
  chromosomes = NULL,
  input.map.distance = "cM",
  reference.map.distance = "cM")#  centimorgan 

samples241_final_parameters <- getIBDparameters(ped.genotypes = samples241_final_genotypes,
                                                number.cores = 8) 

# ibd=1 allele and 2 alleles parameter distribution
moi1_samples<-as.vector(moi1_samples$V1)
ggplot(samples241_final_parameters%>%filter(iid1 %in% moi1_samples,
                                            iid2 %in% moi1_samples))+
  aes(ibd1)+
  geom_histogram()

ggplot(samples241_final_parameters)+
  aes(ibd1)+
  geom_histogram()

saveRDS(samples241_final_genotypes, file = "intermediate/samples241_final_genotypes.rds")

samples241_final_genotypes<-readRDS("intermediate/samples241_final_genotypes.rds")

write.csv(samples241_final_parameters, "intermediate/samples241_final_parameters.csv",row.names=F) 
samples241_final_parameters<-read.csv("intermediate/samples241_final_parameters.csv")
### Detect IBD Segments-------

#Following parameter estimation we can detect IBD segments
#detects genomic regions shared IBD between all pairwise combinations of isolates.
#We have set thresholds on the minimum number of SNPs and length of IBD segments reported in order to reduce false positive IBD calls that are most likely due to population linkage disequilibrium from extremely distant relatedness.

samples241_final_ibd_segments<-
  getIBDsegments(ped.genotypes= samples241_final_genotypes,                    
   parameters = samples241_final_parameters, 
   number.cores = 8, #The number of cores used for parallel execution
   minimum.snps = 20, #the minimum number of SNPs in an IBD segment for it to be reported
   minimum.length.bp = 50000,# The minimum length of a reported IBD segment
   error = 0.001) #The genotyping error rate

# 8249 pairs inferred IBD
# 9501 IBD segments detected
### IBD segment lengths distribution -----------
p1<-ggplot(samples241_final_ibd_segments,aes(length_M*100))+
  geom_histogram(bins=30)+
  geom_vline(xintercept = median(samples241_final_ibd_segments$length_M)*100,linetype=2,col="orange")+
  #geom_density()+
  xlab("Genetic distance of segment (cM)")+
  theme_bw(base_size = 16)

p2<-ggplot(samples241_final_ibd_segments,aes(length_bp/1000))+
  geom_histogram(bins=30)+
  geom_vline(xintercept = median(samples241_final_ibd_segments$length_bp)/1000,linetype=2,col="orange")+
  #geom_density()+
  xlab("Length of segment (Kb)")+
  theme_bw(base_size = 16)

figS13 <- p1 + inset_element(p2, 0.5, 0.5, 0.95, 0.95)
ggsave("output/FigS13_IBD_segments_dist.pdf",figS13,width=8,height=6)
ggsave("output/FigS13_IBD_segments_dist.png",figS13,width=8,height=6)

write.csv(samples241_final_ibd_segments, file = "intermediate/samples241_final_ibd_segments.csv",row.names=F)

samples241_final_ibd_segments<-read.csv("intermediate/samples241_final_ibd_segments.csv")
  
  
getIBDsummary(ped.genotypes = samples241_final_genotypes, 
              ibd.segments = samples241_final_ibd_segments) 

# 10K, 20snp
# 241 / 241 isolates IBD
# 8249 / 28920 pairs of isolates IBD
# total IBD segments detected = 9501 
# average length of segments (bp) = 61364 
# 81 IBD segments detected on chromosome Pf3D7_01_v3 
# 144 IBD segments detected on chromosome Pf3D7_02_v3 
# 203 IBD segments detected on chromosome Pf3D7_03_v3 
# 350 IBD segments detected on chromosome Pf3D7_04_v3 
# 140 IBD segments detected on chromosome Pf3D7_05_v3 
# 271 IBD segments detected on chromosome Pf3D7_06_v3 
# 233 IBD segments detected on chromosome Pf3D7_07_v3 
# 210 IBD segments detected on chromosome Pf3D7_08_v3 
# 38 IBD segments detected on chromosome Pf3D7_09_v3 
# 514 IBD segments detected on chromosome Pf3D7_10_v3 
# 76 IBD segments detected on chromosome Pf3D7_11_v3 
# 178 IBD segments detected on chromosome Pf3D7_12_v3 
# 180 IBD segments detected on chromosome Pf3D7_13_v3 
# 67 IBD segments detected on chromosome Pf3D7_14_v3 
# 6816 IBD segments detected on chromosome Pf3D7_API_v3 

# 50K, 20snp
# 231 / 241 isolates IBD
# 1145 / 28920 pairs of isolates IBD
# total IBD segments detected = 1443 
# average length of segments (bp) = 273755 
# 36 IBD segments detected on chromosome Pf3D7_01_v3 
# 76 IBD segments detected on chromosome Pf3D7_02_v3 
# 178 IBD segments detected on chromosome Pf3D7_03_v3 
# 90 IBD segments detected on chromosome Pf3D7_04_v3 
# 74 IBD segments detected on chromosome Pf3D7_05_v3 
# 240 IBD segments detected on chromosome Pf3D7_06_v3 
# 60 IBD segments detected on chromosome Pf3D7_07_v3 
# 121 IBD segments detected on chromosome Pf3D7_08_v3 
# 28 IBD segments detected on chromosome Pf3D7_09_v3 
# 212 IBD segments detected on chromosome Pf3D7_10_v3 
# 46 IBD segments detected on chromosome Pf3D7_11_v3 
# 166 IBD segments detected on chromosome Pf3D7_12_v3 
# 57 IBD segments detected on chromosome Pf3D7_13_v3 
# 59 IBD segments detected on chromosome Pf3D7_14_v3 

samples241_final_ibd_segments2<-samples241_final_ibd_segments
samples241_final_ibd_segments2<-
  samples241_final_ibd_segments2[,c(3,4,1,2,5:ncol(samples241_final_ibd_segments2))]
colnames(samples241_final_ibd_segments2)[1:4]<-c("fid1","iid1","fid2","iid2")
segmentStats<- rbind(samples241_final_ibd_segments,samples241_final_ibd_segments2)
segmentStats<- segmentStats%>%distinct(iid1,iid2)%>%count(iid1)
  ggplot(segmentStats,aes(n))+geom_histogram()


### IBD patterns----------

# remove Chromosome API from the following calculations
samples241_final_genotypes$genotypes <- samples241_final_genotypes$genotypes %>% 
  filter(chr != "Pf3D7_API_v3")
samples241_final_ibd_segments <- samples241_final_ibd_segments %>%
  filter(chr != "Pf3D7_API_v3")


samples241_final_IBD_matrix<- getIBDmatrix(ped.genotypes = samples241_final_genotypes, 
                                           ibd.segments = samples241_final_ibd_segments)
#write.csv(samples241_final_IBD_matrix,"intermediate/samples241_final_IBD_matrix.csv",row.names=F)


# calculate the proportion of pairs IBD at each SNP for all samples together 

samples241_final_proportion <- getIBDproportion(ped.genotypes = samples241_final_genotypes,
                                                ibd.matrix = samples241_final_IBD_matrix, 
                                                groups = fws_metadata[,1:3])

plotIBDproportions(samples241_final_proportion%>%
                     filter(subpop=="Western.Western"),
                   interval=c("Pf3D7_12_v3",0,1980224),
                   subpop.facet =T)

subpop<-unique(samples241_groups_proportion_final$subpop)
fws_metadata$pid<-fws_metadata$fid
fws_metadata$pid<-gsub("-",".",fws_metadata$pid)
popList<- unique(fws_metadata$pid)
popList<-paste(popList,popList,sep=".")
samples241_groups_proportion_final_sub <- subset(samples241_final_proportion,subpop %in% popList)

plotIBDproportions(samples241_groups_proportion_final_sub,
                   interval=c("Pf3D7_12_v3",0,1980224),
                   subpop.facet =T)

samples241_final_proportion$subpop <- gsub("North.Western","NorthWestern",samples241_final_proportion$subpop, fixed=T)

samples241_final_proportion <- samples241_final_proportion %>% 
  separate(subpop, into=c("pop1","pop2"), sep="[.]",remove=F)


chrList <- unique(samples241_final_genotypes$genotypes$chr)
for (chr1 in chrList){
  print(chr1)
  p <- ggplot(samples241_final_proportion %>% filter(chr==chr1,pop1 != "Other",pop2!="Other"))+
    aes(pos_bp,prop_ibd)+
    geom_point()+
    ggtitle(chr1)+
    facet_grid(pop1~pop2)

  ggsave(paste(chr1,"_pwpop.pdf",sep=""),p)
}


### calculate the significance of IBD sharing------------       

# To identify genomic loci with significant amounts of excess IBD we can apply a transformation to the binary IBD matrix to 
# account for variations in isolate relatedness as well as SNP allele frequencies, then calculate a summary statistic 
#  that can be used to assess significance.

#### all pops together--------

samples241_final_iR <- getIBDiR(ped.genotypes =samples241_final_genotypes,
                                      ibd.matrix =  samples241_final_IBD_matrix, 
                                      groups = NULL)

write.csv(samples241_final_iR,"output/samples241_final_iR.csv",row.names=F)

samples241_final_iR <- read.csv("output/samples241_final_iR.csv")

plotIBDiR(samples241_final_iR,highlight.genes = highlight_genes) 
library(qqman)
data(highlight_genes)
samples241_final_iR$CHR <- str_split_fixed(samples241_final_iR$chr,"_",3)[,2]
samples241_final_iR$P<-10^(-samples241_final_iR$log10_pvalue)
#samples241_final_iR <- samples241_final_iR%>%filter(CHR != "API")
samples241_final_iR$CHR <-as.numeric(samples241_final_iR$CHR)

samples241_final_iR <- samples241_final_iR %>%
  left_join(EFFinfo,by=c("chr","pos_bp"))

#samples241_final_iR <- samples241_final_iR[,1:(ncol(samples241_final_iR)-1)]
#(column= chromosome number 1, 2, 3, etc,column 2= bp number, column3 = snp_id like and Pf3D7_11_v3:592805 and last column= P.value. NB you have to change log10(p.value) result from isoRelate to p.value =1/10^log10(p.value))

write.csv(samples241_final_iR,"output/samples241_final_iR.csv", row.names=F )

write.csv(samples241_final_iR%>%filter(log10_pvalue>5),"output/ibd_sig_sites.csv", row.names=F )


# manhattan function do not annotate all the genes, it only annotate the top one per chromosome
samples241_final_iR%>% filter(CHR==6,log10_pvalue>5)%>%
  select(pos_bp,Effect_Impact,Functional_Class,Amino_Acid_Change,Gene_Name)

manhattan(samples241_final_iR, col = c("black", "royalblue"), annotateTop= T, annotatePval=-log10(1e-5),suggestiveline = F, genomewideline = -log10(1e-5), ylim= c(0,20), chr="CHR", bp="pos_bp", snp="Gene_Name", p="P" , main= "Significant IBD sharing")
debug(manhattan)
undebug(manhattan)
#### separate regions--------

popList<- unique(fws_metadata$pid)
regions <- data.frame(pid = popList,region = c("E","E","W","W","W","E","E","W"))
fws_metadata <- fws_metadata %>% left_join(regions)

samples241_final_iR_group <- getIBDiR(ped.genotypes =samples241_final_genotypes,
                                ibd.matrix =  samples241_final_IBD_matrix, 
                                groups = fws_metadata[,c("fid","iid","region")])
library(cowplot)

samples241_final_iR_group <- samples241_final_iR_group %>% separate(subpop,into=c("pop1","pop2"),sep="[/]",remove=F)
samples241_final_iR_group$CHR <- str_split_fixed(samples241_final_iR_group$chr,"_",3)[,2]
samples241_final_iR_group$P<-10^(-samples241_final_iR_group$log10_pvalue)
#samples241_final_iR <- samples241_final_iR%>%filter(CHR != "API")
samples241_final_iR_group$CHR <-as.numeric(samples241_final_iR_group$CHR)

#(column= chromosome number 1, 2, 3, etc,column 2= bp number, column3 = snp_id like and Pf3D7_11_v3:592805 and last column= P.value. NB you have to change log10(p.value) result from isoRelate to p.value =1/10^log10(p.value))
samples241_final_iR_group<-samples241_final_iR_group%>%left_join(EFFinfo,by=c("chr","pos_bp"))



#samples241_final_iR_group <- samples241_final_iR_group %>% arrange(chr,pos_bp,pop1,pop2)
ggplot(samples241_final_iR_group %>% filter(pop1 != "Other",pop2!="Other"))+
  aes(pos_bp,log10_pvalue)+
  geom_point()+
  facet_grid(subpop~chr,scales="free_x")+
  theme_cowplot()
ggsave("iR.pdf",p1,width=15,height=20)


samples241_final_iR$POSITION<- samples241_final_iR$pos_bp
samples241_final_iR_group$POSITION<- samples241_final_iR_group$pos_bp
figpdf <- 'output/IBD_sig.pdf'
pdf(file = figpdf)
par(mfrow = c(3, 1), mar = c(3, 6, 1.5, 4))
manhattanplot(samples241_final_iR,
              pval = TRUE,
              inset = 10,
              threshold = 5,
              mrk= NULL,
              mrk.cex=2,
              mrk.lab.cex=2,
              cr = ibdRegion,
              ylab = expression(-log[10](italic(p))),
              xlab = "",
              main = expression(paste("Significance of ",X[iR], " in All Provinces",sep="")),
              mgp=c(2,0.5,0)
)
manhattanplot(samples241_final_iR_group%>% filter(subpop=="E/E"),
              pval = TRUE,
              inset = 10,
              threshold = 5,
              mrk= NULL,
              mrk.cex=2,
              mrk.lab.cex=2,
              cr = NULL,
              ylab = expression(-log[10](italic(p))),
              xlab = "",
              main = expression(paste("Significance of ",X[iR], " in Eastern Provinces",sep="")),
              mgp=c(2,0.5,0)
)
manhattanplot(samples241_final_iR_group%>% filter(subpop=="W/W"),
              pval = TRUE,
              inset = 10,
              threshold = 5,
              mrk= NULL,
              mrk.cex=2,
              mrk.lab.cex=2,
              cr = NULL,
              ylab = expression(-log[10](italic(p))),
              xlab = "",
              main = expression(paste("Significance of ",X[iR], " in Western populations",sep="")),
              mgp=c(2,0.5,0)
)
dev.off()

figpdf <- 'output/IBD_sig.pdf'
pdf(file = figpdf)
par(mfrow = c(3, 1), mar = c(5.1, 4.1, 4.1, 2.1))
manhattan(samples241_final_iR, col = c("black", "royalblue"), annotateTop= T, annotatePval=NULL,suggestiveline = F, genomewideline = -log10(1e-5), ylim= c(0,20), chr="CHR", bp="pos_bp", snp="Gene_Name", p="P" , main= "All Provinces")
manhattan(samples241_final_iR_group %>% filter(subpop=="E/E"), col = c("black", "royalblue"), annotateTop= T, annotatePval=NULL,suggestiveline = F, genomewideline = -log10(1e-5), ylim= c(0,20), chr="CHR", bp="pos_bp", snp="Gene_Name", p="P" , main= "Eastern Provinces")
          manhattan(samples241_final_iR_group %>% filter(subpop=="W/W"), col = c("black", "royalblue"), annotateTop= T, annotatePval=NULL,suggestiveline = F, genomewideline = -log10(1e-5), ylim= c(0,20), chr="CHR", bp="pos_bp", snp="Gene_Name", p="P" , main= "Western Provinces")
dev.off()                    

#### MOI=1 vs MOI>1-----------

samples241_final_iR_moi <- getIBDiR(ped.genotypes =samples241_final_genotypes,
                                    ibd.matrix =  samples241_final_IBD_matrix, 
                                    groups = fws_metadata[,c("fid","iid","moi")])


samples241_final_iR_moi <- samples241_final_iR_moi %>% separate(subpop,into=c("pop1","pop2"),sep="[/]",remove=F)
samples241_final_iR_moi$CHR <- str_split_fixed(samples241_final_iR_moi$chr,"_",3)[,2]
samples241_final_iR_moi$P<-10^(-samples241_final_iR_moi$log10_pvalue)
#samples241_final_iR <- samples241_final_iR%>%filter(CHR != "API")
samples241_final_iR_moi$CHR <-as.numeric(samples241_final_iR_moi$CHR)
samples241_final_iR_moi$POSITION<- samples241_final_iR_moi$pos_bp

figpdf <- 'output/IBD_moi1v2_sig.pdf'
pdf(file = figpdf)
#figpdf <- 'output/IBD_moi1v2_sig.png'
#png(file = figpdf,units ="in",width=6,height=15,res=300)
par(mfrow = c(3, 1), mar = c(3, 6, 1.5, 4))
manhattanplot(samples241_final_iR,
              pval = TRUE,
              inset = 10,
              threshold = 5,
              mrk= NULL,
              mrk.cex=2,
              mrk.lab.cex=2,
              cr = ibdRegion,
              ylab = expression(-log[10](italic(p))),
              xlab = "",
              main = expression(paste("Significance of ",X[iR], " in All Samples",sep="")),
              mgp=c(2,0.5,0)
)
manhattanplot(samples241_final_iR_moi%>% filter(subpop=="1/1"),
              pval = TRUE,
              inset = 10,
              threshold = 5,
              mrk= NULL,
              mrk.cex=2,
              mrk.lab.cex=2,
              cr = NULL,
              ylab = expression(-log[10](italic(p))),
              xlab = "",
              main = expression(paste("Significance of ",X[iR], " in Monoclonal Samples",sep="")),
              mgp=c(2,0.5,0)
)
manhattanplot(samples241_final_iR_moi%>% filter(subpop=="2/2"),
              pval = TRUE,
              inset = 10,
              threshold = 5,
              mrk= NULL,
              mrk.cex=2,
              mrk.lab.cex=2,
              cr = NULL,
              ylab = expression(-log[10](italic(p))),
              xlab = "",
              main = expression(paste("Significance of ",X[iR], " in Polyclonal Samples",sep="")),
              mgp=c(2,0.5,0)
)
dev.off()

#### Masking high snp density regions----------------

samples241_downSample_genotypes<-samples241_final_genotypes
samples241_downSample_genotypes$genotypes<-
  samples241_downSample_genotypes$genotypes%>%filter(((chr=="Pf3D7_10_v3")&(pos_bp>1432498)&(pos_bp<1434786))==F,
                                                     ((chr=="Pf3D7_03_v3")&(pos_bp>125687)&(pos_bp<130112))==F
  )


samples241_downSample_parameters <- getIBDparameters(ped.genotypes = samples241_downSample_genotypes,
                                                     number.cores = 8) 

# Infer ibd segments
samples241_downSample_ibd_segments<-
  getIBDsegments(ped.genotypes= samples241_downSample_genotypes,                    
                 parameters = samples241_downSample_parameters, 
                 number.cores = 8, #The number of cores used for parallel execution
                 minimum.snps = 20, #the minimum number of SNPs in an IBD segment for it to be reported
                 minimum.length.bp = 50000,# The minimum length of a reported IBD segment
                 error = 0.001) #The genotyping error rate

getIBDsummary(ped.genotypes = samples241_downSample_genotypes, 
              ibd.segments = samples241_downSample_ibd_segments) 

# calculate IBD matrix
samples241_downSample_IBD_matrix<- getIBDmatrix(ped.genotypes = samples241_downSample_genotypes, 
                                                ibd.segments = samples241_downSample_ibd_segments)
# calculate IBD XiR
samples241_downSample_iR <- getIBDiR(ped.genotypes =samples241_downSample_genotypes,
                                     ibd.matrix =  samples241_downSample_IBD_matrix, 
                                     groups = NULL)

samples241_downSample_iR$CHR <- str_split_fixed(samples241_downSample_iR$chr,"_",3)[,2]
samples241_downSample_iR$P<-10^(-samples241_downSample_iR$log10_pvalue)
#samples241_downSample_iR <- samples241_downSample_iR%>%filter(CHR != "API")
samples241_downSample_iR$CHR <-as.numeric(samples241_downSample_iR$CHR)

samples241_downSample_iR$POSITION<- samples241_downSample_iR$pos_bp
manhattanplot(samples241_downSample_iR,
              pval = TRUE,
              inset = 10,
              threshold = 5,
              mrk= NULL,
              mrk.cex=2,
              mrk.lab.cex=2,
              cr = NULL,
              ylab = expression(-log[10](italic(p))),
              xlab = "",
              main = expression(paste("Significance of ",X[iR], " in All Provinces",sep="")),
              mgp=c(2,0.5,0)
)



## pw isolates network---------
### get species meta data---------

allMeta <- readxl::read_excel("input/Complete_459samples_metadata_revised_byGC_09292022_FINAL.xlsx")
fws_metadata <- fws_metadata %>% left_join(allMeta%>%dplyr::select(Sample_id,province,district,ward,cluster_id,cluster_latitude,cluster_longitude), by=c("iid"="Sample_id"))

head(fws_metadata)
write.csv(fws_metadata,"fws_metadata_geo.csv",row.names=F)
fws_metadata<-read.csv("fws_metadata_geo.csv")
### plot isolate networks------------

## get proportional IBD

#head(samples241_final_genotypes$genotypes[1:10,1:6])
totalLength<- samples241_final_genotypes$genotypes %>% 
  group_by(chr) %>%
  summarize(minbp = min(pos_bp),maxbp=max(pos_bp),bprange = maxbp-minbp+1)

totalSize <- sum(totalLength$bprange)

ibd_stats <- samples241_final_ibd_segments %>% group_by(iid1,iid2) %>%
  summarize(sharedLength = sum(end_position_bp-start_position_bp+1), 
            weight = sharedLength/totalSize) %>%
  ungroup()%>%
  arrange(desc(weight))

write.csv(ibd_stats, "intermediate/ibd_stats_pw_isolates.csv")

ibd_stats_by_chr <- samples241_final_ibd_segments %>% 
  group_by(iid1,iid2,chr) %>%
  summarize(sharedLength = sum(end_position_bp-start_position_bp+1)) %>%
  left_join(totalLength %>% dplyr::select(chr, bprange),by="chr") %>%
  mutate(weight = sharedLength/bprange) %>%
  ungroup()%>%group_by(chr)%>%
  arrange(desc(weight))

write.csv(ibd_stats_by_chr, "intermediate/ibd_stats_pw_isolates_by_chr.csv")


library(igraph)
color_net <-c("red","orange", "tomato4","mediumvioletred","darkgreen", "darkblue", "grey","dodgerblue")

constructIsolateGraph <-function(graphDat, vertices, pop, moi,
                                 weightCutoff = 0.05) {
  if(! "weight" %in% colnames(graphDat)){
    stop("weight column is not defined in the data frame")
  }
  ibdGraph <- graph_from_data_frame(graphDat, directed=F,
                                    vertices = vertices)
  pop<-as.factor(pop)
  V(ibdGraph)$pop <- pop
  V(ibdGraph)$moi <- moi
  
  ibd_05 <- delete.edges(ibdGraph, which(E(ibdGraph)$weight < weightCutoff))
  return(ibd_05)
  
}

plotIsolates <- function(g,col=NULL,
                         title="",plotLegend=T) {
  
  if (is.null(col)){
    col <- brewer.pal(n = length(levels(V(g)$pop)),name = "YlOrRd")
  }
  
  plot(g,vertex.size=(2+(V(g)$moi))*1.2,vertex.color = col[V(g)$pop],
       vertex.label=NA,edge.width = E(g)$weight*5,edge.arrow.size=0, edge.color = "grey20",
       edge.curved = 0.1, layout=layout.fruchterman.reingold,main = title,vertex.frame.color=NA)
 
  if(plotLegend){
    legend(
      x=-1.5,y=-0.2,
      legend = levels(V(g)$pop),
      col  = col,
      pch    = 19,
      cex    = 1,
      bty    = "n",
      title  = "Province",
      title.adj =0.1
    )
    
    legend(
      x=-1.5,y=1,
      legend = c(0.05,0.1,0.25,0.5,1),
      lwd=c(0.05,0.1,0.25,0.5,1)*5, 
      cex=1, 
      col="grey20",
      bty    = "n",
      title  = "IBD Proportion"
    )
    
    legend(
      x=-1.5,y=0.3,
      legend = c("=1",">1"),
      pt.bg  = "grey20",
      pch    = 21,
      pt.cex =c(1,1.5), 
      col="grey20",
      bty    = "n",
      title  = "MOI"
    )
    
  } 
  
}


ibd_stats_by_chr <- ibd_stats_by_chr %>% ungroup()


#### plot network per chr---------
chrList <- unique(samples241_final_genotypes$genotypes$chr)
for (chr1 in chrList){

  tempNet <- constructIsolateGraph(ibd_stats_by_chr%>%filter(chr==chr1)%>%dplyr::select(iid1,iid2,weight), fws_metadata$iid, pop=fws_metadata$pid,moi=samples241_final_genotypes$pedigree$moi)
  
plotIsolates(tempNet,
             col=color_net,title=chr1)

}
ibd_stats<-ibd_stats%>%left_join(fws_metadata%>%select(iid,moi),by=c("iid1"="iid"))
ibd_stats<-ibd_stats%>%left_join(fws_metadata%>%select(iid,moi),by=c("iid2"="iid"))
ibd_stats%>%filter(weight>=0.5)
#### plot all pw network ---------

ibd_02_all <- constructIsolateGraph(ibd_stats%>%dplyr::select(iid1,iid2,weight),
                                    fws_metadata$iid, pop=fws_metadata$pid,
                                    moi=samples241_final_genotypes$pedigree$moi,
                                    weightCutoff = 0.02)


ibd_05_all <- constructIsolateGraph(ibd_stats%>%dplyr::select(iid1,iid2,weight),
             fws_metadata$iid, pop=fws_metadata$pid,
             moi=samples241_final_genotypes$pedigree$moi)


figpdf <- 'output/isolates_network.pdf'
pdf(file = figpdf,width=9.6,height=7)
plotIsolates(ibd_05_all,col=color_net)
dev.off()
plotIsolates(ibd_02_all,col=color_net)


#### plot sig region network ---------
pdf("output/SigRegion_networks.pdf",width=12,height=8)
png("output/SigRegion_networks.png",width=12,height=8,units = "in",res=300)
par(mfrow = c(2, 3),mar=c(0,0,2,0),oma=c(1,1,1,6),xpd=T)
for (i in 1:nrow(ibdRegion)){
  #browser()
  tempChr <- chrList[ibdRegion$CHR[i]]
  tempStart <- ibdRegion$START[i]
  tempEnd <- ibdRegion$END[i]
  
  tempRegionStats <- samples241_final_ibd_segments %>% filter(chr==tempChr)%>%
    group_by(iid1,iid2) %>% rowwise() %>%
    mutate(overlapLen = min(end_position_bp,tempEnd)-max(start_position_bp,tempStart)+1,
           adjustLen = ifelse(overlapLen>0,overlapLen,0))%>%
    summarize(sharedLength = sum(adjustLen))%>%
    mutate(weight = sharedLength/(tempEnd-tempStart+1)) 
  
  tempNet <- constructIsolateGraph(tempRegionStats%>%dplyr::select(iid1,iid2,weight), fws_metadata$iid, pop=fws_metadata$region,moi=samples241_final_genotypes$pedigree$moi, weightCutoff=0.05)
  
  tempName1 <- paste("Chr_",ibdRegion$CHR[i],":",tempStart/1000,"-",tempEnd/1000," Kb",sep="")
  #figpdf <- paste("output/Chr",ibdRegion$CHR[i],tempStart,tempEnd,sep="_")
  #pdf(file = paste(figpdf,"pdf",sep="."),width=9.6,height=7)
  if(i==1){
    plotIsolates(tempNet,col=c("blue","red"),
                 plotLegend=F,title=tempName1)
    
  }else{
    plotIsolates(tempNet,col=c("blue","red"),
                 plotLegend=F,title=tempName1)
    
  }
 
}
legend(
  par('usr')[2], par('usr')[4]-0.5,  xpd=NA,
  legend = levels(V(tempNet)$pop),
  col  = c("blue","red"),
  pch    = 19,
  cex    = 1,
  bty    = "n",
  title  = "Region",
  title.adj =0.1
)

legend(
  par('usr')[2], par('usr')[4]+0.5,  xpd=NA,
  legend = c(0.25,0.5,1),
  lwd=c(0.25,0.5,1)*5, 
  cex=1, 
  col="grey20",
  bty    = "n",
  title  = "IBD Proportion"
)

legend(
  par('usr')[2], par('usr')[4],  xpd=NA,
  legend = c("=1",">1"),
  pt.bg  = "grey20",
  pch    = 21,
  pt.cex =c(1,1.5), 
  col="grey20",
  bty    = "n",
  title  = "MOI"
)

dev.off()



### plot network on geography-----------
library(sf)

V(ibd_05_all)$long <- fws_metadata$cluster_longitude
V(ibd_05_all)$lat <- fws_metadata$cluster_latitude
dv <- degree(ibd_05_all)
ibd_05_geo  <- delete.vertices(ibd_05_all, names(dv[dv==0]))


st_layers("input/gadm41_ZMB.gpkg")
#border<- read_sf("gadm41_ZMB.gpkg",layer = "ADM_ADM_0")
province<-read_sf("input/gadm41_ZMB.gpkg",layer = "ADM_ADM_1")
#plot(province)

## read in dataframe
ibd05_dat <- as_long_data_frame(ibd_05_geo)
#ibd05_dat <- ibd05_dat%>%mutate(flonm = from_long+runif(nrow(ibd05_dat),-0.2,0.2),
#                                tlonm = to_long+runif(nrow(ibd05_dat),-0.2,0.2),
#                                flatm = from_lat+runif(nrow(ibd05_dat),-0.2,0.2),
#                                tlatm = to_lat+runif(nrow(ibd05_dat),-0.2,0.2))
#ibd05_dat_st <- apply(ibd05_dat, 1, function(x) st_linestring(rbind(as.numeric(x[c("flonm","flatm")]),as.numeric(x[c("tlonm","tlatm")]))),simplify=F)
#fromPoints<- apply(ibd05_dat, 1, function(x) st_point(as.numeric(x[c("from_long","from_lat")])),simplify=F)
#fromPoints<-st_sfc(fromPoints,crs=4326)
#toPoints<-apply(ibd05_dat, 1, function(x) st_point(as.numeric(x[c("to_long","to_lat")])),simplify=F)
#toPoints<-st_sfc(toPoints,crs=4326)
#pointDist <- st_distance(fromPoints,toPoints, by_element = T)
#ibd05_dat$geoDist <- pointDist
#ggplot(ibd05_dat, aes(geoDist,weight))+geom_point()

ibd05_dat_st <- apply(ibd05_dat, 1, function(x) st_linestring(rbind(as.numeric(x[c("from_long","from_lat")]),as.numeric(x[c("to_long","to_lat")]))),simplify=F)


ibd05_dat_st <- st_sfc(ibd05_dat_st,crs=4326)
ibd05_dat_st <- cbind(ibd05_dat_st, ibd05_dat)
ibd05_dat_st <- st_sf(ibd05_dat_st)

ibd05_stats <- ibd05_dat_st %>% count(from_long,from_lat,to_long,to_lat) 
ibd05_points<-data.frame(from_long=V(ibd_05_geo)$long, from_lat=V(ibd_05_geo)$lat) %>% distinct() %>%
  left_join(ibd05_stats%>%dplyr::select(from_long,from_lat,n)) 

ibd05_points[is.na(ibd05_points$n),"n"]<-1

#test1 <- st_linestring(rbind(as.numeric(ibd05_dat[1,c("from_long","from_lat")]),as.numeric(ibd05_dat[1,c("to_long","to_lat")])))
#st_read("")
library(cowplot)
color_map <-c("grey","red","orange", "tomato4","grey","mediumvioletred","darkgreen", "darkblue", "grey","dodgerblue")

p2<-ggplot() +
  geom_sf(data = province, aes(fill=NAME_1),alpha=0.6) +
  #geom_sf(data = ibd05_points_st,aes(),fill="grey")+
  geom_point(data = ibd05_points, aes(x=from_long,y=from_lat,size=n))+
  geom_sf_text(data = province,aes(label=NAME_1))+
  geom_sf(data = ibd05_stats%>%filter(st_is(geometry,"LINESTRING")),aes(),
          linewidth=1,linetype=1,color="grey30")+
  #scale_color_viridis_c(inferno)+
  scale_size(name="No.",breaks=c(1,2,3,8))+
  xlab("")+ylab("")+
  scale_fill_manual(values=color_map,name = "Province",guide="none")+
  ggtitle("No. of Isolate Pairs with IBD > 5%")+
  theme_cowplot()+theme(legend.position = c(.8,.15))

ggsave("output/IBD_map.pdf",p2,width=9.6,height=7)

library(patchwork)
#old_par <- par(mar = c(0, 0, 0, 0), mgp = c(1, 0.25, 0), 
#               bg = NA, cex.axis = 0.75, las = 1, tcl = -0.25)
#par(old_par)
#wrap_elements(panel = ~plotIsolates(ibd_05_all,col=color_net), clip = FALSE)+ 
  theme(plot.margin = margin(5.5, 5.5, 5.5, 35))+p2

## SNP density and significance---------

  ndat <- samples241_final_iR[,c("chr","pos_bp","CHR")]
  ndat$sig <- samples241_final_iR$log10_pvalue>5
  ndat$iR <- samples241_final_iR[,"iR"]
  ndat$P <- samples241_final_iR[,"P"]
  colnames(ndat)[2] <- "BP"
  ndat$count<-NA
  for (i in 1:nrow(ndat)){
    dchr<-ndat$chr[i]
    rangel<-ndat$BP[i]-500
    rangeh<-ndat$BP[i]+500
    ndat$count[i]<-ndat%>%filter(chr==dchr,BP<=rangeh,BP>rangel)%>%count()%>%pull(n)
  }
  
  ndat <- ndat %>% mutate(gene = case_when
                          ((chr=="Pf3D7_10_v3")&(BP>1432498)&(BP<1434786)~"pfdblmsp2",
                            (chr=="Pf3D7_03_v3")&(BP>125687)&(BP<130112)~"PF3D7_0302300",
                            .default = "others"
                          ))
  
  p3<-ggplot(ndat,aes(x=count, fill=sig))+
    geom_histogram(aes(y=after_stat(density)),position="dodge")+
    xlab("SNP-density (1Kb vicinity)")+
    scale_fill_manual(values=c("blue","red"),name=expression(P<10^-5))+
    theme_clean(base_size=16)+
    theme(legend.position = c(0.7,0.8))
  p4<-ggplot(ndat%>%filter(sig==T),aes(as.factor(CHR),count))+
    geom_jitter(aes(color=gene), width=0.25)+
    xlab("Chromosome")+ylab("SNP-density (1Kb vicinity)")+
    scale_colour_colorblind()+
    geom_boxplot(width=0.6,color="grey",alpha=0.5, fill=NA,outlier.shape = NA)+
    theme_clean(base_size=16)+
    theme(legend.position = c(.3,.8))
  p4
  Reviewer1 <- p3 + p4+plot_annotation(tag_levels="A")
  Reviewer1
  ggsave("output/SNPdensity_reviewer1.pdf",Reviewer1,width=12,height=6)
  ggsave("output/SNPdensity_reviewer1.png",Reviewer1,width=12,height=6)
  
  
## iHS Calculation-------

# remove moi=2 inds
# load vcf file pf_moi1_final_MAFfilter_002.recode.vcf" 

## bash script:
## converting vcf to transposed haplotype format
# ./bcftools query -f '[%GT ]\n' -o pf_moi1_thap.txt -r Pf3D7_01_v3 pf_moi1_final_MAFfilter_002.recode.vcf
# substitue genotypes from 0(|/)0 to 0, 1/1 to 1, and ./. or 0/1 to .
# awk '{gsub(/0\/0|0\|0/, "0");gsub(/1\/1|1\|1/, "1");gsub(/0\/1|0\|1|\.\/\./, ".");}1' pf_moi1_thap.txt > pf_moi1_thap_ihsrecode.txt

# generate the map file
# ./bcftools query -f '%CHROM\.%POS %CHROM %POS %REF %ALT\n' -o pf_moi1_map.txt pf_moi1_final_MAFfilter_002.recode.vcf

# read in the genotype file to split them into chromosomes

allGen<-read.delim("moi1_ihs/pf_moi1_thap_ihsrecode.txt",header=F, sep=" ")
allGen<-allGen[,1:50]
mapFile<-read.delim("moi1_ihs/pf_moi1_map.txt",header=F, sep="")
mapFile$V2<-str_split_fixed(mapFile$V2,"_",3)[,2]
mapFile <- mapFile %>% filter(V2 != "API")
mapFile$V2<-as.numeric(mapFile$V2)
write.table(mapFile,"moi1_ihs/pf_moi1_map.txt",sep=" ",row.names=F,col.names=F,quote=F)
uniqueChrom <- unique(mapFile$V2)

for (chr in uniqueChrom){
  tempDat <- as.data.frame(t(as.matrix(allGen[mapFile$V2==chr,])))
  tempDat <- cbind(moi1_samples$V1,tempDat)
  write.table(tempDat, paste("moi1_ihs/","pf_moi1_hap_ihsrecode_",chr,".txt",sep=""),quote=F,col.names=F,row.names=F)
}


library(rehh)
moi1<-list()
for (i in 1:14){
  chr <- chrList[i]
  moi1[[i]]<-rehh::data2haplohh(hap_file = paste("moi1_ihs/pf_moi1_hap_ihsrecode_",chr,".txt",sep=""),
                                map_file = "moi1_ihs/pf_moi1_map.txt",
                                allele_coding = '01',
                                min_perc_geno.hap=90,
                                min_perc_geno.mrk =90,
                                min_maf=0.02,
                                chr.name = as.character(i))
  scan <- rehh::scan_hh(moi1[[i]],
                        polarized=F,
                        maxgap = 50000,
                        scalegap = 20000,
                        discard_integration_at_border = T)
  if (i == 1) {
    wgscan <- scan
  } else {
    wgscan <- rbind(wgscan, scan)
  }
  
  
}
# Computing genowide scan 

wgscan_ihs<-ihh2ihs(wgscan, freqbin = 1, verbose = F, standardize = T)
manhattanplot(wgscan_ihs)

wgsRegion<-calc_candidate_regions(wgscan_ihs,threshold = 5, 
                                  ignore_sign = T,
                                  pval = T,
                                  min_n_extr_mrk = 2,
                                  window_size=50000,overlap=10000)
debug(calc_candidate_regions)
undebug(calc_candidate_regions)
wgsRegion$CHR<-str_split_fixed(wgsRegion$CHR,"_",3)[,2]
wgsRegion$CHR<-as.numeric(wgsRegion$CHR)
write.csv(wgsRegion,"output/iHS_sig_region.csv",row.names=F)

distribplot(wgscan_ihs$ihs$IHS, xlab = "iHS")
distribplot(wgscan_ihs$ihs$IHS, 
            xlab = "iHS", 
            qqplot = TRUE)
library(tidyverse)
highMarkers <- wgscan_ihs$ihs %>% filter(LOGPVALUE>5) %>% mutate(cat=cut(POSITION,seq(0,3255610,163))) %>%group_by(CHR, cat) %>%slice_max(LOGPVALUE,n=1)
EFFinfo$chrNum <- as.numeric(str_split_fixed(EFFinfo$chr,"_",3)[,2])
highMarker <- highMarkers %>% left_join(EFFinfo %>% dplyr::select(pos_bp,chrNum,Effect,Amino_Acid_Change,Gene_Name),by=c("CHR"="chrNum","POSITION"="pos_bp"))
highMarker <- as.data.frame(highMarker)
rownames(highMarker)<-paste(highMarker$Gene_Name,highMarker$POSITION, sep="_")
wgscan_ihs$ihs$CHR<-str_split_fixed(wgscan_ihs$ihs$CHR,"_",3)[,2]
wgscan_ihs$ihs$CHR<-as.numeric(wgscan_ihs$ihs$CHR)
pdf("ihs_sig.pdf",width=10,height=4)
manhattanplot(wgscan_ihs,
              pval = TRUE,
              inset = 10,
              threshold = 5,
              mrk= highMarker,
              mrk.cex=2,
              mrk.lab.cex=2,
              cr = wgsRegion,
              main = "p-value of iHS")
dev.off()
wgsResults<- cbind(rownames(wgscan_ihs$ihs),wgscan_ihs$ihs)
colnames(wgsResults)[1]<-"snp_id"
wgsResults$P <- 1/10^(wgsResults$LOGPVALUE)
wgsResults$CHR<-as.numeric(wgsResults$CHR)

wgsSig <- wgsResults %>% filter(LOGPVALUE>5)
wgsSig <- wgsSig %>% left_join(EFFinfo, by=c("CHR"="chrNum","POSITION"="pos_bp"))
write.csv(wgsSig,"output/iHS_sig_sites.csv",row.names=F)

## iHS haplotype plots----------

EFFinfo$variant.id <- paste(EFFinfo$chr, EFFinfo$pos_bp,sep=".")
for (i in 1:nrow(wgsRegion)){
  #longChrName <- sprintf("Pf3D7_%02d_v3",wgsRegion$CHR[i])
  mark<- wgscan_ihs$ihs %>% filter(CHR==wgsRegion$CHR[i],POSITION>= wgsRegion[i,"START"],POSITION<= wgsRegion[i,"END"]) %>%slice_max(LOGPVALUE,n=1,na_rm=T)
  markName <- rownames(mark)[1]
  furcation <- calc_furcation(moi1[[wgsRegion$CHR[i]]],
                              mrk = markName)
  hapDat <- data.frame(iid =hap.names(moi1[[wgsRegion$CHR[i]]]) )
  hapDat <- hapDat %>% left_join(fws_metadata %>%select(iid,province))
  figpdf <- sprintf('iHS_region_%d_%s.pdf',i,markName)
  pdf(file = figpdf,width=8,height=16)
  par(mfrow = c(2, 1), mar = c(5.1, 6, 4.1, 6))
  plot(furcation,xlim = c(wgsRegion$START[i],wgsRegion$END[i]),hap.names = hapDat$province,
       main=EFFinfo$EFF[EFFinfo$variant.id==markName])
  haplen <- calc_haplen(furcation)
  plot(haplen,hap.names = hapDat$province,
       main=EFFinfo$EFF[EFFinfo$variant.id==markName])
  dev.off()
}


## iHS and isoRelate overlap---------

ihs_sig_snp <- wgscan_ihs$ihs %>% filter(LOGPVALUE>5)
ihs_sig_snp$snp_id <- rownames(ihs_sig_snp)
isorelate_sig_snp <- samples241_final_iR %>% filter(log10_pvalue > 5)
isorelate_sig_snp$snp_id <- gsub(":",".",isorelate_sig_snp$snp_id)


install.packages("VennDiagram")
library(VennDiagram)

venn.diagram(
  x = list(ihs_sig_snp$snp_id,isorelate_sig_snp$snp_id),
  category.names = c("iHS","IBD"),
  filename = "output/venn_ihs_IBD.png",
  output=T
  )

par(mfrow = c(2, 1), mar = c(5.1, 6, 4.1, 6))
manhattan(samples241_final_iR, col = c("black", "royalblue"), annotateTop= T, annotatePval=NULL,suggestiveline = F, genomewideline = -log10(1e-5), ylim= c(0,20), chr="CHR", bp="pos_bp", snp="snp_id", p="P" , main= "Significant IBD sharing")
manhattan(wgsResults, col = c("black", "royalblue"), annotateTop= T, suggestiveline = F, genomewideline = -log10(1e-5), ylim= c(0,15), chr="CHR", bp="POSITION", snp="snp_id", p="P" , main= "Significance of |iHS|")

samples241_final_iR$POSITION<- samples241_final_iR$pos_bp

ibdRegion<-calc_candidate_regions(samples241_final_iR,threshold = 5, 
                                  ignore_sign = T,
                                  pval = T,
                                  min_n_extr_mrk = 2,
                                  window_size=50000,overlap=10000)

write.csv(ibdRegion,"output/ibd_regions.csv",row.names=F)

pdf("Comparison_selection_results.pdf",width=10,height=8)
png("Comparison_selection_results.png",width=10,height=8,units = "in",res=300)
par(mfrow = c(2, 1), mar = c(3, 6, 1.5, 4))
manhattanplot(samples241_final_iR,
              pval = TRUE,
              inset = 10,
              threshold = 5,
              mrk= NULL,
              mrk.cex=2,
              mrk.lab.cex=2,
              cr = ibdRegion,
              ylab = expression(-log[10](italic(p))),
              xlab = "",
              main = expression(bold(paste("Significance of ",X[iR]))),
              mgp=c(2,0.5,0)
              )
mtext(substitute(paste(bold('A.'))),adj=0,cex=1.5,line=0.2)
manhattanplot(wgscan_ihs,
              pval = TRUE,
              inset = 10,
              threshold = 5,
              mrk= highMarker,
              mrk.cex=2,
              mrk.lab.cex=2,
              cr = wgsRegion,
              xlab = "Chromosome",
              main = "Significance of |iHS|",
              mgp=c(2,0.5,0))
mtext(substitute(paste(bold('B.'))),adj=0,cex=1.5,line=0.2)
dev.off()


