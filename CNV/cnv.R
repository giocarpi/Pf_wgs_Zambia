library(tidyverse)
pf3D7_featureCounts <- read_tsv("pf3D7_featureCounts.tsv",comment="#")
BasicInfo <- pf3D7_featureCounts[,1:6]

allSamples <- read_tsv("combined_columns.txt",comment="#")

sampName <- colnames(allSamples)

sampName1<- str_split(sampName,"/",simplify = T)[,6]
sampName2 <- str_split(sampName1,"_",simplify = T)[,1]

colnames(allSamples) <- sampName2

moi1SampleList <- read_csv("../fws_metadata_geo.csv")
moi1SampleList <- moi1SampleList %>%filter(moi==1)%>%select(iid)

moi1sample <- allSamples[,sampName2 %in% moi1SampleList$iid ]

finalmoi1 <- cbind(BasicInfo,moi1sample)

#finalmoi1 <- finalmoi1 %>%mutate(avDepth = Sample5027250/Length)
finalmoi1 <- finalmoi1 %>%rowwise()%>%mutate(loc = min(as.numeric(str_split_1(Start,";")),
                                                       as.numeric(str_split_1(End,";"))))
finalmoi1 <- finalmoi1 %>%rowwise()%>%mutate(ChrUniq = str_split_1(Chr,";")[1])
finalmoi1 <-finalmoi1 %>% ungroup() %>%filter(! ChrUniq %in% c("Pf3D7_API_v3","PfD7_MIT_v3"))

finalmoi1 <-finalmoi1 %>% mutate(Strand=NULL,Start=NULL, End = NULL)
finalmoi1 <-finalmoi1 %>% mutate(Chr = NULL)

newDat <- finalmoi1 %>% pivot_longer(cols=3:52, names_to="sample_id",values_to="coverage")

newDat <- newDat %>% mutate(avDepth = coverage/Length)
newDat <- newDat %>% ungroup() %>% group_by(Geneid) %>% mutate(allID = sum(avDepth))%>%filter(allID!=0)
newDat <- newDat %>% ungroup() %>% group_by(sample_id) %>% 
  mutate(medDepth = mean(avDepth),sdDepth = sd(avDepth))
newDat <- newDat %>% ungroup() %>% group_by(sample_id) %>% 
  mutate(medDepth = median(avDepth))
newDat <- newDat %>% ungroup() %>% mutate(allID=NULL, cn = avDepth/medDepth)
newDat <- newDat %>%filter(! ChrUniq %in% c("Pf3D7_API_v3","PfD7_MIT_v3"))
ggplot(newDat%>%filter(sample_id=="PU20"),aes(loc,cn))+
  geom_line()+
  facet_grid(ChrUniq~.)

ggplot(newDat%>%filter(sample_id=="PU20"),aes(cn))+geom_histogram(bins=40)

geneInfo <- read_csv("../IBD/Pf3D7_allSig.csv")
geneInfo <- geneInfo %>% distinct(Geneid,.keep_all = T) %>% select(Geneid, ProductDescription,InputsID)

geneStats <- geneStats[,1:8]

geneInfo <- geneInfo%>%left_join(geneStats)

geneInfo <- geneInfo[,1:3]
selSample <- newDat %>% inner_join(geneInfo, by="Geneid")

geneStats <- selSample%>% group_by(ChrUniq, loc,Geneid,Length,  InputsID, ProductDescription) %>% 
  summarize(meanCP = mean(cn),
            medianCP = median(cn),
            sdCP = sd(cn))%>%
  mutate( cvCP = sdCP/meanCP)  

write_csv(geneStats,"selected_gene_cnv.csv")

ggplot(newDat%>%filter(Geneid=="PF3D7_0523000"), 
       aes(cn,fct_reorder(Geneid,loc)))+
  geom_boxplot()+
  geom_vline(xintercept=1,linetype=2,col="red")+
  xlab("copy number")+
  ylab("")

ggplot(selSample%>%filter(ChrUniq=="Pf3D7_10_v3", !InputsID %in% c("Pf11-1", "LSA1")), 
       aes(cn,fct_reorder(InputsID,loc)))+
  geom_boxplot()+
  geom_vline(xintercept=1,linetype=2,col="red")+
  xlab("copy number")+
  ylab("")

ggplot(selSample%>%filter(ChrUniq=="Pf3D7_03_v3", !InputsID %in% c("Pf11-1", "LSA1")), 
       aes(cn,fct_reorder(InputsID,loc)))+
  geom_boxplot()+
  geom_vline(xintercept=1,linetype=2,col="red")+
  xlab("copy number")+
  ylab("")


ggplot(selSample%>%filter(ChrUniq=="Pf3D7_06_v3", !InputsID %in% c("Pf11-1", "LSA1")), 
       aes(cn,fct_reorder(InputsID,loc)))+
  geom_boxplot()+
  geom_vline(xintercept=1,linetype=2,col="red")+
  xlab("copy number")+
  ylab("")

ggplot(selSample%>%filter(ChrUniq %in% c("Pf3D7_08_v3","Pf3D7_12_v3"), !InputsID %in% c("Pf11-1", "LSA1")), 
       aes(cn,fct_reorder(InputsID,loc)))+
  geom_boxplot()+
  geom_vline(xintercept=1,linetype=2,col="red")+
  xlab("copy number")+
  ylab("")


library(scales)
ggplot(newDat%>%filter(Geneid=="PF3D7_1036300"),aes(cn))+geom_histogram(bins=40)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="b")+
  theme_classic()

ggplot(newDat%>%filter(Geneid=="PF3D7_1223400"),aes(cn))+geom_histogram(bins=40)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="b")+
  theme_classic()

ggplot(newDat%>%filter(Geneid=="PF3D7_1224000"),aes(cn))+geom_histogram(bins=40)+
  scale_x_log10()+annotation_logticks(sides="b")+
  theme_classic()

geneStats <- newDat%>% group_by(Geneid,Length, loc, ChrUniq) %>% 
  summarize(meanCP = mean(cn),
            medianCP = median(cn),
            sdCP = sd(cn))%>%
  mutate( cvCP = sdCP/meanCP)           
           
geneStats %>% ungroup() %>%slice_max(cvCP, n=20)
geneStats <- geneStats %>%mutate(Transcript_ID = str_c(Geneid,".1"))
eff <- read_csv("../processingCode/input/EffInfo.csv")

eff <- eff %>% distinct(Transcript_ID,.keep_all = T)

geneStats <- geneStats %>% arrange(desc(cvCP))


