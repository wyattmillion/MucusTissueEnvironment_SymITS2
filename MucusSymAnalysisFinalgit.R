#Mucus vs. Tissue Symbiont Community Analysis
library(tidyverse)
library(naniar)
library(miaViz)
library(cowplot)
library(plyr)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("phylosmith-package")
#devtools::install_github('schuyler-smith/phylosmith')
#install.packages("remotes")
remotes::install_github("Russel88/MicEco")
#remotes::install_github("schuyler-smith/phyloschuyler")
library(phyloseq)
library(phylosmith) #helps work with phyloseq objects
library(MicEco)
library(edgeR)
library(MCMC.OTU)
library(microbial)
library(ggpubr)
library(grid)
# Reading in the sampe metadata
info<-read.csv("info_300718.csv",header=T) #read in the file that links sample_ID to colony info (i.e. plot, ID, species, type of sample)
info$environ.type<-ifelse(grepl("Coral",info$Sample.ID),"tissue",
                                 ifelse(grepl("Mucus",info$Sample.ID),"mucus",
                                        info$environ.type))
info<-drop_na(info,plot..) #remove samples that are not part of the organized sampling design (should removes non coral species & negative controls)
#missing 1 mucus sample from coral 29
info$plot..<-as.factor(info$plot..)
info$environ.type<-as.factor(info$environ.type)

#####Sequencing summary
########
#reads per sample before Symportal
hist(info$no..of.reads)  #few low count samples, and few with >500,000 reads
#may want to remove samples with low read # (less than 10,000?)- this threshold would only remove 2 turf samples
#contains 2 negative control samples

summary(info$no..of.reads)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2816  116210  169374  223428  256836 1174688
#three lowest counts come from turf samples (18TAB.2, 7TAB and 16TAB) which are dramatically lower than the next lowest read sample and both negative controls
#will remove these low count samples
info<-info[-c(195,190,173),]
info$x<-1
summarise(info,mean(no..of.reads),sd(no..of.reads))
ggplot(info,aes(x=environ.type,y=no..of.reads))+
  geom_point(position = position_jitter(.1))+
  geom_boxplot(outlier.shape = 24,alpha=.7)

####ending this step with 196 samples which represent plot-associated samples from all environ types only
#######

#####Symportal summary
####### 
#should use post-med sequences for the comparison between environmental sources and coral sources
#this is because there are no ITS profiles for env. samples and post-med also reduced seq by abundance
#however, it is recommended to use the pre-med sequences for assessments of alpha diversity 

cseq<-read.delim("coral_mucus_sp_results/post_med_seqs/167_20210827_DBV_20210830T011943.seqs.absolute.abund_and_meta.txt",header=T) #use absolute abundance for this, can always recalculate relative abundances
#cpre<-read.csv("coral_mucus_sp_results/pre_med_seqs/pre_med_absolute_abundance_df.csv",header=T) #use absolute abundance for this, can always recalculate relative abundances
eseq<-read.delim("non_coral_mucus_sp_results/post_med_seqs/20210827T042640.seqs.absolute.abund_and_meta.txt",header=T)
#epre<-read.csv("non_coral_mucus_sp_results/pre_med_seqs/pre_med_absolute_abundance_df.csv",header=T)


#first remove the samples we identified in the first step
cseq<-cseq[cseq$sample_name %in% info$Sample.ID,]
eseq<-eseq[eseq$sample_name %in% info$Sample.ID,]

hist(cseq$post_med_absolute)
hist(eseq$post_med_absolute) #large discrepency in # of sequences per sample for coral vs non coral samples

test<-filter(eseq,post_med_absolute>500)

#will want to do outlier removal and normalization across sample
##Purging outliers (based on https://ryaneckert.github.io/website/ITS2_data/)
#remove low abundance (< 0.01%) sequences
#normalize once all samples are together in merge phyloseq object
cseq.purge = purgeOutliers(cseq, count.columns = 38:length(cseq), otu.cut = 0.0001, sampleZcut = -5) #reduces 781 DIVs to 208
eseq.purge = purgeOutliers(eseq, count.columns = 38:length(eseq), otu.cut = 0.0001, sampleZcut = -5) #reduces 1298 DIVs to 861


##3. Phyloseq analysis
#####

#statistical comparison (unifrac, bray curtis, , review methods of sym community differences)
#asking how similar or different are the sources of symbionts (within a plot) or if sources are more similar regardless of plot

#PERMANOVA (adonis2 in vegan) to test for differences between sources
#getting data into a phyloseq format

#You can merge phyloseq datasets together easily so keep env and coral seqs separate for now
cp_otu<-as.matrix(cseq.purge[,38:245])
c_otu<-as.matrix(cseq[,38:781]) #use this version to compare effect of removing outliers
ep_otu<-as.matrix(eseq.purge[,38:448])
e_otu<-as.matrix(eseq[,38:1298]) 


row.names(c_otu)<-cseq[,2] #adding sample names to the matrix
row.names(e_otu)<-eseq[,2]
row.names(cp_otu)<-cseq[,2] #adding sample names to the matrix
row.names(ep_otu)<-eseq[,2]


c_otu<-t(c_otu) #transpose so that OTU=rows
e_otu<-t(e_otu)
cp_otu<-t(cp_otu) #transpose so that OTU=rows
ep_otu<-t(ep_otu)


#now need a taxanomy table for both env and coral data
c_tax<-as.data.frame(row.names(c_otu))
e_tax<-as.data.frame(row.names(e_otu))
cp_tax<-as.data.frame(row.names(cp_otu))
ep_tax<-as.data.frame(row.names(ep_otu))


rownames(c_tax)<-c_tax$`row.names(c_otu)` #rownames of tax table have to match row names in OTU table
rownames(e_tax)<-e_tax$`row.names(e_otu)`
rownames(cp_tax)<-cp_tax$`row.names(cp_otu)` #rownames of tax table have to match row names in OTU table
rownames(ep_tax)<-ep_tax$`row.names(ep_otu)`


colnames(c_tax)<-c("DIV") #rename the column to DIV (unique section identifier from Symportal)
colnames(e_tax)<-c("DIV")
colnames(cp_tax)<-c("DIV") #rename the column to DIV (unique section identifier from Symportal)
colnames(ep_tax)<-c("DIV")


#adding another column that gives the Genus, taken from the DIV name (optional)
cp_tax$SymGenus[grepl("A",cp_tax$DIV)]<-"Symbiodinium"
cp_tax$SymGenus[grepl("B",cp_tax$DIV)]<-"Breviolum"
cp_tax$SymGenus[grepl("C",cp_tax$DIV)]<-"Cladocopium"
cp_tax$SymGenus[grepl("D",cp_tax$DIV)]<-"Durisdinium"
cp_tax$SymGenus[grepl("E",cp_tax$DIV)]<-"Effrenium"
cp_tax$SymGenus[grepl("G",cp_tax$DIV)]<-"Gerakladium"
cp_tax$SymGenus[grepl("F",cp_tax$DIV)]<-"Fugacium"
cp_tax$SymGenus[grepl("H",cp_tax$DIV)]<-"Clade H"
cp_tax$SymGenus[grepl("I",cp_tax$DIV)]<-"Clade I"

ep_tax$SymGenus[grepl("A",ep_tax$DIV)]<-"Symbiodinium"
ep_tax$SymGenus[grepl("B",ep_tax$DIV)]<-"Breviolum"
ep_tax$SymGenus[grepl("C",ep_tax$DIV)]<-"Cladocopium"
ep_tax$SymGenus[grepl("D",ep_tax$DIV)]<-"Durisdinium"
ep_tax$SymGenus[grepl("E",ep_tax$DIV)]<-"Effrenium"
ep_tax$SymGenus[grepl("G",ep_tax$DIV)]<-"Gerakladium"
ep_tax$SymGenus[grepl("F",ep_tax$DIV)]<-"Fugacium"
ep_tax$SymGenus[grepl("H",ep_tax$DIV)]<-"Clade H"
ep_tax$SymGenus[grepl("I",ep_tax$DIV)]<-"Clade I"

c_tax<-as.matrix(c_tax) #convert to a matrix for phyloseq
e_tax<-as.matrix(e_tax)
cp_tax<-as.matrix(cp_tax) #convert to a matrix for phyloseq
ep_tax<-as.matrix(ep_tax)


#finally, a sample table that gives metadata for the sample names
sample_data<-info[,c(1:6,8,11)] #getting just the data we want from info
row.names(sample_data)<-sample_data$Sample.ID #have  to have names that match the sample names in OUT table 

#seperating out the samples (135 for coral 102 for env)
csample_data<- sample_data %>% filter(environ.type=="mucus"|environ.type=="tissue")
esample_data<-sample_data %>% filter(environ.type!="mucus" & environ.type!="tissue")

##converting to phyloseq format
cPS = phyloseq(otu_table(c_otu, taxa_are_rows = TRUE),tax_table(c_tax),sample_data(csample_data)) #making th phyloseq object
ePS=phyloseq(otu_table(e_otu, taxa_are_rows = TRUE),tax_table(e_tax),sample_data(esample_data))
cpPS = phyloseq(otu_table(cp_otu, taxa_are_rows = TRUE),tax_table(cp_tax),sample_data(csample_data)) #making th phyloseq object
epPS=phyloseq(otu_table(ep_otu, taxa_are_rows = TRUE),tax_table(ep_tax),sample_data(esample_data))

#ePS.full=phyloseq(otu_table(e_otu.full, taxa_are_rows = TRUE),tax_table(e_tax.full),sample_data(esample_data))

####For pre-med data
########
cpre<-read.csv("coral_mucus_sp_results/pre_med_seqs/pre_med_absolute_abundance_df.csv",header=T) #use absolute abundance for this, can always recalculate relative abundances
epre<-read.csv("non_coral_mucus_sp_results/pre_med_seqs/pre_med_absolute_abundance_df.csv",header=T)

#row.names(cpre_otu)<-cpre[1:135,2]
#row.names(epre_otu)<-epre[c(1:71,73:81,83:102),2]

#cpre_otu<-t(cpre_otu)
#epre_otu<-t(epre_otu)

#now need a taxanomy table for both env and coral data
#cpre_tax<-as.data.frame(row.names(cpre_otu))
#epre_tax<-as.data.frame(row.names(epre_otu))

#rownames(cpre_tax)<-cpre_tax$`row.names(cpre_otu)`
#rownames(epre_tax)<-epre_tax$`row.names(epre_otu)`

#colnames(cpre_tax)<-c("DIV")
#colnames(epre_tax)<-c("DIV")

#convert to a matrix for phyloseq
#cpre_tax<-as.matrix(cpre_tax)
#epre_tax<-as.matrix(epre_tax)

#cprePS=phyloseq(otu_table(cpre_otu, taxa_are_rows = TRUE),tax_table(cpre_tax),sample_data(csample_data))
#eprePS=phyloseq(otu_table(epre_otu, taxa_are_rows = TRUE),tax_table(epre_tax),sample_data(esample_data))
######

#final phyloseq objects
cPS #no outlier removal: 744 taxa, 99 samples
cpPS #Purged outliers: 208 taxa, 99 samples

ePS # no outlier removal: 1261 taxa, 97 samples
epPS #purged outliers: 824 taxa, 97 samples 

#cprePS #135 samples, 143890 seq
#eprePS #100 samples, 15426 seq

# make both into an RDS and only run this once
saveRDS(cpPS,"cpPS.purged_3.6.23.RDS")
saveRDS(cPS,"cPS.unpurged_3.6.23.RDS")
saveRDS(epPS,"epPS.purged_3.6.23.RDS")
saveRDS(ePS,"ePS.unpurged_3.6.23.RDS")
#saveRDS(eprePS,"eprePS.RDS")
######

# read in RDS if you need to
cpPS <- readRDS("cpPS.purged_3.6.23.RDS")
epPS <- readRDS("epPS.purged_3.6.23.RDS")
cPS<-readRDS("cPS.unpurged_3.6.23.RDS")
ePS<-readRDS("ePS.unpurged_3.6.23.RDS")
#eprePS <- readRDS("eprePS.RDS")
#cprePS <- readRDS("cprePS.RDS")

#if you want to merge the coral and environmental phyloseq data together
allPS<-merge_phyloseq(cPS,ePS)
allpPS<-merge_phyloseq(cpPS,epPS)
#allprePS<-merge_phyloseq(cprePS,eprePS)


########## Look at Diversity ##########
#######################################
#Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset***
#total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)
#Shannon:Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).
#Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, λ obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of λ that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/λ) and the Gini–Simpson index (1 − λ).

#alpha diversity metrics must be run on raw, untrimmed data because many metrics rely on singletons
#however, even the full env. seq data set gives errors from phyloseq about singletons
# plot Shannon and Simpson diversity
plot_richness(allprePS, x="environ.type") + theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#even when including all taxa from premed sequences, still receive the warning of "no singletons, this is highly suspicious" and asks for untrimmed data

#prune step  to low abundance taxa,  worth doing for premed sequences but post med sequences already underwent consolidation

plot_bar(trimPS, "Sample", fill="DIV",facet_grid=~environ.type)



#######

### How to dealing with large differences in counts per sample
######
#1. Rarefying: probably not the best option as recommended by Phyloseq people and the fact that this is the Symportal output
#but if we want to try it, we can set a threshold and go from there
rarePS.T<-rarefy_even_depth(allPS,sample.size = 500, replace=T) #removed 37 samples and 150 OTUs
rarePS.F<-rarefy_even_depth(allPS,sample.size = 500, replace=F) #removes 37 samples and 133 OTUs

#2. Normalization techniques: the main thought is that the scale of coral vs noncoral samples will make any ordination misleading
#goal is to try and get things on same scale

#relative normalization (relative abundances rather than absolute abundance)
allPS.rel  = transform_sample_counts(allPS, function(x) x / sum(x) ) #transformed to relative abundance
allpPS.rel  = transform_sample_counts(allpPS, function(x) x / sum(x) ) #transformed to relative abundance

#allPS.rel = normalize(allPS, method="relative") #does the same thing as line above

####**best to just use a relative transformation
##Normalization using DESeq2 varianceStabilizingTransformation method
allPS.vst<-normalize(allPS,group = "environ.type" ,method="vst")
##Normalization using TMM
allPS.tmm<-normalize(allPS, method="TMM")
##Normalization with log2
allPS.log<-normalize(allPS, method = "log2")




#############################################
##################  PCoA  ###################
#############################################

#Bray Curtis PCoA All collapsed types with PCA plot

### Setting up distance matrix
#not sure which method to use, lets compare
dist = "bray"
ord_meths = c("PCoA","NMDS","MDS") #pick which methods to include
object=normalize(epPS,method = "relative")
type="environment only, purged, relative tranformation"
plist = llply(as.list(ord_meths), function(i, allPS, dist){
  ordi = ordinate(object, method=i, distance=dist)
  plot_ordination(object, ordi, type="sample", color="environ.type")
}, object, dist)

names(plist) <- ord_meths #add methods names for organization

#new dataframe containing PC componet values
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

#plotting the mutliple ordination plots
p = ggplot(pdataframe, aes(Axis_1, Axis_2, color=environ.type, shape=plot..))+
  geom_point(size=2)+
  facet_wrap(~method, scales="free")+
  labs(title=type)+
  scale_color_manual(values=c("lightblue", "goldenrod", "salmon","seagreen"))
p

#NMDS uses non-metric multidimensional scaling (converts distances to rank)
#MDS and PCoA use metric multidim scaling so these two will look more similar

#going with PCoA because I want to use actual distances not rank-based
target=normalize(cpPS,method="relative")
ps.ord <- ordinate(target,"PCoA",distance="bray")
ps.ord
#plotting the PCoA
quartz()
plot_ordination(target, ps.ord,type="sample", color ="environ.type")+ #color by site to see the grouping we care about here
  scale_color_manual(values=c("lightblue", "goldenrod", "salmon","seagreen"))+
  scale_shape_manual(values=c(15,16,17,18,20,8,11,0,6))+
  geom_point(alpha=0.5, size=5)+
  stat_ellipse(aes(linetype=environ.type))+
  theme_cowplot()

#####################   Stats   ########################
######
library(vegan)
#install.packages("remotes")
#remotes::install_github("Jtrachsel/funfuns")
library(funfuns)
library(edgeR)
library("remotes")

#For all data (coral + environmenal)
#not sure if I should use the normalized or raw data...
seq.ps <- data.frame(normalize(allpPS,method="relative")@otu_table) #sample and DIV info
  #removing rows with all zero
  #seq.ps$count<-apply(seq.ps,1,function(x){sum(x>0)})
  #seq.ps<-filter(seq.ps,count>0)
  #seq.ps<-seq.ps[,1:235]
seq.ps<-t(seq.ps)
samdf.ps <- data.frame(normalize(allpPS,method="relative")@sam_data) #samples with metadata
dist.ps <- vegdist(seq.ps)#computes dissimilarity indeces
#dist.ps<-phyloseq::distance(seq.ps, method = "bray")
bet.ps <- betadisper(dist.ps,as.vector(samdf.ps$environ.type)) #assess betadiversity metrics

#Effects of environment type on betadiversity (does dispersion within an eviron.type differ)
anova(bet.ps) 
permutest(bet.ps,pairwise=TRUE,permutations=999, model = "full") #if beta diversity is significant by environtype, use this to see where the differences occur
plot(bet.ps)

#Permanova to see how different communities are from eachother
results<-adonis(seq.ps ~ environ.type*plot.., data=samdf.ps, permutations=999,method="bray")
results$aov.tab

pairwise.adonis(seq.ps, factors=samdf.ps$environ.type,sim.method = "bray", p.adjust.m = "BH", permutations=999)
#shows every environ group is significantly different from each other


##How do the other data sets compare? Just verifying type of transformation and filtering do not impact general conculsion
######
PS<-normalize(ePS,method="log")
seq.ps <- data.frame(PS@otu_table) #sample and DIV info
seq.ps<-t(seq.ps)
samdf.ps <- data.frame(PS@sam_data) #samples with metadata
dist.ps <- vegdist(seq.ps)#computes dissimilarity indeces
#dist.ps<-phyloseq::distance(seq.ps, method = "bray")
bet.ps <- betadisper(dist.ps,as.vector(samdf.ps$environ.type)) #assess betadiversity metrics
anova(bet.ps) 
permutest(bet.ps,pairwise=TRUE,permutations=999) #if beta diversity is significant by environtype, use this to see where the differences occur
plot(bet.ps)
results<-adonis(seq.ps ~ environ.type*plot.., data=samdf.ps, permutations=999)
results$aov.tab

pairwise.adonis(seq.ps, factors=samdf.ps$environ.type, permutations=999)

#######end of verification

####### end of stats on seq-based analysis

####Visualizing sequences shared across sources
################

evenn<-inner_join(eseq.purge,filter(info, !grepl('coral|mucus', Sample.no.))[,c(4,8)], by=c("sample.1"="Sample.ID"))
evenn<-pivot_longer(evenn,38:448, names_to = "Taxa",values_to = "counts")
evenn<-evenn[,-(3:37)]

cvenn<-inner_join(cseq.purge,filter(info, grepl('coral|mucus', Sample.no.))[,c(4,8)], by=c("sample.1"="Sample.ID"))
cvenn<-pivot_longer(cvenn,38:245, names_to = "Taxa",values_to = "counts" )
cvenn<-cvenn[,-(3:37)]

allvenn<-rbind(cvenn,evenn)

#make a list of vectors that can be fed into UpSet
tissue<-unique(allvenn %>% filter(environ.type == "tissue" & counts>0) %>% pull(4))
mucus<-unique(allvenn %>% filter(environ.type == "mucus" & counts>0) %>% pull(4))
seawater<-unique(allvenn %>% filter(environ.type == "seawater" & counts>0) %>% pull(4))
turf<-unique(allvenn %>% filter(environ.type == "turf algae" & counts>0) %>% pull(4))
sedclose<-unique(allvenn %>% filter(environ.type == "sed close" & counts>0) %>% pull(4))
sedfar<-unique(allvenn %>% filter(environ.type == "sed far" & counts>0) %>% pull(4))

vlist<-list(tissue,mucus,seawater,turf, sedclose, sedfar)
names(vlist)<-c("tissue","mucus","seawater","turf algae", "sediment close", "sediment far")

#time to make the UpSet plot
library(UpSetR)

intersections <- which(rowSums(df) >= 10)

# Keep only those intersections in the data frame
df_filtered <- df[, intersections]

# Plot the UpSet chart
upset(df_filtered, sets = names(df_filtered))
upset(fromList(vlist), sets= c(names(vlist)),
      order.by = c("freq"),
      set_size.show = F,
      mb.ratio = c(0.6, 0.4),
      queries =  list(list(query = intersects, params = list("tissue"), color = "rosybrown", active = T),
           list(query = intersects, params = list("mucus"), color = "paleturquoise4", active = T),
           list(query = intersects, params = list("turf algae"), color = "seagreen", active = T),
           list(query = intersects, params = list("seawater"), color = "lightblue", active = T),
           list(query = intersects,params = list("sediment far"),color="salmon", active = T),
           list(query = intersects,params = list("mucus","tissue"),color="purple", active = T)))

##identifying the sequences in sml or tissue but not in both
length(tissue) #205
length(mucus) #192

length(unique(c(tissue,mucus))) #208 unique sequences

notshared<-tissue[!tissue %in% mucus]
notshared<-c(notshared,mucus[!mucus %in% tissue])
length(notshared)

relcpPS<-as.data.frame(otu_table(normalize(cpPS,method = "relative")))
relcpPS<-relcpPS[notshared,]
relcpPS$ITS2<-row.names(relcpPS)
test<-pivot_longer(relcpPS,cols = 1:99,names_to = 'sample',values_to = "relative abundance")
hist(test$`relative abundance`)
#the majority of the sequences occur at below 1% if occuring at all in any sample
max(test$`relative abundance`) #highest abundance sequence is 3.7% in a single mucus sample

######

#this step will group low abundance DIVs into "other X" based on genus
sw <- epPS
sw.df <- psmelt(sw)

unique(sw.df$SymGenus)
Sym<-sw.df %>% filter(SymGenus=="Symbiodinium") %>%
  group_by(DIV) %>%
  mutate(Rel.Abundance=sum(Abundance)/283353)
Sym$DIV[Sym$Rel.Abundance < 0.01] <- "Other A"

Clad<-sw.df %>% filter(SymGenus=="Cladocopium") %>%
  group_by(DIV) %>%
  mutate(Rel.Abundance=sum(Abundance)/318431)  
Clad$DIV[Clad$Rel.Abundance < 0.01] <- "Other C"

Dur<-sw.df %>% filter(SymGenus=="Durisdinium") %>%
  group_by(DIV) %>%
  mutate(Rel.Abundance=sum(Abundance)/39616)  
Dur$DIV[Dur$Rel.Abundance < 0.01] <- "Other D"

Ger<-sw.df %>% filter(SymGenus=="Gerakladium") %>%
  group_by(DIV) %>%
  mutate(Rel.Abundance=sum(Abundance)/6301)  
Ger$DIV[Ger$Rel.Abundance < 0.01] <- "Other G"

Fug<-sw.df %>% filter(SymGenus=="Fugacium") %>%
  group_by(DIV) %>%
  mutate(Rel.Abundance=sum(Abundance)/381)  
Fug$DIV[Fug$Rel.Abundance < 0.01] <- "Other F"

ClH.<-sw.df %>% filter(SymGenus=="Clade H") %>%
  group_by(DIV) %>%
  mutate(Rel.Abundance=sum(Abundance)/93)  
ClH.$DIV[ClH.$Rel.Abundance < 0.01] <- "Other H"

ClI.<-sw.df %>% filter(SymGenus=="Clade I") %>%
  group_by(DIV) %>%
  mutate(Rel.Abundance=sum(Abundance)/66)  
ClI.$DIV[ClI.$Rel.Abundance < 0.01] <- "Other I"

sw.df<-rbind(Sym,Clad,Dur,Ger,Fug,ClH.,ClI.)

re.sw<-normalize(ePS,method="relative")
glom <- tax_glom(re.sw, taxrank = 'DIV')
re.sw.df <- psmelt(glom)

re.sw.df$Sample_Taxa<-paste(re.sw.df$Sample,"_",re.sw.df$OTU)
sw.df$Sample_Taxa<-paste(sw.df$Sample,"_",sw.df$OTU)
re.sw.df<-re.sw.df %>% 
  select(-12) %>%
  inner_join(sw.df[,c(12,15)],by="Sample_Taxa")
#glom <- tax_glom(sw, taxrank = 'DIV')
#glom # should list # taxa as # phyla
#data <- psmelt(glom) # create dataframe from phyloseq object
#data$DIV <- as.character(data$DIV) #convert to character

#reorder DIV legend so that colors follow according to Genus
re.sw.df$DIV<-factor(re.sw.df$DIV,levels=c("Other A","Other C","Other D","Other G",
  "A1","A1bk","A1j",   "A4at","X4488_A","X57113_A","X61508_A","X61510_A", "X61511_A", "X6322_A",
  "C1","C15","C15h","C1b","C3","C3z","C41","C41b","C91e","C91f","C91i","C92a","X1490_C",
  "G3b", "X57128_G","X57129_G","X57130_G","X57131_G","X57132_G", "X57133_G","X57867_G", "X58181_G","X1869317_G" ,
  "F3a","F5a","F5b","F5c","X63307_F",
  "D1","D1c","D3","D3a","D4","D4f","D5","D6", "X56861_D",
  "X61691_H","X61697_H",
  "X42459_I","X57807_I","X58909_I","X738829_I","X1869401_I"))                                
                         
sw<-filter(re.sw.df,environ.type=="seawater")
sw$Sample<-paste(sw$plot..,"_",sw$Sample)

#make the figures first then stitch together after
pSW <- ggplot(data=sw, aes(x=Sample, y=Abundance, fill=DIV))+
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual("DIV",values=c(gray.colors(4),paletteer_dynamic("cartography::blue.pal", 10),
                                     paletteer_dynamic("cartography::orange.pal", 13),
                                     paletteer_dynamic("cartography::red.pal", 10),
                                     paletteer_dynamic("cartography::green.pal",5),
                                     paletteer_dynamic("cartography::purple.pal",9),
                                     paletteer_dynamic("cartography::pink.pal",2),
                                     paletteer_dynamic("cartography::sand.pal",5)))+
  theme(legend.position="none",axis.text.x=element_blank()) + guides(fill=guide_legend(ncol = 2))

sc<-filter(re.sw.df,environ.type=="sed close")
sc$Sample<-paste(sc$plot..,"_",sc$Sample)
pClose<-ggplot(data=sc, aes(x=Sample, y=Abundance, fill=DIV))+
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual("DIV",values=c(gray.colors(4),paletteer_dynamic("cartography::blue.pal", 10),
                                   paletteer_dynamic("cartography::orange.pal", 13),
                                   paletteer_dynamic("cartography::red.pal", 10),
                                   paletteer_dynamic("cartography::green.pal",5),
                                   paletteer_dynamic("cartography::purple.pal",9),
                                   paletteer_dynamic("cartography::pink.pal",2),
                                   paletteer_dynamic("cartography::sand.pal",5)))+
  theme(legend.position="none",axis.text.x=element_blank())+
  guides(fill=guide_legend(ncol = 2))

sf<-filter(re.sw.df,environ.type=="sed far")
sf$Sample<-paste(sf$plot..,"_",sf$Sample)
  pFar<-ggplot(data=sf, aes(x=Sample, y=Abundance, fill=DIV))+
    geom_bar(aes(), stat="identity", position="stack") +
    scale_fill_manual("DIV",values=c(gray.colors(4),paletteer_dynamic("cartography::blue.pal", 10),
                                     paletteer_dynamic("cartography::orange.pal", 13),
                                     paletteer_dynamic("cartography::red.pal", 10),
                                     paletteer_dynamic("cartography::green.pal",5),
                                     paletteer_dynamic("cartography::purple.pal",9),
                                     paletteer_dynamic("cartography::pink.pal",2),
                                     paletteer_dynamic("cartography::sand.pal",5)))+
    theme(legend.position="none",axis.text.x=element_blank())+
  guides(fill=guide_legend(ncol = 2))
    
  
ta<-filter(re.sw.df,environ.type=="turf algae")
ta$Sample<-paste(ta$plot..,"_",ta$Sample)
pTurf<-ggplot(data=ta, aes(x=Sample, y=Abundance, fill=DIV))+
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual("DIV",values=c(gray.colors(4),paletteer_dynamic("cartography::blue.pal", 10),
                                   paletteer_dynamic("cartography::orange.pal", 13),
                                   paletteer_dynamic("cartography::red.pal", 10),
                                   paletteer_dynamic("cartography::green.pal",5),
                                   paletteer_dynamic("cartography::purple.pal",9),
                                   paletteer_dynamic("cartography::pink.pal",2),
                                   paletteer_dynamic("cartography::sand.pal",5)))+
  theme(legend.position="none",axis.text.x=element_blank()) + guides(fill=guide_legend(ncol = 2))

pSed<-plot_grid(pClose,pFar,labels = c("b","c"))

quartz()
plot_grid(pSW,pSed,pTurf,ncol = 1, labels = c("a","","d"))
ggarrange(pSW,pClose,pFar,pTurf,nrow = 3, common.legend = T,legend = "right")

##Finding the number of sequences per genera in each sample
sf_no0<-filter(sf,Abundance>0)
sc_no0<-filter(sc,Abundance>0)
sw_no0<-filter(sw,Abundance>0)
ta_no0<-filter(ta,Abundance>0)

length(unique(sw_no0$DIV)) #41
length(unique(sw_no0$OTU)) #446
  length(unique(sw_no0$OTU)[(grepl("A",unique(sw_no0$OTU)))]) #96
  length(unique(sw_no0$OTU)[(grepl("C",unique(sw_no0$OTU)))]) #184
  length(unique(sw_no0$OTU)[(grepl("D",unique(sw_no0$OTU)))]) #127
  length(unique(sw_no0$OTU)[(grepl("F",unique(sw_no0$OTU)))]) #0
  length(unique(sw_no0$OTU)[(grepl("G",unique(sw_no0$OTU)))]) #38
  length(unique(sw_no0$OTU)[(grepl("H",unique(sw_no0$OTU)))]) #1
  length(unique(sw_no0$OTU)[(grepl("I",unique(sw_no0$OTU)))]) #0
  
length(unique(sc_no0$DIV)) #38
length(unique(sc_no0$OTU)) #181
  length(unique(sc_no0$OTU)[(grepl("A",unique(sc_no0$OTU)))]) #74
  length(unique(sc_no0$OTU)[(grepl("C",unique(sc_no0$OTU)))]) #71
  length(unique(sc_no0$OTU)[(grepl("D",unique(sc_no0$OTU)))]) #25
  length(unique(sc_no0$OTU)[(grepl("F",unique(sc_no0$OTU)))]) #3
  length(unique(sc_no0$OTU)[(grepl("G",unique(sc_no0$OTU)))]) #1
  length(unique(sc_no0$OTU)[(grepl("H",unique(sc_no0$OTU)))]) #2
  length(unique(sc_no0$OTU)[(grepl("I",unique(sc_no0$OTU)))]) #5
  
length(unique(sf_no0$DIV)) #41
length(unique(sf_no0$OTU)) #243
  length(unique(sf_no0$OTU)[(grepl("A",unique(sf_no0$OTU)))]) #58
  length(unique(sf_no0$OTU)[(grepl("C",unique(sf_no0$OTU)))]) #114
  length(unique(sf_no0$OTU)[(grepl("D",unique(sf_no0$OTU)))]) #53
  length(unique(sf_no0$OTU)[(grepl("F",unique(sf_no0$OTU)))]) #4
  length(unique(sf_no0$OTU)[(grepl("G",unique(sf_no0$OTU)))]) #12
  length(unique(sf_no0$OTU)[(grepl("H",unique(sf_no0$OTU)))]) #1
  length(unique(sf_no0$OTU)[(grepl("I",unique(sf_no0$OTU)))]) #1
  
length(unique(ta_no0$DIV)) #32
length(unique(ta_no0$OTU)) #333
  length(unique(ta_no0$OTU)[(grepl("A",unique(ta_no0$OTU)))]) #121
  length(unique(ta_no0$OTU)[(grepl("C",unique(ta_no0$OTU)))]) #137
  length(unique(ta_no0$OTU)[(grepl("D",unique(ta_no0$OTU)))]) #75
  length(unique(ta_no0$OTU)[(grepl("F",unique(ta_no0$OTU)))]) #0
  length(unique(ta_no0$OTU)[(grepl("G",unique(ta_no0$OTU)))]) #0
  length(unique(ta_no0$OTU)[(grepl("H",unique(ta_no0$OTU)))]) #0
  length(unique(ta_no0$OTU)[(grepl("I",unique(ta_no0$OTU)))]) #0

unique(sf$OTU)

unique(sc$OTU) %in% unique(ta$OTU)

# boxplots
dat %>% ggplot(aes(fill=DIV, y=Abundance, x=Sample)) + 
  geom_bar(position="stack", stat="identity") +
  labs(y="ITS2 Type Profile Relative Abundance")+
  xlab(NULL)

plot_bar(dat, fill= "SymGenus")+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position = "right",legend.text = element_text(size = 4))

#5.  Comparison of mucus to tissue symbiont communities using ITS2 type profiles

its<-read.delim("coral_mucus_sp_results/its2_type_profiles/167_20210827_DBV_20210830T011943.profiles.absolute.abund_and_meta.txt")

its_meta<-its[1:6,] #metadata for ITS2 type profiles
its_meta<-as.data.frame(t(its_meta),header=T)

its2<-its[6:143,] #ITS2 type profile abundances per coral sample
colnames(its2)<-c("UID","Sample.ID",its2[1,3:57])
its2<-its2[-which(its2$Sample.ID==""),] #remove rows not associated with a sample 
its2<-its2[which(its2$Sample.ID %in% info$Sample.ID),]
library(funrar)
its2 <- its2 %>% mutate_at(colnames(its2[3:57]), as.numeric)

#calculating relative abundances
rits<-make_relative(as.matrix(its2[,3:57]))
rits2<-cbind(its2[1:2],rits)


cinfo<-filter(info, grepl('coral|mucus', Sample.no.)) #remove samples that aren't from the mucus or coral, just trying to focus on this comparison
its2<-inner_join(cinfo,its2)
rits2<-inner_join(cinfo,rits2)

its2<-filter(its2,Sample.no.!="coral 29") #removing the coral sample that doesn't have a coresponding mucus sample
rits2<-filter(rits2,Sample.no.!="coral 29")

its2$coral.genus<-ifelse(its2$coral.genus=="","Unknown",its2$coral.genus)
rits2$coral.genus<-ifelse(rits2$coral.genus=="","Unknown",rits2$coral.genus)
rits2$coral.genus<-ifelse(rits2$coral.genus=="Potites","Porites",rits2$coral.genus)

long_rits<-rits2 %>% pivot_longer(14:68,names_to = "Type_Pro", values_to = "counts")

#Investigating the relationship between species and ITS2 type profile
genera<-unique(long_rits$coral.genus)
results_list<-list()
for (genus in genera) {
  a<-unique(filter(long_rits,coral.genus == genus,counts>.5)$Type_Pro)
  results_list[[genus]]<-a
}
print(results_list)

result_df <- results_list %>%
  enframe(name = "genus", value = "Type_Pro") %>%
  unnest(Type_Pro)

# Add a column to indicate presence (1) for the dot plot
result_df <- result_df %>%
  mutate(presence = 1)

type_pro_uniqueness <- result_df %>%
  group_by(Type_Pro) %>%
  dplyr::summarise(unique_count = n()) %>%
  mutate(is_unique = ifelse(unique_count == 1, "unique", "non-unique"))

# Merge this information back into the result_df
result_df <- result_df %>%
  left_join(type_pro_uniqueness, by = "Type_Pro")
# Print the data frame to check the structure
print(result_df)

ggplot(result_df, aes(x = genus, y = Type_Pro)) +
  geom_point(aes(size = presence, color = is_unique), shape = 16) + # Color points based on uniqueness
  scale_size_continuous(range = c(3, 3)) +       # Ensure all dots are the same size
  scale_color_manual(values = c("unique" = "red", "non-unique" = "black")) + # Color mapping
  labs(x = "Genus", y = "ITS2 Type Profile", title = "Presence of ITS2 Type Profile in Each Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#3. plotting 


long_rits<-long_rits %>% 
  mutate(Genus=case_when(startsWith(Type_Pro, "A") ~ "Symbiodinium",
                         startsWith(Type_Pro, "C")~ "Cladocopium",
                         startsWith(Type_Pro, "D")~"Durisdinium"))

long_rits<- long_rits %>%
  separate(col = Sample.no., into = c("Source", "Col_Num"), sep = " ") #Split Sample number into type of sample (coral vs mucus) and colony number

long_rits$Sample<-paste(long_rits$Col_Num,"_",long_rits$Source)


library(paletteer)
long_rits<- long_rits %>% arrange(coral.genus) #reorder rows so coral/mucus samples are together per colony
long_rits$plot..<-as.factor(long_rits$plot..)

#ordering so that samples are grouped by dominant type
long_rits$DomType<-ifelse(long_rits$counts > .5,long_rits$Type_Pro,NA)
for (i in unique(long_rits$Sample)) {
  t<-unique(na.omit(filter(long_rits,Sample==i)$DomType))
  long_rits$DomType<-ifelse(long_rits$Sample==i,t,long_rits$DomType)
}

coralist<-c("Acropora","Goniastrea","Platygyra","Echinopora","Galaxea","Favia","Fungia","Ctenactis","Porites","Xenia","Pavona","Pocillopora","Millepora", "Stylophora","Unknown")
long_rits$coral.genus<- factor(long_rits$coral.genus,levels = coralist)

long_rits<-long_rits %>% arrange((coral.genus)) %>% 
  mutate(
    Sample = factor(Sample, levels = unique(Sample)))
    
p1<-long_rits %>% filter(environ.type =="tissue") %>% 
  ggplot(aes(fill=Type_Pro, y=counts, x=Sample,group=counts)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label=coral.genus),y=.25,alpha=.1, angle=90,size=2)+
  labs(y="ITS2 Type Profile Relative Abundance")+
  xlab(NULL)+
  scale_fill_manual("Genus",values=c(hsv(.6, .6,seq(.4,1,length.out = 11)),
                                      hsv(.5, .6,seq(.4,1,length.out = 10)),
                                       hsv(0.05, .6,seq(.6,1,length.out = 12)),
                                       hsv(.1, .6,seq(.2,1,length.out = 13)),
                                       hsv(0.3, .6,seq(.4,.75,length.out = 9))))+
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1),
                      legend.position = "none",legend.text = element_text(size = 4))
legend<-get_legend(p1)

#now we want to order mucus samples to match tissue samples
ord<-levels(filter(long_rits,environ.type =="tissue")$Sample)
cord <- ord[grep("coral",ord)]

#setting the mucus order to mirror the tissue order
#if ordering samples by Sym type
mord1 <- c("3 _ mucus",  "39 _ mucus", "34 _ mucus","19 _ mucus","20 _ mucus","23 _ mucus", 
          "8 _ mucus" ,"38 _ mucus","4 _ mucus","40 _ mucus","42 _ mucus","14 _ mucus","16 _ mucus",
          "18 _ mucus","24 _ mucus","26 _ mucus","31 _ mucus","32 _ mucus", "35 _ mucus","5 _ mucus",
          "6 _ mucus","10 _ mucus", "27 _ mucus","46 _ mucus","47 _ mucus","1 _ mucus","21 _ mucus","36 _ mucus",  
          "48 _ mucus", "30 _ mucus","41 _ mucus" ,"7 _ mucus",  "15 _ mucus","25 _ mucus","13 _ mucus", "33 _ mucus",
          "43 _ mucus",  "49 _ mucus","50 _ mucus","22 _ mucus","17 _ mucus","12 _ mucus","45 _ mucus",
          "2 _ mucus",  "37 _ mucus", "11 _ mucus" , "28 _ mucus", "44 _ mucus","9 _ mucus")
#if ordering samples by genus
mord2<-c("1 _ mucus" , "5 _ mucus",  "14 _ mucus", "16 _ mucus", "21 _ mucus" ,"26 _ mucus",
         "31 _ mucus" ,"32 _ mucus", "38 _ mucus" ,"40 _ mucus", "42 _ mucus", "46 _ mucus",
         "10 _ mucus", "27 _ mucus", "36 _ mucus", "47 _ mucus", "4 _ mucus" , "18 _ mucus",
         "8 _ mucus" , "24 _ mucus", "35 _ mucus", "6 _ mucus" , "33 _ mucus", "43 _ mucus",
         "13 _ mucus" ,"50 _ mucus", "25 _ mucus", "7 _ mucus" , "15 _ mucus", "30 _ mucus",
         "34 _ mucus", "41 _ mucus", "3 _ mucus",  "39 _ mucus" ,"23 _ mucus", "9 _ mucus" ,
         "12 _ mucus", "28 _ mucus", "44 _ mucus", "17 _ mucus" ,"22 _ mucus", "2 _ mucus" ,
         "11 _ mucus", "37 _ mucus", "19 _ mucus", "20 _ mucus" ,"45 _ mucus", "48 _ mucus","49 _ mucus")

p2<-long_rits %>% filter(environ.type =="mucus") %>%
  mutate(Sample=factor(Sample, levels = mord2)) %>%
  ggplot(aes(fill=Type_Pro, y=counts, x=Sample,group=counts)) + 
  geom_bar(position="stack", stat="identity") +
  labs(y="ITS2 Type Profile Relative Abundance")+
  geom_text(aes(label=coral.genus),y=.25,alpha=.1, angle=90,size=2)+
  xlab(NULL)+
  scale_fill_manual("Genus",values=c(hsv(.6, .6,seq(.4,1,length.out = 11)),
                                     hsv(.5, .6,seq(.4,1,length.out = 10)),
                                     hsv(0.05, .6,seq(.6,1,length.out = 12)),
                                     hsv(.1, .6,seq(.2,1,length.out = 13)),
                                     hsv(0.3, .6,seq(.4,.75,length.out = 9))))+
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1),
                     legend.position = "none",legend.text = element_text(size = 4))

#Making combined figure
quartz()
ggarrange(p1+ rremove("ylab"),
          p2+ rremove("ylab"),
          ncol = 1, nrow = 2)%>%
  annotate_figure(left = textGrob("ITS2 Type Profile Relative Abundance", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                  bottom = textGrob("Sample", gp = gpar(cex = 1.3)))
as_ggplot(legend)

##for this plot, we want to show the DIV bar plot for colonies that don't have the same type-profile between tissue and mucus
coralDom<-long_rits %>% filter(Source=="coral", counts>.5)
mucusDom<-long_rits %>% filter(Source=="mucus", counts>.5)

nonsim<-coralDom[which(coralDom$DomType!=mucusDom$DomType),2]


cw <- cpPS
cw.df <- psmelt(cw)


#plotting the DIV absensce/precense per sample
adj<-cp_otu

adj2<-ifelse(adj>0,1,0)
library(reshape2)
adf<-melt(adj2)

colnames(adf)<-c("Taxa","Sample","Presence")
adf2<-filter(adf,str_detect(Sample,"Mucus"))
detach(package:plyr)
DIVs<-adf2 %>% group_by(Taxa) %>% summarise(NofOccur=sum(Presence))
DIVs<-DIVs %>% arrange(desc(NofOccur))
DIVs$Taxa<-as.character(DIVs$Taxa)
DIVord<-DIVs$Taxa

#want to set the order of the Taxa so most abundance taxa are at bottom
adf2$Taxa<-factor(adf2$Taxa, levels = DIVord)

#coloring DIVs but genus
adf2<-adf2 %>% mutate(SymGenus=case_when(grepl("A",adf2$Taxa)~"Symbiodinium",
                                          grepl("B",adf2$Taxa)~"Breviolum",
                                          grepl("C",adf2$Taxa)~"Cladocopium",
                                          grepl("D",adf2$Taxa)~"Durisdinium",
                                          grepl("E",adf2$Taxa)~"Effrenium",
                                          grepl("G",adf2$Taxa)~"Gerakladium",
                                          grepl("F",adf2$Taxa)~"Fugacium",
                                          grepl("H",adf2$Taxa)~"Clade H",
                                          grepl("I",adf2$Taxa)~"Clade I"))
adf2$SymGenus<-as.factor(adf2$SymGenus)
quartz()
ggplot(subset(adf2, Presence == 1), aes(x = Sample, y = Taxa,fill=Taxa))+
  geom_tile()+ylim(levels(adf2$Taxa))+
  scale_fill_manual("SymGenus",values=c(hsv(.58, .7,seq(.45,1,length.out = 48)),
                                        hsv(.5, .7,seq(.5,1,length.out = 111)),
                                        hsv(0, .7,seq(.5,1,length.out = 49))))+
  theme_classic2()+theme(legend.position = "none",
                   #axis.title.x = element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   #axis.title.y = element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),)
  


long_reldat[6818:13635,] %>% 
  #filter(Col_Num!="29") %>% 
  ggplot(aes(fill=Type_Pro, y=counts, x=Sample,group=counts)) + 
  geom_bar(position="stack", stat="identity") +
  labs(subtitle="ITS2 Profile Type",y="ITS2 Type Profile Relative Abundance")+
  xlab(NULL)+
  +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1),legend.position = "none")


mucus<-long_reldat %>% filter(Source=="mucus") %>% 
  ggplot(aes(fill=Type_Pro, y=counts, x=Col_Num,group=counts)) + 
  geom_bar(position="stack", stat="identity") +
  labs(subtitle="ITS2 Profile Type",y="ITS2 Type Profile Relative Abundance")+
  xlab(NULL)+
  scale_fill_manual("Genus",values=c(paletteer_dynamic("cartography::blue.pal", 13),
                                     paletteer_dynamic("cartography::orange.pal", 20),
                                     paletteer_dynamic("cartography::red.pal", 16),
                                     paletteer_dynamic("cartography::green.pal",9),
                                     gray.colors(59)))+
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1),legend.position = "none")

ggarrange(coral,mucus,ncol=1, nrow = 2,common.legend = F)

long_dat %>% filter(Source=="coral") %>% 
  ggplot(aes(fill=Type_Pro, y=counts, x=Col_Num)) + 
  geom_bar(position="stack", stat="identity") +
  labs(subtitle="ITS2 Profile Type",y="reads mapping to ITS2 Type Profile")+
  xlab(NULL)+
  #scale_fill_manual(values=c("steelblue3","chartreuse4","darkorchid3","orangered"))+
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))










#Comparing communities of different sources (turf, sediment, water column, coral, mucus )
#1. seq files have 700+ seq from each sample 
#if we want to remove low abundance sequences it looks like other papers used a variety of methods
#for now we will just move forward with all sequences maintained (jump to line 51)
rownames(dat)<-dat$sample_name
dat2<-dat[,3:142892]
rownames(reldat)<-reldat$sample_name
reldat2<-reldat[,3:142892]

tad<-as.data.frame(t(dat2))
tadler<-as.data.frame(t(reldat2))

#1.1 removing types with low abundance in a certain number of samples
tadler$low<-apply(tad,1,function(x){sum(x<=.0001)})
tadler0<-tadler[-which(tadler$low>60),] #60 is 44.4% of 135. So we are removing types that are very low abundance (.01%) in 44% of samples
#results in 35 high abundance types

#1.2 another paper (Howells et al 2020, Mol Eco) says type profiles at <1% proportional abundance within samples were excluded from downstream analyses.
tadler %>% replace_with_na_all(condition = ~.x < .01) #replaces values less than .01 with NA
tad[rowSums(is.na(tadler)) != ncol(tadler), ] #removes rows (type profiles) that have NA in every column (sample)

#1.3 another paper (Eckert et al. 2020 Fron.MicroBio) says it showed only the most abundant (>0.01% of all reads) sequences are displayed (n = ) 
#our n=135, so only sequences that represent at least .01% of all reads.
tad<-tad %>% mutate(
  sum=rowSums(.,na.rm=T)
)
sum(tad$sum)*.001
tad0<-tad[-which(tad$sum<18343),] #produces 62 types that have at least .1% of reads
ITSgood<-rownames(tad0)
tadler0<-tadler[rownames(tadler) %in% ITSgood,]

dat2<-as.data.frame(t(tad0))
dat2$Sample.ID<-rownames(dat2)

reldat2<-as.data.frame(t(tadler0))
reldat2$Sample.ID<-rownames(reldat2)