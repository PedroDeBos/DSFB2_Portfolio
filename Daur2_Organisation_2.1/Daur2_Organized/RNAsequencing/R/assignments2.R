library(Rsubread)
library(tidyverse)
library(DESeq2)
library(pheatmap)
#Making a read-counts table ####
##Reden: laat zien relatief hoe vaak genen worden geactiveerd. Als de experimentele conditie de activatie van een
##gen verhoogd, zullen de genen vaker worden geactiveerd. Met dit kan je dat zien (en als input voor iets toekomstigs gebruiken)


#Uitvoeren van een count
bam_dir<-"/home/daur2/rnaseq/rnaseq_airway/bam"
counts_dir<-"~/temp_daur2/lesson2/counts/"
bam_files<-list.files(bam_dir, pattern = ".bam$", full.names = TRUE)
featureCounts(
  files=bam_files,
  annot.inbuilt = "hg38",
  useMetaFeatures = TRUE,
  strandSpecific = 0,
  isPairedEnd = TRUE,
  countReadPairs = TRUE,
  nthreads = 10
)

#Readin RDS attempt 1 ####
read_counts<-readRDS("/home/daur2/rnaseq/rnaseq_airway/counts/read_counts.rds")
#Inspecting lists
str(read_counts)
#Making it a tibble
read_counts_tb<-read_counts[[4]] %>% t() %>% as_tibble()
#Adding names
names<-list.files(bam_dir, pattern = ".bam$", full.names = FALSE)
read_counts_tb<-read_counts_tb %>% mutate(names)
#Renaming important headers and converting to numeric
read_counts_tb<-read_counts_tb %>% select(names, V1, V2)
read_counts_tb$V1<-read_counts_tb$V1 %>% as.numeric()
read_counts_tb$V2<-read_counts_tb$V2 %>% as.numeric()
read_counts_tb<-read_counts_tb %>% rename(assigned="V1", unassigned_unmapped="V2")
#Fuck it what I'm doing is wayy to complicated, I'll just look at the answer
#Reading RDS attempt 2 ####
read_counts<-readRDS("/home/daur2/rnaseq/rnaseq_airway/counts/read_counts.rds")
count_stats<-read_counts$stat

rownames(count_stats) <- count_stats$Status
count_stats$Status <- NULL
#god this little part makes it SO much easier holy fuck

names<-list.files(bam_dir, pattern = ".bam$", full.names = FALSE)

count_stats_sum<-count_stats %>% 
  t %>% 
  as_tibble() %>%
  mutate(names) %>% 
  mutate(total=colSums(count_stats)) %>% 
  mutate(perc_assigned=Assigned/total*100) %>%
  select(names, Assigned, total, perc_assigned)

count_stats_sum %>% ggplot(aes(x=names, y=perc_assigned, fill=names))+
  geom_col(colour="black")+
  labs(
    title="Percentage assigned DNA",
    x="",
    y="Percentage assigned (%)"
  )+
  theme_classic()+
  theme(axis.text.x=element_blank())

#Deseq 2, analysing RNA data ####
BiocManager::install("DESeq2")

#What we need for a DESeq2 data set:
##Count matrix,
##Dataframe with metadata
##Design formula

#Making count matrix
count_matrix<-read_counts$counts

#Obtaining meta data (how is this done in a real lab?) and rediting it to a dataframe
metadata<-read.csv("/home/daur2/rnaseq/rnaseq_airway/airway_sampledata.csv")
metadata<-as.data.frame(metadata)
rownames(metadata)<-paste0(metadata$Run, ".bam")
colnames(count_matrix) == rownames(metadata)

#Creating a factor (DO UNTREATED BEFORE TREATED, the program compares 1 to 2)
metadata<-metadata %>% mutate(treatment=str_replace(dex, "trt", "treated"))
metadata$treatment<-metadata$treatment %>% factor(levels=c("untreated", "treated"))

#Combining this all into a DESeq
dds<-DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ treatment
)

#Normalising RNA-seq count data ####
quantile(count_matrix[,"SRR1039509.bam"])
dds_normalized<-rlog(dds)

#Setting up a PCA ####
pca<-dds_normalized %>% assay() %>% t() %>% prcomp(.)
pca_summary<-summary(pca)$importance

#PCA voor PC1 en PC2
pca_plotting<-cbind(metadata,pca$x)

PCA_1<-round(pca_summary["Proportion of Variance", "PC1"]*100, digits = 1)
PCA_2<-round(pca_summary["Proportion of Variance", "PC2"]*100, digits = 1)

pca_plotting %>% ggplot()+
  geom_point(aes(x=PC1, y=PC2, color=treatment, shape=cell_line), size=4)+
  ggtitle("PCA voor aiway study")+
  xlab(paste0('PC 1(',PCA_1,"%)"))+
  ylab(paste0("PC 2(",PCA_2,"%)"))+
  theme_bw()

#PCA voor PC3 en PC4
pca_plotting<-cbind(metadata,pca$x)

PCA_3<-round(pca_summary["Proportion of Variance", "PC3"]*100, digits = 1)
PCA_4<-round(pca_summary["Proportion of Variance", "PC4"]*100, digits = 1)

pca_plotting %>% ggplot()+
  geom_point(aes(x=PC3, y=PC4, color=treatment, shape=cell_line), size=4)+
  ggtitle("PCA voor aiway study")+
  xlab(paste0('PC 1(',PCA_3,"%)"))+
  ylab(paste0("PC 2(",PCA_4,"%)"))+
  theme_bw()


#Setting up a heatmap ####
dds_normalized_matrix<-assay(dds_normalized)

airway_cor<-cor(dds_normalized_matrix)

pheatmap(airway_cor, annotation=metadata["treatment"], cluster_rows = TRUE, cluster_cols = TRUE)
?pheatmap()
