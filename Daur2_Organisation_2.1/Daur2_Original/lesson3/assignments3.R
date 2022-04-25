library(Rsubread)
library(tidyverse)
library(DESeq2)
library(pheatmap)
#Making the DEseq2 data set ####
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
#Performing a differential gene analysis ####
airway_dge<-DESeq(dds)
airway_dge_results<-results(airway_dge)
airway_dge_results
summary(airway_dge_results)

airway_dge_results_2<-results(airway_dge, lfcThreshold = 1, alpha=0.05) %>% summary()

#Mijn data lijkt vekeerd om?? als ik treated en untreated als factor omdraai lijkt het wel goed te zijn maar
#weet nu niet welke goed en welke slecht is.

#Graphing ####

#Selecting significant genes
sign_genes<-airway_dge_results[which(airway_dge_results$padj < 0.05),]

#Getting the top gene
topGene<-sign_genes[which.max(sign_genes$log2FoldChange),]
topGene_name<-rownames(topGene)
#Getting all the data from the top gene
geneCounts<-plotCounts(dds, gene=topGene_name, intgroup = c("treatment"), returnData=TRUE)
#Plotting it
geneCounts %>% ggplot(aes(x=treatment, y=count, colour=treatment))+
  scale_y_log10()+
  geom_point(position=position_jitter(width=0.1,height=0))+
  labs(
    title=paste0("Gen ", topGene_name),
    y="log count"
  )

#Selecting the bottom gene
botGene<-sign_genes[which.min(sign_genes$log2FoldChange),]
botGene_name<-rownames(botGene)
#Getting al the data from the bot gene
geneCounts_bot<-plotCounts(dds, gene=botGene_name, intgroup=c("treatment"), returnData=TRUE)

geneCounts_bot %>% ggplot(aes(x=treatment, y=count, colour=treatment))+
  scale_y_log10()+
  geom_point(position=position_jitter(width=0.1, height=0))+
  labs(
    title=paste0("Gen ", botGene_name),
    y="Log count"
  )

#Volcano graph ####
#Removing all p=NA data (p=NA means that there is no difference in regulation at all)
airway_dge_plotting<-data.frame(airway_dge_results) %>% filter(!is.na(padj))
#Create column which states if it differs significantly, and by how much
airway_dge_plotting_signif<-airway_dge_plotting %>% mutate(
  signif=if_else(padj<0.05,"padj < 0.05", "Not significant")
)
#Plot that shit
airway_dge_plotting_signif %>% ggplot(aes(x=log2FoldChange, y=-log10(padj), colour=signif))+
  geom_point()+
  labs(
    title="Gene change",
    y="-10 log adjusted p-value",
    x="Log 2-fold change"
  )+
  theme_bw()+
  scale_colour_manual(values=c("grey", "darkblue"), name="Significance")+
  annotate("text",x=topGene$log2FoldChange, y=-log10(topGene$padj)*0.8, label=topGene_name, colour="Red")

#Do it yourself volcano
airway_dge_plotting_LFC<-airway_dge_plotting %>% mutate(
  signif=if_else(padj<0.01,"padj < 0.01", "Not significant"),
  LFC=if_else(log2FoldChange>1,"LFC>1", "LFC<1"),
  signif_LFC=if_else(padj<0.01 & abs(log2FoldChange)>1,"Significant","Not significant")
)

airway_dge_plotting_LFC %>% ggplot(aes(x=log2FoldChange, y=-log10(padj), colour=signif_LFC))+
  geom_point()+
  labs(
    title="Signif + LFC",
    x="Log 2-fold change",
    y="-Log10 of adjusted P"
  )+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_vline(xintercept = -1, linetype="dashed")+
  geom_hline(yintercept = -log10(0.01), linetype="dashed")+
  scale_colour_manual(values=c("darkgrey", "darkred"), name="Significance+LFC")+
  theme_bw()

airway_dge_plotting %>% mutate(
  signif_LFC=if_else(padj<0.01 & log2FoldChange>1,"1","2")
)

#Heat map ####

top_10_genes<-rownames(airway_dge_results[order(airway_dge_results$padj)[1:10],])
count_values<-assay(dds)[top_10_genes,]
colnames(count_values)<-colData(dds)$treatment
pheatmap(count_values,show_rownames = TRUE, scale="row")

bot_10_genes<-rownames(airway_dge_results[order(airway_dge_results$padj, decreasing = TRUE)[1:10],])
count_values_bot<-assay(dds)[bot_10_genes,]
colnames(count_values_bot)<-colData(dds)$treatment
pheatmap(count_values_bot,show_rownames = TRUE, scale="row")
