#VERTAAL ALLES NOG KNULLO

#Libraries ####
library(tidyverse)
#Package voor het inscannen en presenteren van afbeeldingen
library(png)
library(grid)
library(gridExtra)
#For creating the DESeq2 object
library(DESeq2)
#For creating heatmaps
library(pheatmap)
#For transforming enterezid into gene identifiers
library("org.Hs.eg.db")
#For transforming enterezid -> GO terms -> gene activity
library(GO.db)
#For creating a GOHyperG dataset
library(GOstats)

#(Opdracht 2)Inscannen en tentoonstellen van afbeeldingen fastqc PerBaseQuality ####
#Inlezen van afbeeldingen die handmatig uit de fastqc bestanden zijn gehaald
SRR7866699_PBQ<-readPNG("~/temp_daur2/eindopdracht/afbeeldingen/SRR7866699_1_perbasequality.png") %>% as.raster() %>% rasterGrob()
SRR7866703_PBQ<-readPNG("~/temp_daur2/eindopdracht/afbeeldingen/SRR7866703_2_perbasequality.png") %>% as.raster() %>% rasterGrob()
#Arrangen van deze afbeeldingen, toevoegen van titel
grid.arrange(SRR7866699_PBQ, SRR7866703_PBQ, ncol=2, top=textGrob("Per base kwaliteitsdistributie voor hoogste kwaliteit read (links, SRR7866699_1) en laagste kwaliteit read (rechts, SRR7866703_2).", gp=gpar(fontsize=8, font=8)))

#Inlezen van abeeldingen
SRR7866699_1_PSQ<-readPNG("~/temp_daur2/eindopdracht/afbeeldingen/SRR7866699_1_persequencequality.png") %>% as.raster() %>% rasterGrob()
SRR7866703_2_PSQ<-readPNG("~/temp_daur2/eindopdracht/afbeeldingen/SRR7866703_2_persequencequality.png") %>% as.raster() %>% rasterGrob()

grid.arrange(SRR7866699_1_PSQ, SRR7866703_2_PSQ, ncol=2)

#(Opdracht 3) Genereren van een count table

#Creeëren van een directory
dir.create("~/temp_daur2/eindopdracht/counts")

#Selecteren van pad naar directories voor code
bam_dir<-"/home/daur2/rnaseq/rnaseq_onecut/bam/"
counts_dir<-"~/temp_daur2/eindopdracht/counts/"
bam_dir_regex<-regex(pattern="SRR7866[67][90][0349].bam$")
list.files(bam_dir, pattern = bam_dir_regex, full.names=TRUE)

#(Opdracht 3)Echt genereren van een count table ####
read_counts<-featureCounts(
  files=bam_files,
  annot.inbuilt = "hg38",
  useMetaFeatures = TRUE,
  strandSpecific = 1,
  countReadPairs = TRUE,
  isPairedEnd = TRUE,
  nthreads = 10
)
#Reading the count directory (was already pre-generated, otherwise would be from counts_dir)
read_counts<-readRDS("/home/daur2/rnaseq/rnaseq_onecut/counts/read_counts_OC2.rds")
#Control to see if this is the proper count-table
read_counts$targets

#(Opdracht 4)Making the DESeq2 library ####

#Extracting the count matrix
count_matrix<-read_counts$counts

#Loading in metadata
metadata<-read_csv("/home/daur2/rnaseq/rnaseq_onecut/onecut_sampledata_OC2.csv") %>% as.data.frame()

#Transforming the row-names so metadata can properly combine with the count matrix
rownames(metadata)<-paste0(metadata$Run, ".bam")

#Checking if the names are truly the same
colnames(count_matrix)==rownames(metadata)

#Creating a factor-value in metadata, to be used to distinguish the monsters in graphs
metadata<-metadata %>% mutate(Treatment=c("Bclxl", "Bclxl", "OC2", "OC2"))
metadata$Treatment<-metadata$Treatment %>% factor(levels=c("Bclxl", "OC2"))

#Making sure treatment is a proper factor
is.factor(metadata$Treatment)

#Creating the DESeq2 object
dds<-DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ Treatment
)

#(Opdracht 5)Performing a PCA analysis:PropOfVariance Barchart ####

#Normalising of DESeq2 library
dds_normalised<-rlog(dds)

#Creating a PCA (Assay transforms the dds into a single vector-like object, T creates a transposure, )
PCA<-dds_normalised %>% assay() %>% t() %>% prcomp()

#Creating a summary of the PCA
PCA_summary<-summary(PCA)$importance
PCA_summary_variance<-PCA_summary %>% t() %>% as.data.frame() %>% mutate(PC=colnames(PCA_summary))
PCA_summary_variance$`Proportion of Variance`



#Creating a bar chart to summarise the proportions of variance
PCA_summary_variance %>% ggplot(aes(x=PC, y=`Proportion of Variance`, fill=PC))+
  geom_col(colour="black")+
  labs(
    title="Proportion of variance of each principal component",
    x="Principal component",
    y="Proportion of variance"
  )+
  theme_classic()+
  theme(legend.position = "none")

#(Opdracht 5)Performing a PCA analysis:PCA-graph ####

#Sticking the metadata-results onto the PCA results
PCA_plotting<-cbind(metadata, PCA$x)

#Creating variables containing the exact PC proportion of variance percentages
PC1<-round(PCA_summary["Proportion of Variance","PC1"]*100, digits=1)
PC2<-round(PCA_summary["Proportion of Variance","PC2"]*100, digits=1)

#Creating a plot 
PCA_plotting %>% ggplot(aes(x=PC1, y=PC2))+
  geom_point(aes(colour=Treatment), size=3)+
  labs(
    title="Principal component analysis of fibroblasts treated by Bclxl or OC2",
    x=paste0("PC1 (",PC1,"%)"),
    y=paste0("PC2 (",PC2,"%)")
  )+
  theme_bw()

#(Opdracht 6)DGE analyse + volcano plot ####

#Performing DGE analysis
neuron_dge<-DESeq(dds)

#Storing the results with proper p and LFC value in a seperate RDS and inspecting the results
neuron_dge_results<-results(neuron_dge, alpha = 0.01, lfcThreshold = 1)
summary(neuron_dge_results)

#Filtering all NA padj values
neuron_dge_plotting<-neuron_dge_results %>% data.frame() %>% filter(!is.na(padj))
neuron_dge_plotting %>% head()

#Adding a variable which tracks if a gene has a padj<0.05 and a log>|1|
neuron_dge_plotting<-neuron_dge_plotting %>% mutate(
  Significance=if_else(padj<0.05 & abs(log2FoldChange)>1,"Significant","Not significant")
)

#seperating the total amount of up- and downregulated genes (Is not working so far, just took the amount of genes by looking at summary for now)
upregulated<-neuron_dge_results %>% as.data.frame() %>% filter(padj<0.01 & log2FoldChange>1) %>% nrow()
downregulated<-neuron_dge_results %>% as.data.frame() %>% filter(padj<0.01 & log2FoldChange< -1) %>% nrow()

summary(neuron_dge_results)

#Plotting in graph
neuron_dge_plotting %>% ggplot(aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(colour=Significance), size=0.7)+
  labs(
    title = "All genes showing a change in expression between Bclxl and OC2",
    subtitle = "Significance is defined as padj<0.05 and |LFC|>1",
    x="Log2 fold change in gene expression",
    y="-10log adjusted P.value"
  )+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_vline(xintercept = -1, linetype="dashed")+
  geom_hline(yintercept = -log10(0.01), linetype="dashed")+
  scale_colour_manual(values=c("Dark blue", "Dark orange"))+
  annotate("text", x=-5, y=170, label=c("Downregulated:\n\n",downregulated))+
  annotate("text", x=11, y=170, label=c("Upregulated:\n\n",upregulated))+
  theme_bw()

#(Opdracht 6)Heatmap  ####

#Obtaining the top 5 upregulated/downregulated genes, putting them into 1 list
neuron_dge_results_p_0.01<-neuron_dge_results %>% as.data.frame() %>% filter(padj<0.01)

bot_5_genes<-rownames(neuron_dge_results_p_0.01[order(neuron_dge_results_p_0.01$log2FoldChange)[1:5],])
top_5_genes<-rownames(neuron_dge_results_p_0.01[order(neuron_dge_results_p_0.01$log2FoldChange, decreasing = TRUE)[1:5],])
tb_10_genes<-c(top_5_genes, bot_5_genes)

#Obtaining the original data associated with these genes
count_values<-assay(dds)[tb_10_genes,]

#Replacing column names with treament
colnames(count_values)<-colData(dds)$Treatment

#Replacing row name entrez numbers with gene names
tb_10_genes_geneID<-data.frame(entrez=tb_10_genes)
tb_10_genes_geneID$Symbol<-mapIds(org.Hs.eg.db, 
       keys=tb_10_genes_geneID$entrez,
       column = "SYMBOL",
       keytype = "ENTREZID",
       multiVals = "first")
rownames(count_values)<-tb_10_genes_geneID$Symbol

#Generating a heatmap
pheatmap(count_values, scale="row", main = "Heatmap of top 5 up & downregulated genes")

#(Opdracht 7)Functie voor DGE analyse ####
#The function TRY 1 
goToSymbol<-function(GOTerm){
  GO_list<-select(org.Hs.eg.db,
                  keys=GOTerm,
                  columns = "SYMBOL",
                  keytype = "GO",
                  multiVals = "list")
  GO_list_vector<-GO_list$SYMBOL %>% as.vector()
  print(GO_list_vector)
  message("Is the output a vector ánd a character?")
  c(is.vector(GO_list_vector),is.character(GO_list_vector))
}

goToSymbol("GO:0036003")
#The function TRY 2
GOidToSymbol<-function(GOid){
  GO_list<-select(org.Hs.eg.db,
         keys=GOid,
         columns = "SYMBOL",
         keytype = "GO",
         multiVals = "list")
  GO_list$SYMBOL %>% as.vector()
}
#Run function
goToSymbol2(GOid = c("GO:0036002", "GO:0036003"))

#The function, now featuring controls!
GOidToSymbol<-function(GOid){
  if(str_detect(string = paste(GOid), pattern = "GO:[0-9]{7}$")){
    GO_list<-select(org.Hs.eg.db,
                    keys=GOid,
                    columns = "SYMBOL",
                    keytype = "GO",
                    multiVals = "list")
    GO_list$SYMBOL %>% as.vector()  
  }
  else("Input GO-term is not of the correct type: has to start 'GO:' followed by 7 numbers")
}

testfunction("GO:0033600")
#
#(Opdracht 8)GO term enrichment analysis ####

#Selecteren van upregulated, downregulated and all genes
neuron_dge_upregulated<-neuron_dge_results %>% data.frame() %>% filter(log2FoldChange>1 & padj<0.01) %>% rownames()
neuron_dge_downregulated<-neuron_dge_results %>% data.frame() %>% filter(log2FoldChange< -1 & padj<0.01) %>% rownames()
neuron_dge_all<-neuron_dge_results %>% data.frame() %>% rownames()

# ^ UPREGULATED ####
#GO-term analysis object creation
test_object<-new("GOHyperGParams",
    geneIds=neuron_dge_upregulated,
    universeGeneIds=neuron_dge_all,
    annotation="org.Hs.eg.db",
    ontology="BP",
    pvalueCutoff=1,
    testDirection="over")

#Hyper G test
neuron_dge_hyperG<-hyperGTest(test_object)
neuron_dge_hyperG_result<-summary(neuron_dge_hyperG)

#Padjust for the "BH" method, Bejmani hochberg which compensates for the false discovery rate
neuron_dge_hyperG_result$padj<-p.adjust(neuron_dge_hyperG_result$Pvalue, method="BH")

#Selecting gene sets larger than 5 but smaller than 500, to prevent too generic/too specific results from influencing the results
neuron_dge_hyperG_result<-neuron_dge_hyperG_result %>% filter(Count>5 & Count<500)

#Selecting the top 20 GO-terms, to prevent the graph from become unreadable
hyperG_top20<-neuron_dge_hyperG_result[order(neuron_dge_hyperG_result$padj)[1:20],]

#Creating a factor so that the graph goes the proper way
hyperG_top20$Term<-factor(hyperG_top20$Term,
                          levels = hyperG_top20$Term[
                            order(hyperG_top20$padj,decreasing = TRUE)
                          ])

#Plotting the p-values of the top 20
hyperG_top20 %>% ggplot(aes(x=Term, y=-log10(padj)))+
  geom_point()+
  coord_flip()+
  labs(
    title="Top 20 most upregulated \nG-terms",
    x="GO terms",
    y="-log10 adjusted P-value"
  )+
  theme_bw()

# ^ DOWNREGULATED ####

test_object_down<-new("GOHyperGParams",
                 geneIds=neuron_dge_downregulated,
                 universeGeneIds=neuron_dge_all,
                 annotation="org.Hs.eg.db",
                 ontology="BP",
                 pvalueCutoff=1,
                 testDirection="over")

#Hyper G test
neuron_dge_hyperG_down<-hyperGTest(test_object_down)
neuron_dge_hyperG_result_down<-summary(neuron_dge_hyperG_down)

#Padjust for the "BH" method, Bejmani hochberg which compensates for the false discovery rate
neuron_dge_hyperG_result_down$padj<-p.adjust(neuron_dge_hyperG_result_down$Pvalue, method="BH")

#Selecting gene sets larger than 5 but smaller than 500, to prevent too generic/too specific results from influencing the results
neuron_dge_hyperG_result_down<-neuron_dge_hyperG_result_down %>% filter(Count>5 & Count<500)

#Selecting the top 20 GO-terms, to prevent the graph from become unreadable
hyperG_top20_down<-neuron_dge_hyperG_result_down[order(neuron_dge_hyperG_result_down$padj)[1:20],]

#Making a factor so the plot goes in the proper direction
hyperG_top20_down$Term <- factor(hyperG_top20_down$Term,
                                 levels = hyperG_top20_down$Term[
                                   order(hyperG_top20_down$padj, decreasing = TRUE)
                                 ])

#Plotting the p-values of the top 20
hyperG_top20_down %>% ggplot(aes(x=Term, y=-log10(padj)))+
  geom_point()+
  coord_flip()+
  labs(
    title="Top 20 most downregulated \nG-terms",
    x="GO terms",
    y="-log10 adjusted P-value"
  )+
  theme_bw()
