library(tidyverse)
library(Rsubread)
library(DESeq2)
library(pheatmap)
library("org.Hs.eg.db")
library(dplyr)
# Formatieve opdracht ####
#d)
#Index maken
hg38_ref<-"/home/daur2/rnaseq/hg38_genome/GRCh38.primary_assembly.genome.fa"
hg38_index<-"~/temp_daur2/formatieve_opdracht1/hg38_index"
buildindex(basename = hg38_index,
           reference = hg38_ref,
           gappedIndex = FALSE,
           indexSplit = FALSE)

#Alignment maken
##directories voorbereiden
fastq_dir<-"/home/daur2/rnaseq/rnaseq_ipsc/fastq/"
hg38_index_2<-"/home/daur2/rnaseq/hg38_index/hg38_index"
bam_dir<-"~/temp_daur2/formatieve_opdracht1/bam_dir"
samples<-list.files(fastq_dir, pattern="_[12].fastq.gz") %>% str_remove(pattern="_[12].fastq.gz") %>% unique()

alignments_statistics_ipsc<-align(
  index=hg38_index_2,
  readfile1=paste0(fastq_dir,samples,"_1.fastq.gz"),
  readfile2=paste0(fastq_dir,samples,"_2.fastq.gz"),
  type="rna", input_format = "gzFASTQ", output_format = "BAM",
  output_file = paste0(bam_dir,samples,".bam"),
  
  unique = TRUE,
  
  nthreads = 10
)
saveRDS(alignments_statistics_ipsc, file=paste0(bam_dir,"alignment_statistics.rds"))

#e)
#Data inlezen
alignment_statistics<-read_rds("/home/daur2/rnaseq/rnaseq_ipsc/bam/alignment_statistics.rds")
#Als tibble opslaan
alignment_statistics<-alignment_statistics %>% t() %>% as_tibble()
#Toevoegen van monster nummers
alignment_statistics<-alignment_statistics %>% mutate(monster=paste0("SRR",seq(786687,786694)))
#Volgorde veranderen
alignment_statistics<-alignment_statistics %>% select(monster, Total_fragments:Indels )

#Summary maken
alignment_statistics_sum<-alignment_statistics %>% group_by(monster) %>% summarise(
  perc_mapped=Mapped_fragments/Total_fragments*100,
  perc_unmapped=100-perc_mapped
)
alignment_statistics_pivot<-alignment_statistics_sum %>% pivot_longer(cols=c(perc_mapped, perc_unmapped), names_to = "map_status", values_to = "percentage")

#GGplot maken
alignment_statistics_pivot %>% ggplot(aes(x=monster, y=percentage, fill=map_status))+
  geom_col(colour="black")+
  theme_classic()+
  labs(
    title="Percentage mapped fragments experiment IPSC",
    x="Monster",
    y="Percentage (%)"
  )+
  theme(axis.text.x=element_text(angle=90))

#Making van count directory ####

bam_dir<-"/home/daur2/rnaseq/rnaseq_ipsc/bam/"
counts_dir<-"~/temp_daur2/formatieve_opdracht1/counts_dir/"
bam_files<-list.files(bam_dir, pattern=".bam$", full.names = TRUE)
featureCounts(
  files=bam_files,
  annot.inbuilt = "hg38",
  useMetaFeatures = TRUE,
  strandSpecific = 1,
  isPairedEnd = TRUE,
  nthreads = 10
)
#Reading the count directory (was already pre-generated, otherwise would be from counts_dir)
count_ipsc<-readRDS("/home/daur2/rnaseq/rnaseq_ipsc/counts/read_counts.rds")
count_matrix_ipsc<-count_ipsc$counts

#2b making a graph from $stat ####
count_ipsc_stats<-count_ipsc$stat

rownames(count_ipsc_stats)<-count_ipsc_stats$Status
count_ipsc_stats$Status<-NULL

count_ipsc_stats_sum<-count_ipsc_stats %>% t() %>% as_tibble() %>%
  mutate(bamfile=colnames(count_ipsc_stats)) %>% 
  mutate(total=colSums(count_ipsc_stats)) %>% 
  mutate(perc_bound=Assigned/total*100) %>%
  select(bamfile,total, perc_bound)

count_ipsc_stats_sum %>% dplyr::select(bamfile,total, perc_bound)

count_ipsc_stats_sum %>% ggplot(aes(x=bamfile, y=perc_bound, fill=bamfile))+
  geom_col(colour="black")+
  theme_classic()+
  theme(axis.text.x=element_blank())+
  labs(
    x="Bamfile",
    y="Percentage bound (%)",
    title="Percentage alligned RNA in ipsc"
  )

#Making a DSeq2 library ####
#Obtaining count matrix
count_ipsc<-readRDS("/home/daur2/rnaseq/rnaseq_ipsc/counts/read_counts.rds")
count_matrix_ipsc<-count_ipsc$counts
colnames(count_matrix_ipsc)

#Obtaining metadata
metadata_ipsc<-read_csv("/home/daur2/rnaseq/rnaseq_ipsc/ipsc_sampledata.csv")
metadata_ipsc<-as.data.frame(metadata_ipsc)
rownames(metadata_ipsc)<-paste0(metadata_ipsc$Run,".bam")

#Checking of row/colnames are the same
rownames(metadata_ipsc) == colnames(count_matrix_ipsc)

#Creating factor in metadata
metadata_ipsc$Cell_type
is.factor(metadata_ipsc$Cell_type)
metadata_ipsc$Cell_type<-metadata_ipsc$Cell_type %>% str_replace("Skin derived fibroblast", "skin_derived_fibroblast")
metadata_ipsc$Cell_type<-metadata_ipsc$Cell_type %>% factor(levels=c("skin_derived_fibroblast", "iPSC"))
is.factor(metadata_ipsc$Cell_type)

#Finally, creating the library
dds_ipsc<-DESeqDataSetFromMatrix(
  countData = count_matrix_ipsc,
  colData = metadata_ipsc,
  design = ~Cell_type
)

dds_ipsc
#Normalising it
dds_ipsc_normalised<-rlog(dds_ipsc)
dds_ipsc_normalised

#Maken van PCA ####
pca_ipsc<-dds_ipsc_normalised %>% assay() %>% t() %>% prcomp()
pca_summary<-summary(pca_ipsc)$importance
pca_ipsc$x

#plotting klaarzetten
pca_ipsc_plotting<-cbind(metadata_ipsc,pca_ipsc$x)
#Variabelen 1 en 2 klaarzetten
PCA_1<-round(pca_summary["Proportion of Variance","PC1"]*100,digits=1)
PCA_2<-round(pca_summary["Proportion of Variance","PC2"]*100,digits=1)
#In een ggplot
pca_ipsc_plotting %>% ggplot()+
  geom_point(aes(x=PC1, y=PC2, colour=Cell_type, shape=source_name), size=2)+
  labs(
    x=paste0("PC1 (",PCA_1,"%)"),
    y=paste0("PC2 (",PCA_2,"%)"),
    title="PCA tussen fibroblasten en iPSC"
  )+
  theme_bw()

#Heatmap maken ####
dds_ipsc_normalised_matrix<-assay(dds_ipsc_normalised)
ipsc_cor<-cor(dds_ipsc_normalised_matrix)
pheatmap(ipsc_cor, annotation = metadata_ipsc["Cell_type"])
?pheatmap()

#Les 3: DGE analyse ####
dds_ipsc_normalised
ipsc_dge<-DESeq(dds_ipsc)
ipsc_results<-results(ipsc_dge, alpha=0.05, lfcThreshold = 1)
summary(ipsc_results)
#LFC > 1 = 3222 (14%)
#LFC < -1 = 2510 (11%)
#in vergelijking met de luchtweg studie:
#LFC > 1 = 68 (0.32%)
#LFC < 1 = 36 (0.17%)
#Er is veeeeel meer significant verschil dan bij de luchtwegen

#les 3:
#Les 3: Volcano plot ####
ipsc_results_volcano<-data.frame(ipsc_results) %>% filter(!is.na(padj))

ipsc_results_volcano_2<-ipsc_results_volcano %>% mutate(
  significance=if_else(padj<0.05,"significant","not significant"),
  log2foldchange=if_else(abs(log2FoldChange)>1,"L2FC>1","L2FC<1"),
  signif_LFC=if_else(padj<0.05 & abs(log2FoldChange)>1,"Significant","Not significant")
)

ipsc_results_volcano_2 %>% ggplot(aes(x=log2FoldChange,y=-log10(padj),colour=signif_LFC))+
  geom_point()+
  labs(
    title="Gene expression IPSC",
    x="Log 2 fold change",
    y="-10log adjusted p-value"
  )+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_vline(xintercept = -1, linetype="dashed")+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  scale_colour_manual(values=c("darkgrey", "darkblue"), name="Significance+LFC")+
  theme_bw()

#Les 3: Heatmap top 15 ####
top_15_ipsc_names<-rownames(ipsc_results[order(ipsc_results$padj)[1:15],])
count_values_ipsc<-assay(dds_ipsc)[top_15_ipsc_names,]
colnames(count_values_ipsc)<-colData(dds_ipsc)$source_name
pheatmap(count_values_ipsc, show_rownames = TRUE, scale="row")


#Les 4: Heatmap als 3c, maar dan met gensymbool i.p.v. Entrez ####

#Code voor pheatmap
top_15_ipsc_names<-rownames(ipsc_results[order(ipsc_results$padj)[1:15],])
count_values_ipsc<-assay(dds_ipsc)[top_15_ipsc_names,]
colnames(count_values_ipsc)<-colData(dds_ipsc)$source_name
pheatmap(count_values_ipsc, show_rownames = TRUE, scale="row")

#Gen-symbolen tevoorschijn toveren
#Eerst de Entrez-waardes opslaan als data frame
top_15_ipsc_genidentifiers<-rownames(count_values_ipsc)
top_15_ipsc_genidentifiers<-top_15_ipsc_genidentifiers %>% data.frame()
#Met map-IDS de symbols krijgen
top_15_symbols<-mapIds(org.Hs.eg.db,
       keys = top_15_ipsc_genidentifiers$.,
       column = "SYMBOL",
       keytype = "ENTREZID",
       multiVals = "first") %>% data.frame()
rownames(count_values_ipsc)<-top_15_symbols$.

#Nogmaals heatmap maken
pheatmap(count_values_ipsc, show_rownames = TRUE, scale="row")


#Les 4: Top 20 genen opgereguleerd ####

#upregulated genen opslaan
ipsc_gene_upregulated<-ipsc_results %>% data.frame() %>% filter(log2FoldChange > 1, padj < 0.01) %>% rownames()

#ALLE genen opslaan
ipsc_gene_total<-ipsc_results %>% data.frame() %>% rownames()

#Data als object voor analyse opslaan
ipsc_testobject_upregulated<-new("GOHyperGParams",
                 geneIds=ipsc_gene_upregulated,
                 universeGeneIds=ipsc_gene_total,
                 annotation="org.Hs.eg.db",
                 ontology="BP",
                 pvalueCutoff=1,
                 testDirection="over"
)

#Analyse uitvoeren
ipsc_gotermanalysis_upregulated<-hyperGTest(ipsc_testobject_upregulated)

#Resultaten opslaan als summary
ipsc_goterm_sum_upregulated<-summary(ipsc_gotermanalysis_upregulated)

#Adjusted P-value toevoegen
ipsc_goterm_sum_upregulated$padj<-p.adjust(ipsc_goterm_sum_upregulated$Pvalue, method = "BH")

#Filteren genen kleiner dan 5, groter dan 500
ipsc_goterm_sum_upregulated<-ipsc_goterm_sum_upregulated %>% filter(Count>5 & Count<500)

#Selecteren van top 20 genen
ipsc_goterm_sum_upregulated_top20<-ipsc_goterm_sum_upregulated[order(ipsc_goterm_sum_upregulated$padj)[1:20],]

#Factor maken van "Term" en daarop sorteren
ipsc_goterm_sum_upregulated_top20$Term<-factor(
  ipsc_goterm_sum_upregulated_top20$Term,
  levels=ipsc_goterm_sum_upregulated_top20$Term[order(ipsc_goterm_sum_upregulated_top20$padj, decreasing = TRUE)]
)
is.factor(ipsc_goterm_sum_upregulated_top20$Term)

#Plotten
ipsc_goterm_sum_upregulated_top20 %>% ggplot(aes(x=Term, y=-log10(padj)))+
  geom_point()+
  coord_flip()+
  labs(
    title="Top 20 upregulated genes in iPSC's \n compared to fibroblasts",
    y="-log10 adjusted P-value",
    x="GO-terms"
  )+
  theme_bw()

#Les 4: Top 20 downregulated genes ####
#downregulated genen opslaan
ipsc_gene_downregulated<-ipsc_results %>% data.frame() %>% filter(log2FoldChange < -1, padj < 0.01) %>% rownames()

#ALLE genen opslaan
ipsc_gene_total<-ipsc_results %>% data.frame() %>% rownames()

#Data als object voor analyse opslaan
ipsc_testobject_downregulated<-new("GOHyperGParams",
                                 geneIds=ipsc_gene_downregulated,
                                 universeGeneIds=ipsc_gene_total,
                                 annotation="org.Hs.eg.db",
                                 ontology="BP",
                                 pvalueCutoff=1,
                                 testDirection="over"
)

#Analyse uitvoeren
ipsc_gotermanalysis_downregulated<-hyperGTest(ipsc_testobject_downregulated)

#Resultaten opslaan als summary
ipsc_goterm_sum_downregulated<-summary(ipsc_gotermanalysis_downregulated)

#Adjusted P-value toevoegen
ipsc_goterm_sum_downregulated$padj<-p.adjust(ipsc_goterm_sum_downregulated$Pvalue, method = "BH")

#Filteren genen kleiner dan 5, groter dan 500
ipsc_goterm_sum_downregulated<-ipsc_goterm_sum_downregulated %>% filter(Count>5 & Count<500)

#Selecteren van top 20 genen
ipsc_goterm_sum_downregulated_top20<-ipsc_goterm_sum_downregulated[order(ipsc_goterm_sum_downregulated$padj)[1:20],]

#Factor maken van "Term" en daarop sorteren
ipsc_goterm_sum_downregulated_top20$Term<-factor(
  ipsc_goterm_sum_downregulated_top20$Term,
  levels=ipsc_goterm_sum_downregulated_top20$Term[order(ipsc_goterm_sum_downregulated_top20$padj, decreasing = TRUE)]
)
is.factor(ipsc_goterm_sum_downregulated_top20$Term)

#Plotten
ipsc_goterm_sum_downregulated_top20 %>% ggplot(aes(x=Term, y=-log10(padj)))+
  geom_point()+
  coord_flip()+
  labs(
    title="Top 20 downregulated genes in iPSC's \n compared to fibroblasts",
    y="-log10 adjusted P-value",
    x="GO-terms"
  )+
  theme_bw()

#Les 4: Maak er een functie van ####

#De functie
DESeqEnrich<-function(DESeq2_results,regulated,LFC,p_val){
  #Controleren of input TRUE (upregulated) of FALSE (downregulated) is en gebaseerd hierop goed op de LFC filteren
  gene_regulated<-DESeq2_results %>% data.frame()
  if(regulated){
    gene_regulated<-gene_regulated %>% filter(log2FoldChange > LFC, padj < p_val)
  }
  else{
    gene_regulated<-gene_regulated %>% filter(log2FoldChange < -LFC, padj < p_val)
  }
  gene_regulated<-gene_regulated %>% rownames()
  
  #Totale hoeveelheid genen opslaan
  gene_total<-DESeq2_results %>% data.frame() %>% rownames()
  
  #Data als object voor analyse opslaan
  testobject<-new("GOHyperGParams",
                                   geneIds=gene_regulated,
                                   universeGeneIds=gene_total,
                                   annotation="org.Hs.eg.db",
                                   ontology="BP",
                                   pvalueCutoff=1,
                                   testDirection="over"
  )
  #Analyse uitvoeren
  gotermanalysis<-hyperGTest(testobject)
  
  #Summary verkrijgen
  summary(gotermanalysis)
  
}
#Controle: Deze twee zouden hetzelfde moeten zijn
summary(ipsc_gotermanalysis_upregulated) %>% head()
DESeqEnrich(ipsc_results,TRUE,1,0.01) %>% head()

#Wat de functie uiteindelijk moet doen:
#upregulated genen opslaan
ipsc_gene_upregulated<-ipsc_results %>% data.frame() %>% filter(log2FoldChange > 1, padj < 0.01) %>% rownames()
view(ipsc_gene_upregulated)
view(ipsc_gene_downregulated)
view(ipsc_gene_total)
view(ipsc_goterm_sum_downregulated)


#ALLE genen opslaan
ipsc_gene_total<-ipsc_results %>% data.frame() %>% rownames()

#Data als object voor analyse opslaan
ipsc_testobject_upregulated<-new("GOHyperGParams",
                                 geneIds=ipsc_gene_upregulated,
                                 universeGeneIds=ipsc_gene_total,
                                 annotation="org.Hs.eg.db",
                                 ontology="BP",
                                 pvalueCutoff=1,
                                 testDirection="over"
)

#Analyse uitvoeren
ipsc_gotermanalysis_upregulated<-hyperGTest(ipsc_testobject_upregulated)

#Resultaten opslaan als summary
ipsc_goterm_sum_upregulated<-summary(ipsc_gotermanalysis_upregulated)
