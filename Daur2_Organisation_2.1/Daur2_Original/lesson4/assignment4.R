library(tidyverse)
library(ggplot2)
library(DESeq2)
library(Rsubread)
library(pheatmap)
library(here)
library("org.Hs.eg.db")
library(GOstats)
#Gen-namen krijgen ####
top10_genes <- airway_dge_results[order(airway_dge_results$padj)[1:10],] %>% data.frame()
BiocManager::install("org.Hs.eg.db")
columns(org.Hs.eg.db)
help("symbol")
top10_genes<-top10_genes %>% mutate(entrezid=rownames(top10_genes))
top10_genes$symbol<-mapIds(org.Hs.eg.db,
       keys=top10_genes$entrezid,
       column = "SYMBOL",
       keytype = "ENTREZID",
       multiVals = "first")

#Opdracht 1: Wat betekent mapIds()? ####
?mapIds

#Opdracht 3 ####
top10_genes$chromosonal_location<-mapIds(org.Hs.eg.db,
       keys=top10_genes$entrezid,
       column = "MAP",
       keytype = "ENTREZID",
       multiVals = "first")

# GO-term vinden ####
top_gene<-top10_genes[which.max(top10_genes$log2FoldChange),"entrezid"]
top_gene_go<-select(org.Hs.eg.db,
       keys=top_gene,
       column=c("GO","ONTOLOGY"),
       keytype = "ENTREZID",
       multiVals = "list"
       ) %>% filter(ONTOLOGY=="BP")
unique(top_gene_go$GO)

library(GO.db)
GOterms_description_airway<-select(GO.db, keys=unique(top_gene_go$GO), columns="DEFINITION",keytype="GOID")
GOterms_description_airway$DEFINITION

# Opdracht 4.5 ####
# Vinden van de GO:term gerelateerd aan het gen
gene_1843_go<-select(org.Hs.eg.db,
       keys="1843",
       columns = c("GO","ONTOLOGY"),
       keytype="ENTREZID",
       multiVals="list"
       ) %>%
  filter(ONTOLOGY=="BP")
gene_1843_description<-select(GO.db, keys=unique(gene_1843_go$GO), columns="DEFINITION",keytype="GOID")
gene_1843_description_noNA<-gene_1843_description[!is.na(gene_1843_description$DEFINITION),]
gene_1843_description_noNA$DEFINITION[str_detect(gene_1843_description_noNA$DEFINITION, pattern = "gluco")]

# GO-term enrichtment assay ####
#Vinden van significante genen
upregulated_genes_airway<-airway_dge_results %>% data.frame() %>% filter(log2FoldChange>1 & padj < 0.01) %>% rownames()
#Vinden van ALLE genen
all_genes_airway<-airway_dge_results %>% data.frame() %>% rownames()
#analyse uitvoeren
test_object<-new("GOHyperGParams",
    geneIds=upregulated_genes_airway,
    universeGeneIds=all_genes_airway,
    annotation="org.Hs.eg.db",
    ontology="BP",
    pvalueCutoff=1,
    testDirection="over"
    )
goterm_analysis<-hyperGTest(test_object)
#Resultaten opslaan als summary
go_term_results<-summary(goterm_analysis)
#P.value adjusten, omdat de we vorige test dat niet automatisch hebben laten doen
go_term_results$padj<-p.adjust(go_term_results$Pvalue,method="BH")
#Filteren van genen kleiner dan 5 of groter dan 500
go_term_results<-go_term_results %>% filter(Count>5 & Count<500)
#Selecteren van top 20 genen
go_term_top20<-go_term_results[order(go_term_results$padj)[1:20],]
#Dit plotten, eerst een factor maken van "term"
go_term_top20$Term<-factor(go_term_top20$Term, levels=go_term_top20$Term[order(go_term_top20$padj, decreasing=TRUE)])

go_term_top20 %>% ggplot(aes(x=Term, y=-log10(padj)))+
  geom_point()+
  coord_flip()+
  labs(
    title="Top 20 genen \nin upregulated genes",
    y="-log10 adjusted P-value",
    x="GO-terms"
  )+
  theme_bw()

# Opdracht 6: Hetzelfde doen maar dan voor downregulated genes ####
#alle downregulated genes krijgen
downregulated_genes_airway<-airway_dge_results %>% data.frame() %>% filter(log2FoldChange < -1 & padj < 0.01) %>% rownames()
#Downregulated object maken
downregulated_genes_testobject<-new("GOHyperGParams",
                 geneIds=downregulated_genes_airway,
                 universeGeneIds=all_genes_airway,
                 annotation="org.Hs.eg.db",
                 ontology="BP",
                 pvalueCutoff=1,
                 testDirection="over"
)
#GO-term analyse
downregulated_goterm_analysis<-hyperGTest(downregulated_genes_testobject)
#Resultaten eruit halen
goterm_analysis_down_results<-summary(downregulated_goterm_analysis)
#padj toevoegen
goterm_analysis_down_results$padj<-p.adjust(goterm_analysis_down_results$Pvalue, method="BH")
#selecteren GO-terms kleiner dan 500, groter dan 5
goterm_analysis_down_results<-goterm_analysis_down_results %>% filter(Count>5 & Count<500)
#Selecteren van top 20
goterm_analysis_down_results<-goterm_analysis_down_results[order(goterm_analysis_down_results$padj)[1:20],]
#Tot factor veranderen
goterm_analysis_down_results$Term<-factor(goterm_analysis_down_results$Term, levels=goterm_analysis_down_results$Term[order(goterm_analysis_down_results$padj, decreasing=TRUE)])
#Plotten
goterm_analysis_down_results %>% ggplot(aes(x=Term, y=-log10(padj)))+
  geom_point()+
  coord_flip()+
  labs(
    title="Top 20 genen \nin downregulated genes",
    y="-log10 adjusted P-value",
    x="GO-terms"
  )+
  theme_bw()
#Antwoord kopieren ####
# Create a list of upregulated genes
downregulated_genes <- airway_dge_results %>% data.frame() %>% 
  filter(log2FoldChange < -1, padj < 0.01) %>% rownames()

# Create a list of all genes in the dataset
all_genes <- airway_dge_results %>% data.frame() %>% rownames()

# Perform GO term enrichment analysis
test_object2 <- new("GOHyperGParams",
                    geneIds = downregulated_genes,
                    universeGeneIds = all_genes, 
                    annotation = "org.Hs.eg.db", 
                    ontology = "BP", 
                    pvalueCutoff = 1,
                    testDirection = "over")
goterm_analysis2 <- hyperGTest(test_object2)

# Obtains dataframe with results of GO term analysis
goterm_analysis_results2 <- summary(goterm_analysis2) 

# Adjust the p values for multiple testing
goterm_analysis_results2$padj <- p.adjust(goterm_analysis_results2$Pvalue, method = "BH")

# Select only gene sets with 5 < gene count < 500
goterm_analysis_results2 <- goterm_analysis_results2 %>% filter(Count > 5) %>% filter(Count < 500)

# Select the top 20 GO terms
goterm_analysis_top20 <- goterm_analysis_results2[order(goterm_analysis_results2$padj)[1:20],]

# Plot the p-values of the top 20 GO terms
goterm_analysis_top20$Term <- factor(goterm_analysis_top20$Term, 
                                     levels = goterm_analysis_top20$Term[
                                       order(goterm_analysis_top20$padj, decreasing = TRUE)])
goterm_analysis_top20 %>% ggplot(aes(x = Term, y = -log10(padj))) +
  geom_point() +
  coord_flip() +
  ylab(expression(-log[10](adjusted~italic(P)~value))) + 
  xlab("GO terms") +
  ggtitle("Top 20 enriched GO terms\n for downregulated genes") +
  theme_bw()
