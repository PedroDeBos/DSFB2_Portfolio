#Oh mijn fucking god het komt omdat ik één keer een onnodige %>% stap had toegevoed ik wil dood
#Selecteren van upregulated, downregulated and all genes
neuron_dge_upregulated<-neuron_dge_results %>% data.frame() %>% filter(log2FoldChange>1 & padj<0.01) %>% rownames()
neuron_dge_downregulated<-neuron_dge_results %>% data.frame() %>% filter(log2FoldChange< -1 & padj<0.01) %>% rownames()
neuron_dge_all<-neuron_dge_results %>% data.frame() %>% rownames()


#UPREGULATED FOUT ####
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
hyperG_top20<-neuron_dge_hyperG_result %>% head(20)

#Creating a factor so that the graph goes the proper way
hyperG_top20$Term<-hyperG_top20$Term %>% factor(hyperG_top20$Term, levels=(hyperG_top20$Term[order(hyperG_top20$padj, decreasing=TRUE)]))

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

# DOWNREGULATED FOUT ####
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
hyperG_top20_down$Term<-hyperG_top20$Term %>% factor(hyperG_top20$Term, levels=(hyperG_top20$Term[order(hyperG_top20$padj, decreasing=TRUE)]))

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

