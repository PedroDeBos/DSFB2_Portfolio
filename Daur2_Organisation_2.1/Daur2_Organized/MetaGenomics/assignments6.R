library(tidyverse)
library(dplyr)
library(png)
library(grid)
library(gridExtra)
library(reshape2)

basequality_img_1<-readPNG("~/temp_daur2/lesson6/fastqc_waternet/HU1_MOCK1_L001_R1_001_fastqc/Images/per_base_quality.png") %>%
  as.raster() %>% rasterGrob()
basequality_img_2<-readPNG("~/temp_daur2/lesson6/fastqc_waternet/HU1_MOCK1_L001_R2_001_fastqc/Images/per_base_quality.png") %>%
  as.raster() %>% rasterGrob()
grid.arrange(basequality_img_1, basequality_img_2, ncol=2, top=textGrob("Per base quality distribution of forward (right) and reverse (left) reads", gp=gpar(fontsize=10,font=8)))

#Lesson 8.1: Editing the biom data so we can properly use it #####
#Data inlezen
data<-"~/temp_daur2/lesson6/mock1/mock1_bracken_species.biom"
merged_metagenomes<-import_biom(data)
#Data inspecteren
##View(merged_metagenomes@tax_table@.Data)
#Eerste 3 (onnodige) characters verwijderen
merged_metagenomes@tax_table@.Data<-substring(merged_metagenomes@tax_table@.Data, 4)
#kolomn-namen duidelijker maken
colnames(merged_metagenomes@tax_table@.Data) <- 
  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#Soorteren op dataset van alleen bacteriën
merged_metagenomes_b<-subset_taxa(merged_metagenomes, Kingdom=="Bacteria")
sample_names(merged_metagenomes_b) <- "Bacteria"

#Opdracht 8.2: Hoeveelheid andere taxa #####
merged_metagenomes_a<-subset_taxa(merged_metagenomes, Kingdom=="Archaea") #11 Archaea
sample_names(merged_metagenomes_a) <- "Archaea"
merged_metagenomes_v<-subset_taxa(merged_metagenomes, Kingdom=="Viruses") #10 virussen
sample_names(merged_metagenomes_v) <- "Viruses"
merged_metagenomes_e<-subset_taxa(merged_metagenomes, Kingdom=="Eukaryota") #1 Eukaryoot
sample_names(merged_metagenomes_e) <- "Eukaryota"
#De reads reflecteren  dit doordat het voornamelijk bacteriën zijn: 800+ bacteriën vs 22 niet-bacteriën totaal
#Daaruit kan je duidelijk zien dat het experiment zich focust op bacteriën

#Deze merged metagenomes omzetten tot logischere data frames
data_b <- data.frame(Samples = sample_names(merged_metagenomes_b),
                     Reads = sample_sums(merged_metagenomes_b))

data_e <- data.frame(Samples = sample_names(merged_metagenomes_e),
                     Reads = sample_sums(merged_metagenomes_e))

data_a <- data.frame(Samples = sample_names(merged_metagenomes_a),
                     Reads = sample_sums(merged_metagenomes_a))

data_v <- data.frame(Samples = sample_names(merged_metagenomes_v),
                     Reads = sample_sums(merged_metagenomes_v))

data_t<-rbind(data_b,data_a,data_v,data_e)
#Maken van een grafiek gebaseerd op de data om goed te zien wat de schaal in groteverschil is

data_t %>% ggplot(aes(x=Samples, y=Reads, fill=Samples))+
  geom_col()+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(
    title="Hoeveelheid reads van taxa in waternet sample",
  )

#Lesson 8.4: Visualiseren van geidentificeerde soorten ####
#Zelfde stappen als hierboven uitvoeren: importeren, onnodige characters verwijderen, colnames toevoegen, otu-table mock1 noemen
merged_metagenomes<-import_biom(data)
merged_metagenomes@tax_table@.Data<-substring(merged_metagenomes@tax_table@.Data,4)
#view(merged_metagenomes@tax_table@.Data)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(merged_metagenomes@otu_table) <- c("mock1")

#checken voor lege labels
summary(merged_metagenomes@tax_table@.Data=="")

#Data klaarzetten voor transformatie tot plot
glom<-tax_glom(merged_metagenomes, taxrank = "Species")
mock1_metagenome_species<-psmelt(glom)

#Combining species and genus for plot
mock1_metagenome_species$Species <- as.character(mock1_metagenome_species$Species)
mock1_metagenome_species$Species<-paste(mock1_metagenome_species$Genus, mock1_metagenome_species$Species)

#Making the plot
mock1_metagenome_species %>% ggplot(aes(x=Sample, y=Abundance, fill=Species))+
  geom_col()

#Selecting for abundant species
mock1_metagenome_species$Species[mock1_metagenome_species$Abundance < 160000] <- "Species < 160.000 abund."
#Making the plot again
mock1_metagenome_species %>% ggplot(aes(x=Sample, y=Abundance, fill=Species))+
  geom_col(aes(), position="stack")      

#Glomming again, this time for percentage graphs
glom<-tax_glom(merged_metagenomes, taxrank="Species")
mock1_metagenome_species_percentage<-psmelt(glom)
mock1_metagenome_species_percentage$Abundance<-mock1_metagenome_species_percentage$Abundance/sum(mock1_metagenome_species_percentage$Abundance)*100
#Selecting all species with a abundance higher than 0.5
mock1_metagenome_species_percentage$Species<-paste(mock1_metagenome_species_percentage$Genus, mock1_metagenome_species_percentage$Species)
mock1_metagenome_species_percentage$Species[mock1_metagenome_species_percentage$Abundance<0.5] <- "Species < 0.5% abund."

mock1_metagenome_species_percentage %>% ggplot(aes(x=Sample, y=Abundance, fill=Species))+  
geom_col(aes(), position="stack") 

#Lesson 8.5: Vergelijken met de werkelijke inhoud van het monster
mock1_composition<-as.data.frame(read.csv("/home/daur2/metagenomics/reader_data/HU_waternet_MOCK1_composition.csv", sep=";", row.names=1))
#$Amount omzetten tot numeric (wat het al was??? maar ach ja)
mock1_composition$amount<-gsub(",", ".", mock1_composition$amount) %>% as.numeric()
#Colnames herkenbaar maken
colnames(mock1_composition) <- c("name","amount","amountP","sample_name","total_volume")

#Species bestanden intersecten
mock1_and_composition_intersect<-mock1_metagenome_species_percentage[mock1_metagenome_species_percentage$Species %in% mock1_composition$name,]
`%!in%` <- Negate(`%in%`)
comp_not_in_mock1<-mock1_composition[mock1_composition$name %!in% mock1_metagenome_species_percentage$Species,]
unique(mock1_and_composition_intersect$Species)
unique(comp_not_in_mock1$name)

mock1_and_composition_intersect$amountP <- NA
for (m1_label in mock1_and_composition_intersect$Species){
  for (m1c_label in mock1_composition$name){
    if(m1_label == m1c_label){
      mock1_and_composition_intersect$amountP[mock1_and_composition_intersect$Species == m1_label] <- mock1_composition$amountP[mock1_composition$name == m1c_label]
    }
  }
}

mock1_and_comp_plotting_data           <- mock1_and_composition_intersect[,c(10,3,11)]
mock1_and_comp_plotting_data$amountP<-gsub(",", ".", mock1_and_comp_plotting_data$amountP) %>% as.numeric()
colnames(mock1_and_comp_plotting_data) <- c("species", "k_abundance", "c_abundance")

mock1_and_comp_plotting_data       <- melt(mock1_and_comp_plotting_data, id.var = "species")
mock1_and_comp_plotting_data$value <- as.numeric(mock1_and_comp_plotting_data$value)

ggplot(mock1_and_comp_plotting_data, aes(x = species, y = value, fill = variable)) + 
  geom_bar(aes(), stat="identity", position="dodge") +
  theme_classic() +
  ylab("Abundance (%)") +
  xlab("") +
  ggtitle("Abundance comparison between Kraken2 results and composition") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0,25) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1)) +
  scale_fill_manual(values=c("skyblue", "orangered"))
