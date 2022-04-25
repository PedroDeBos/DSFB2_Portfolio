library(tidyverse)
#Opdracht 1.1 ####

#a) HASM-cellen werden verkregen van 4 witte patiënten:
  #1: onder geen behandeling
  #2: behandeling met beta-2 agonist
  #3: behandeling met glucocortisosteroid
  #4: behandeling van zowel beta-2 als gluco

#b) Voor het experiment werd Dexamethasone als glucocortisosteroid gebruikt, voor 18 uur

#c) Voor sequencing werd de illumina truSeq gebruikt, een vorm van short-read sequencing

#d) Voor het experiment is paired-end sequencing gebruikt

#Opdracht 1.2-4 ####
#Zie het bash script "importing_fastq.sh" voor de code om alle fastq bestanden te downloaden

#Opdracht 1.5 ####
##Data lijkt inderdaad in het algemeen goed, maar op sommige plekken is het wel rood, vooral in "gene duplication".
#is dat niet een significant probleem?

#Instaleren van Rsubread
install.packages("BiocManager")
BiocManager::install("Rsubread")
browseVignettes("Rsubread")
#BrowseVignettes doet het niet, sadge. Gewoon de link gebruikt


#Opdracht 1.6: ####
##Gewoon op internet opgezocht

#Code van alignment uitvoeren begrijpen
list.files("/home/daur2/rnaseq/rnaseq_airway/fastq/", pattern="_[12].fastq.gz") %>%
  str_remove(pattern="_[12].fastq.gz") %>% unique

list.files(fastq_dir, pattern = "_[12].fastq.gz") %>% 
  str_remove(pattern = "_[12].fastq.gz") %>%
  unique()

#Opdracht 1.8 ####
#Dataset tidy maken
alignment_statistics_tidy<-t(alignment_statistics) %>% as_tibble() %>% 
  mutate(
  monster=list.files("/home/daur2/rnaseq/rnaseq_airway/fastq/", pattern="_[12].fastq.gz") %>%
    str_remove(pattern="_[12].fastq.gz") %>% 
    unique
  ) %>%
  select(monster, Total_fragments:Inversed_mapping)
alignment_statistics_tidy
#Summary maken voor bar-grafiek
total<-alignment_statistics_tidy$Total_fragments %>% sum()

alignment_sum<-alignment_statistics_tidy %>% group_by(monster) %>%
  summarise(
    relative_fragments=Total_fragments/total*100
  )

alignment_sum %>% ggplot(aes(x=monster, y=relative_fragments, fill=monster))+
  geom_col(colour="black")+
  labs(
    title="Relatieve hoeveelheid DNA fragmenten tussen experimenten",
    subtitle="Ten opzichte van totaal hoeveelheid fragmenten",
    x="Monster",
    y="Percentage fragmenten (%)"
  )+
  theme_classic()+
  theme(axis.text.x=element_blank())

#FUCK ZE WOUDEN DE HOEVEELHEID UNIQUELY MAPPED KRIJGEN
alignment_sum_uniq<-alignment_statistics_tidy %>% group_by(monster) %>%
  summarise(
    relative_fragments_uniq=Uniquely_mapped_fragments/Total_fragments*100
  )

alignment_sum_uniq %>% ggplot(aes(x=monster, y=relative_fragments_uniq, fill=monster))+
  geom_col(colour="black")+
  labs(
    title="Relatieve hoeveelheid unieke DNA fragmenten",
    subtitle="Ten opzichte van totaal hoeveelheid fragmenten",
    x="Monster",
    y="Percentage fragmenten (%)"
  )+
  theme_classic()+
  theme(axis.text.x=element_blank())

