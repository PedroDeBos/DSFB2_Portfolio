library(here)

#Installing DRC for RMD 1, to properly determine LD50's 
install.packages("drc")

#Instaling FS for RMD 2, creating directory trees
install.packages("fs")

#Installing vitae for RMD 3.1, creating a CV
install.packages("vitae")

#Installing pagedown for RMD 3.1, creating a cv
install.packages("pagedown")

#Installing devtools to install rbbt
install.packages("devtools")

#Installing rbbt for citation
devtools::install_github("paleolimbot/rbbt")

#Installing available too check for package names
install.packages("available")

#Installing class for machine learning
install.packages("class")

install.packages("dslabs")

install.packages("RPostgres")
install.packages("RPostgreSQL")

#Installing pokemonAnalyse to analyse data
install_github("PedroDeBos/PokemonAnalyse", build_vignettes = TRUE)
