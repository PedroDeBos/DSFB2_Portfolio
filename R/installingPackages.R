library(here)

#Installing DRC for RMD 1, to properly determine LD50's 
install.packages("drc")

#Instaling FS for RMD 2, creating directory trees
install.packages("fs")

#Installing vitae for RMD 3.1, creating a CV
install.packages("vitae")

#Installing pagedown for RMD 3.1, creating a cv
install.packages("pagedown")

#Installing rbbt for citation
devtools::install_github("paleolimbot/rbbt")

#Installing available too check for package names
install.packages("available")