library(here)

#Installing DRC for RMD 1, to properly determine LD50's 
install.packages("drc")

#Instaling FS for RMD 2, creating directory trees
install.packages("fs")

fs::dir_tree(path = "Daur2_Organisation_2.1/Daur2_Original/")

fs::dir_tree(path = here("Daur2_Organisation_2.1/Daur2_Organized/"))
