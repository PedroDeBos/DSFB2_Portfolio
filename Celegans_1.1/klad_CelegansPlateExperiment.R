#Basic functions
library(tidyverse)
#Library to create dose-response curves
library(drc)
#Library to read excel documents
library(readxl)

#Reading the dataset ####

#storing the location of the dataset
excel_location<-"/Users/pedro/Downloads/CE.LIQ.FLOW.062_Tidydata.xlsx"

#Reading the excel file using read_excel()
Celegans_data<-read_excel(excel_location, sheet = 1)

#Dataset inspection ####

#Inspecting the full data set
Celegans_data %>% view()
#We want to study de effect of different components in different concentrations on the amount of offspring created by the
#Elegans after exposure. 

#Selecting data for this goal
Celegans_data_select<-Celegans_data %>% select(RawData, compName, compConcentration)

#Comp concentration is stored as a character, transforming into a numeric variable
Celegans_data_select$compConcentration %>% as.numeric()

#Leads to creation of an NA, inspectin data to figure out why
Celegans_data_select$compConcentration[259]

#Datapoint 259 has a comma instead of a point. Transforming value via str_replace
Celegans_data_select$compConcentration<-Celegans_data_select$compConcentration %>% str_replace(",", ".")

#checking if datapoint still has a comma
Celegans_data_select$compConcentration[259]

#Now properly transforming the compConcentration data into numeric
Celegans_data_select$compConcentration<-Celegans_data_select$compConcentration %>% as.numeric()

#Checking of transformation went properly
Celegans_data_select

#Transforming compName into a factor variable
levels_Celegans<-unique(Celegans_data_select$compName) #Storing all levels
parse_factor(Celegans_data_select$compName, levels = levels_Celegans) #Checking of the factor gives any errors
Celegans_data_select$compName<-factor(Celegans_data_select$compName, levels = levels_Celegans) #Transforming compName into factor

#Checking if transformation went properly
Celegans_data_select
