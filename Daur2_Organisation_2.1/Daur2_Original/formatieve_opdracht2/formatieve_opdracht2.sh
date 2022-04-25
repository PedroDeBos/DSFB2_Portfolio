#!/bin/bash
#Onderdeel 1: For-do-done loop voor fastqc uittesten 
#Testing for-do-done loop
for fastq in /home/daur2/metagenomics/formative_data/HU2_*
do
  echo "fastqc -o ~/temp_daur2/lesson6/fastqc_waternet $fastq"
done