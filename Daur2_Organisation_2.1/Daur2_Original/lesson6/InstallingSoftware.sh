#!/bin/bash

#touch setup_meta_env.yml

#Setting an alias for meta
#conda activate meta

#Trying to set up a for-do-done loop
##for fastq in /home/daur2/metagenomics/formative_data/HU2_*
##do
##echo $fastq
##done

#Trying to set up a for-do-done loop #2
for fastqc in /home/1762403/temp_daur2/lesson6/fastqc_waternet/*.zip
do
unzip -d /home/1762403/temp_daur2/lesson6/fastqc_waternet $fastqc
done
#Testing out GZIP
