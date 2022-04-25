#!/bin/bash
#Maken van een conda in temp_daur2
##conda create -n temp_daur2
##conda activate temp_daur2
##fastq-dump --help

#Split-3 lijkt het bestand in 3 files te splitsen? Splitst in _1.fastq en _2.fastq als er 2 reads zijn
#, blijft als origineel .fastq als er maar 1 read is

#--outdir zorgt ervoor dat het bestand in een andere locatie neerzet dan in de directory waar de originele
#Fastq bestanden in staan

#--gzip zorgt ervoor dat de output ge-compressed is

#Bepalen hoe het command werkt
##fastq-dump --split-3 --outdir '/home/1762403/temp_daur2/lesson1.' --gzip SRR1039508

#For-do-done loop opstellen
##for SRR in $(cat srr_identifier.txt)
##do
  ##echo "fastq-dump --split-3 --outdir '/home/1762403/temp_daur2/lesson1.' --gzip $SRR"
##done

#zonder terminal dat cat bestandje maken
echo "SRR7866687 SRR7866688 SRR7866689 SRR7866690 SRR7866691 SRR7866692 SRR7866693 SRR7866694" > SRR_identifier.txt

for SRR in (cat SRR_identifier.txt)
do
  echo $SRR
done

