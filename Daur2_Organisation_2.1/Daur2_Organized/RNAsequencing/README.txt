In this directory, all R code, bash code and local data files for performing an RNA sequencing experiment have been collected.

"Bash" contains all code that's to be performed in bash, specifically for downloading all necessary fastq files

"R" contains all code for performing the RNAsequencing analysis. They are named after assignments described in the online course "daur2-reader". 

Assignments1 describes checking FASTQC-quality files and performing an RNA allignment
Assignments2 describes generating a count table, creating a DESEQ2 object, normalising the data, performing a PCA analysis and creating a heatmap
Assignments3 describes performing a DGE-analysis and visualising the DGE-analysis via count plots, volcano plots and heatmaps
Assignment4 describes how to couple entrez gene-identifiers to gene symbols, obtaining GO-terms and performing a GO-enrichtment analysis.

"Data" contains all local data on which the R code depends to properly copy the RNA-sequencing data from the server.
  srr_identifier.txt is used to store the SRR-identiifers which's data needs to be analysed.
  test_srr_identifiers.txt was used to test out in which way the SRR-identifiers had to be structured: in the same row or in a column underneath eachother
  it is unkown for what "randomdata.txt" was used.