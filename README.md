# ZK_CICDUX4_POLE

The code in this repository accompanies the manuscript "The CIC::DUX4 oncoprotein maintains DNA integrity through direct regulation of the catalytic subunit of DNA polymerase epsilon (POLE)." by Kosibaty, Luck, and Okimoto (2024). These files are not plug-and-play but serve to provide transparency for how we analyzed the data. They may need to be modified for use in different environments or with different data. With questions please contact the corresponding author at ross.okimoto@ucsf.edu. RNA-seq raw and processed (counts, log2cpm) data files are available at GEO, series [GSE279467](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279467).

There are two files, with their purposes described below:

***zeinab_documentation_siPOLE_github.txt*** A .txt file containing a description of how .fastq files were processed with STAR and other tools both locally and on Wynton, the UCSF HPC.

***X1_Zeinab_siPOLE_processing.R*** An R script describing how the data was analyzed, including differential gene expression analysis with edgeR.
