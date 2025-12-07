# BIOLM0051-Introduction-to-Bioinformatics-Short-Analysis-Pipeline

Short Analysis Pipeline for analysing unknown samples by analysing their FASTQ files and constructing a phylogenetic tree for deducing the organism in each sample

## Steps involved in Analysis of unkown Samples

1. Navigate to the "01-Sample Data" sub-folder and you will find the sample FASTQ Files. Alongwith the FASTQ files there are 3 Bash Scripts, which are numbered. Run them accordingly. These scripts will clean the FASTQ files, convert them to FASTA files and combine the FASTA files by sample into another directory name "combine_files"

2. Now using these FASTA files, go to NCBI BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome), to "blastn suite". Paste or upload the file into the query textbox. Select "Core nucleotide database (nt). Under Genreal Parameters select the sequence threshold as "500" and word size as "28". BLAST the sequence data. Do this for the 4 samples and select the top hits of phylogentically different species looking through the description tab. Download the "FASTA(alignedsequences)" file from the BLAST window. These files are in "02-BLAST Data"

3. Now take these BLAST fasta files and paste them in the "03-Translated Data" and also paste the converted sample FASTA files in the same directory. There are two python scripts there named "translating.py" and "fasta_cleaning.py". The "translating.py" script changes the BLAST FASTA files from DNA to amino acid sequences and saves as "translated.fasta". The "fasta_cleaning.py" script cleans the "translated.fasta" files and gives a "cleantranslated.fasta".

4. Now take "cleantranslated.fasta" and go to https://www.ebi.ac.uk/jdispatcher/ for running Multiple Sequence Alignment using MUSCLE which uses ClustalW in a more precise manner. This will give an alignment file in ".fa" format and a phylogenetic tree which can be downloaded and viewed in the web viewer. Save these files in "04-Alignment Data". 

5. Now using this alignment data we can construct a phylogenetic tree. Inside the "05-Phylogenetic Tree Data", theres a script named "tree.py". This generates a tree with all the samples. This generates a tree with all the samples. Take the files from the "04-Alignemt Data" folder and paste it in this directory to get results.

## NOTE: 

The folder "00-FASTQC Quality Analysis" contains the Quality Analysis reports of the FASTQ files given done for preliminary analysis. This is purely for references. It also contains a ".slurm" script to run the FASTQ analysis on an HPC. 
