#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=0:5:00
#SBATCH --account=BISC036945

# Load the necessary module

module load fastqc/0.12.1

# Run FastQC on all the data files with the help of the wildcard "*"

fastqc Data_practical/mystery_data_*.fastq.gz
