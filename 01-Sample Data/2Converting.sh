#!/bin/bash
#Last line in SampleA_part1 has an unnecessary symbol
sed -i '$d' sampleA_part1.fastq
# Enable nullglob so patterns with no matches are skipped
shopt -s nullglob



# Collect all FASTQ files (case-insensitive)
files=( *.fastq *.FASTQ )

# Extract unique sample prefixes
samples=()
for f in "${files[@]}"; do
    # Remove _part_X.fastq or _part_X.FASTQ suffix to get sample name
    sample=$(echo "$f" | sed -E 's/(_part_[0-9]+)?\.(fastq|FASTQ)//')
    samples+=("$sample")
done

# Get unique sample names
unique_samples=($(echo "${samples[@]}" | tr ' ' '\n' | sort -u))

# Loop over each sample
for sample in "${unique_samples[@]}"; do
    combined_file="${sample}.fasta"
    > "$combined_file"  # empty/create output file

    # Find all matching files for this sample
    for fq in "${sample}"*.fastq "${sample}"*.FASTQ; do
        [ -f "$fq" ] || continue  # skip if file doesn't exist

        # Convert FASTQ/FASTQ to FASTA and clean header
        awk -v sample="$sample" 'NR%4==1 {print ">" sample} NR%4==2 {print}' "$fq" >> "$combined_file"
    done
done
