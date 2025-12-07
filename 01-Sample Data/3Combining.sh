outdir="combined_fasta"
mkdir -p "$outdir"

for sample in sampleA sampleB sampleC sampleD; do
    {
        # write a single header
        echo ">${sample}"

        # concatenate all files, remove all lines starting with ">"
        cat ${sample}_part1.fasta ${sample}_part2.fasta ${sample}_part3.fasta \
            | grep -v '^>' \
            | tr -d '[:space:]'
    } > "${outdir}/${sample}.fasta"
done


