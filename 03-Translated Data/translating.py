from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob
import warnings
from Bio import BiopythonWarning
warnings.filterwarnings("ignore", category=BiopythonWarning)
import os


# Load sequences from multiple FASTA files
sequences = [record
             for fasta in glob.glob("*.fasta")
             for record in SeqIO.parse(fasta, "fasta")]

print(sequences)


aa_sequences = []   # List to store output amino-acid SeqRecord objects

for record in sequences:
    dna_seq = record.seq   # Nucleotide sequence from the FASTA record

    best_orf = ""          # Will store the longest ORF found across all frames
    best_frame = None      # Will store which frame (0,1,2) contains the longest ORF

    # Translate sequence in all three forward reading frames
    for frame in range(3):
        # Translate starting at this frame; keep stop codons (to_stop=False)
        protein = dna_seq[frame:].translate(to_stop=False)

        # Split translation into fragments separated by stop codons ("*")
        # Each fragment is a potential ORF (no internal stops)
        fragments = str(protein).split("*")

        # Find the longest continuous amino-acid stretch (longest ORF) in this frame
        longest_in_frame = max(fragments, key=len)

        # If this ORF is longer than the best one we've seen so far, store it
        if len(longest_in_frame) > len(best_orf):
            best_orf = longest_in_frame
            best_frame = frame

    # Convert the longest ORF string back into a Seq object
    best_protein = Seq(best_orf)

    # Extract species name (Latin binomial) from the FASTA description
    # Assumes format: >ID Genus species
    start = record.description.find("[organism=") + len("[organism=")
    end = record.description.find("]", start)
    latin_name = record.description[start:end].replace(" ", "_")  # replace space with underscore


    # Create a new SeqRecord containing the longest ORF for this sequence
    aa_record = SeqRecord(
        best_protein,
        id=record.id,                      # Keep the same FASTA ID
        name=latin_name,                   # Store species name
        description=record.description,    # Preserve full original header
        annotations={
            "type": "longest ORF",
            "frame": best_frame            # Which reading frame the ORF came from
        }
    )

    aa_sequences.append(aa_record)         # Add to output list

# Write all translated longest-ORF sequences to a FASTA file
SeqIO.write(aa_sequences, "translated.fasta", "fasta")

