import re


def clean_fasta_header(header):
    """
    Extract accession number and scientific name from FASTA header.
    Handles common formats like:
    >ACC123456.1 Homo sapiens [additional info]
    >gi|123456|ref|ACC123.1| Genus species [organism info]
    """
    # Remove the leading '>'
    header = header.lstrip('>')

    # Split by whitespace to get parts
    parts = header.split()

    if len(parts) < 2:
        # If not enough parts, return as is with underscores
        return '>' + '_'.join(parts)

    # First part is typically the accession
    # Handle piped formats like gi|123|ref|ACC
    accession = parts[0].split('|')[-1]

    # Scientific name is typically the next two words (Genus species)
    if len(parts) >= 3:
        scientific_name = f"{parts[1]}_{parts[2]}"
    else:
        scientific_name = parts[1]

    # Combine accession and scientific name
    return f">{accession}_{scientific_name}"


def process_fasta(input_file, output_file):
    """Process FASTA file and clean headers"""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Clean the header
                cleaned = clean_fasta_header(line.strip())
                outfile.write(cleaned + '\n')
            else:
                # Write sequence lines as is
                outfile.write(line)

    print(f"Successfully processed {input_file}")
    print(f"Output written to {output_file}")


# Example usage:
if __name__ == "__main__":
    # Simply change these file paths to your input and output files
    input_fasta = "translated.fasta"
    output_fasta = "cleantranslated.fasta"

    process_fasta(input_fasta, output_fasta)
