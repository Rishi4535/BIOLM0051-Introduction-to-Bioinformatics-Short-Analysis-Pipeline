from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob
import warnings
from Bio import BiopythonWarning
warnings.filterwarnings("ignore", category=BiopythonWarning)
import os
from Bio import Phylo

alignment = AlignIO.read("aa_aligned.fa", "fasta")

# Replace '?' with 'X' in all sequences
cleaned_records = []
for record in alignment:
    cleaned_seq = Seq(str(record.seq).replace("?", "X"))
    cleaned_record = SeqRecord(
        cleaned_seq,
        id=record.id,
        name=record.name,
        description=record.description
    )
    cleaned_records.append(cleaned_record)

# Create a new MultipleSeqAlignment object
cleaned_alignment = MultipleSeqAlignment(cleaned_records)


calculator = DistanceCalculator("blosum62")
dm = calculator.get_distance(cleaned_alignment)

print("Distance matrix:")
print(dm)




constructor = DistanceTreeConstructor()
nj_tree = constructor.nj(dm)
Phylo.draw_ascii(nj_tree)

import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
Phylo.draw(nj_tree, do_show=True)

print([term.name for term in nj_tree.get_terminals()])
