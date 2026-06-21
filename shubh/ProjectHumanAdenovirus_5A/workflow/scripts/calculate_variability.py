from Bio import AlignIO
import numpy as np
import sys

# Redirect standard output and error to the log file defined in the Snakemake rule
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

MSA = AlignIO.read(snakemake.input.msa, "fasta")


def shannon_entropy(column):
    from collections import Counter
    import math
    
    counts = Counter(column)
    total = sum(counts.values())
    
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    
    return entropy



entropies = []
for i in range(MSA.get_alignment_length()):
    column = MSA[:, i]
    variability = shannon_entropy(column)
    entropies.append(variability)
    print(f"Position {i+1}: Variability = {variability}")

# Save the variability values to a text file
np.savetxt(snakemake.output.variability, entropies, fmt="%.6f", delimiter="\t")