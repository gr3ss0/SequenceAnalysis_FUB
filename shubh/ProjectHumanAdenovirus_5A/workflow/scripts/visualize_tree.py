from Bio import Phylo
import matplotlib.pyplot as plt

# Fix: Use input[0] and output[0] natively instead of '{input}' strings
tree = Phylo.read(snakemake.input[0], 'newick')
tree.ladderize()

fig, ax = plt.subplots(figsize=(10, 8))
Phylo.draw(tree, axes=ax, do_show=False)

plt.savefig(snakemake.output[0], dpi=300, bbox_inches='tight')
plt.close()