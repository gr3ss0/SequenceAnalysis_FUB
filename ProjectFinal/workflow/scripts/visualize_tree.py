from Bio import Phylo
import matplotlib.pyplot as plt

# Fix: Use input[0] and output[0] natively instead of '{input}' strings
tree = Phylo.read(snakemake.input[0], 'newick')
tree.ladderize()

num_leaves = tree.count_terminals()
logging.info(f"Tree contains {num_leaves} leaf nodes.")

# Adjust height proportionally so labels don't overlap (minimum 6 inches)
fig_height = max(6, num_leaves * 0.25)
fig_width = 10

logging.info(f"Setting figure size to: ({fig_width}, {fig_height})")
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

Phylo.draw(tree, axes=ax, do_show=False)

plt.savefig(snakemake.output[0], dpi=300, bbox_inches='tight')
plt.close()