import os
import pandas as pd
from Bio import SeqIO

dfs_list = []
for tsv in snakemake.input.tsvs:
    df = pd.read_csv(tsv, sep='\t',names=['qseqid', 'sseqid', 'pident'])
    dfs_list.append(df)

df_all = pd.concat(dfs_list)

protein_dict = {}

# Iterate over each file path in snakemake.input.protein_pools
for protein_pool in snakemake.input.protein_pools:
    # Update the master dictionary with records from the current fasta file
    protein_dict.update(SeqIO.to_dict(SeqIO.parse(protein_pool, "fasta")))
    
threshold = snakemake.params.threshold
split_dir = snakemake.output.splitted_dir

os.makedirs(split_dir, exist_ok=True)

filter_homologs = df_all[df_all['pident'] >= threshold]

group_df = filter_homologs.groupby('qseqid')

for query, group in group_df:
    output_file = os.path.join(split_dir, f"{str(query).split('|')[1]}.fasta")
    with open(output_file, "w") as f:
        for homolog in group['sseqid']:
            SeqIO.write(protein_dict[homolog], f, "fasta")