import os
import pandas as pd
from Bio import SeqIO

# -------------------------------------------------------------------
# 1. Load Alignment TSVs
# -------------------------------------------------------------------
dfs_list = []
for tsv in snakemake.input.tsvs:
    df = pd.read_csv(tsv, sep='\t', names=['qseqid', 'sseqid', 'pident'])
    dfs_list.append(df)

if dfs_list:
    df_all = pd.concat(dfs_list, ignore_index=True)
else:
    df_all = pd.DataFrame(columns=['qseqid', 'sseqid', 'pident'])

# -------------------------------------------------------------------
# 2. Index Protein Pools (Memory-Efficient)
# -------------------------------------------------------------------
protein_dict = {}

# SeqIO.to_dict loads full objects into RAM; for large fasta files,
# updating a standard dictionary with Parsed records is fine, but
# doing it in a single dictionary comprehension or keeping indices avoids overhead.
for protein_pool in snakemake.input.protein_pools:
    protein_dict.update(SeqIO.to_dict(SeqIO.parse(protein_pool, "fasta")))

# -------------------------------------------------------------------
# 3. Filter and Split by Query
# -------------------------------------------------------------------
threshold = float(snakemake.params.threshold)
split_dir = snakemake.output.splitted_dir

os.makedirs(split_dir, exist_ok=True)

# Filter hits by sequence identity threshold
filter_homologs = df_all[df_all['pident'] >= threshold]

# Group hits by query protein ID
grouped = filter_homologs.groupby('qseqid')

for query, group in grouped:
    query_str = str(query)
    
    # Safely extract query name if formatted like "sp|P12345|NAME"
    clean_query_id = query_str.split('|')[1] if '|' in query_str else query_str
    
    output_file = os.path.join(split_dir, f"{clean_query_id}.fasta")
    
    # Deduplicate subject IDs to avoid writing duplicate sequences
    unique_homologs = set(group['sseqid'])
    
    # Collect records and write out
    records_to_write = [
        protein_dict[homolog] 
        for homolog in unique_homologs 
        if homolog in protein_dict
    ]
    
    if records_to_write:
        SeqIO.write(records_to_write, output_file, "fasta")

# import os
# import pandas as pd
# from Bio import SeqIO

# dfs_list = []
# for tsv in snakemake.input.tsvs:
#     df = pd.read_csv(tsv, sep='\t',names=['qseqid', 'sseqid', 'pident'])
#     dfs_list.append(df)

# df_all = pd.concat(dfs_list)

# protein_dict = {}

# # Iterate over each file path in snakemake.input.protein_pools
# for protein_pool in snakemake.input.protein_pools:
#     # Update the master dictionary with records from the current fasta file
#     protein_dict.update(SeqIO.to_dict(SeqIO.parse(protein_pool, "fasta")))
    
# threshold = snakemake.params.threshold
# split_dir = snakemake.output.splitted_dir

# os.makedirs(split_dir, exist_ok=True)

# filter_homologs = df_all[df_all['pident'] >= threshold]

# group_df = filter_homologs.groupby('qseqid')

# for query, group in group_df:
#     output_file = os.path.join(split_dir, f"{str(query).split('|')[1]}.fasta")
#     with open(output_file, "w") as f:
#         for homolog in group['sseqid']:
#             SeqIO.write(protein_dict[homolog], f, "fasta")