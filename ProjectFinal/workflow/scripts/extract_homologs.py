import os
import pandas as pd
from Bio import SeqIO

alignment_df = pd.read_csv(
    snakemake.input.tsv, 
    sep='\t', 
    names=['qseqid', 'sseqid', 'pident']
)

prot_dict = SeqIO.to_dict(SeqIO.parse(snakemake.input.protein_pool, "fasta"))
threshold = snakemake.params.threshold
split_dir = snakemake.output.splitted_dir

os.makedirs(split_dir, exist_ok=True)

filter_homologs = alignment_df[alignment_df['pident'] >= threshold]

group_df = filter_homologs.groupby('qseqid')

for query, group in group_df:
    output_file = os.path.join(split_dir, f"{str(query).split('|')[1]}.fasta")
    with open(output_file, "w") as f:
        for homolog in group['sseqid']:
            SeqIO.write(prot_dict[homolog], f, "fasta")