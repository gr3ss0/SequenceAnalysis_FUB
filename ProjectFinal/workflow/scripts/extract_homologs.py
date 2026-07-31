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
output_files = {
    os.path.basename(f).replace(".fasta", ""): f
    for f in snakemake.output.fasta_files
}

alignment_df["sample"] = (
    alignment_df["sseqid"]
    .str.rsplit("_", n=1)
    .str[0]
)

alignment_df = alignment_df.sort_values(
    by=["qseqid", "sample", "pident"],
    ascending=[True, True, False]
)


top_hits = (
    alignment_df
    .groupby(["qseqid", "sample"])
    .head(2)
) #only taking the top 2 hits for each query and sample

for query, group in top_hits.groupby('qseqid'):
    output_file = output_files[str(query)]
    with open(output_file, "w") as f:
        for protein_id in group['sseqid']:
            SeqIO.write(prot_dict[protein_id], f, "fasta")









#filter_homologs = alignment_df[alignment_df['pident'] >= threshold].copy()



#group_df = filter_homologs.groupby('qseqid')

#for query, group in group_df:
    #output_file = output_files[str(query)]
    #with open(output_file, "w") as f:
        #for homolog in group['sseqid']:
            #SeqIO.write(prot_dict[homolog], f, "fasta")