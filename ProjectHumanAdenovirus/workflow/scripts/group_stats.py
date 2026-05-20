import pandas as pd
import os
import re
'''
idxstats returns
reference sequence name, sequence length, # mapped reads and # unmapped 
'''
gamma = pd.DataFrame(columns=['ref_name','seq_len', 'mapped', 'unmapped'])
for stat in snakemake.input.list_of_stats:
    new = pd.read_csv(stat, sep='\t', names=['ref_name','seq_len', 'mapped', 'unmapped'])
    gamma = pd.concat([gamma, new])

out = gamma.groupby('ref_name').agg({'seq_len': 'max','mapped':'sum'})
out.to_csv(snakemake.output[0], sep='\t', header=False)