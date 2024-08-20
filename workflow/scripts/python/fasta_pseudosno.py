#!/usr/bin/python3
import pandas as pd

# Load dfs
dfs = []
for path in snakemake.input:
    df = pd.read_csv(path, sep='\t')
    dfs.append(df)

# Create dict of sequences and fasta
all_pseudosno = pd.concat(dfs)
seq_dict = dict(zip(all_pseudosno.gene_id, all_pseudosno.sequence))

with open(snakemake.output.fa, 'w') as f:
    for gene_id, seq in seq_dict.items():
        f.write(f'>{gene_id}\n{seq}\n')
