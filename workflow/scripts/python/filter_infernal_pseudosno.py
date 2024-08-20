#!/usr/bin/python3
import pandas as pd
import subprocess as sp

dfs = []
for path in snakemake.input:
    if 'infernal' not in path:
        df = pd.read_csv(path, sep='\t')
        dfs.append(df)

sno_df = pd.concat(dfs)
infernal = snakemake.input.infernal

# Get the Rfam family id with awk (pandas does not work)
sp.call("awk 'NR>2 {print $2,$3,$15}' "+infernal+" | sed '/^[^RF]/d' > temp_infernal_pseudogene.txt", shell=True)
infernal_df = pd.read_csv('temp_infernal_pseudogene.txt', sep='\s', 
                        names=['rfam_family_id', 'gene_id', 'score'], 
                        engine='python')

# Keep only the Rfam family with the highest score if multiple families per snoRNA pseudogene
infernal_df = infernal_df.loc[infernal_df.groupby('gene_id')['score'].idxmax()]

# Merge rfam family id to all snoRNA pseudogenes (across mouse, drosophila and human)
# Some snoRNAs could not be associated with a Rfam family
sno_df = sno_df.merge(infernal_df[['gene_id', 'rfam_family_id']], 
                        how='left', on='gene_id')


print(sno_df)
sno_df.to_csv(snakemake.output.df, sep='\t', index=False)

sp.call('rm temp_infernal_pseudogene.txt', shell=True)