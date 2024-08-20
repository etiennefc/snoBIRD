#!/usr/bin/python3
import pandas as pd
import subprocess as sp

sno_df = pd.read_csv(snakemake.input.sno_df, sep='\t')
infernal = snakemake.input.infernal

# Get the Rfam family id with awk (pandas does not work)
sp.call("awk 'NR>2 {print $2,$3,$15}' "+infernal+" | sed '/^[^RF]/d' > temp_infernal.txt", shell=True)
infernal_df = pd.read_csv('temp_infernal.txt', sep='\s', 
                        names=['rfam_family_id', 'gene_id', 'score'], 
                        engine='python')

# Keep only the Rfam family with the highest score if multiple families per snoRNA
infernal_df = infernal_df.loc[infernal_df.groupby('gene_id')['score'].idxmax()]

# Merge rfam family id to all expressed C/D
# Some snoRNAs could not be associated with a Rfam family
sno_df = sno_df.merge(infernal_df[['gene_id', 'rfam_family_id']], 
                        how='left', on='gene_id')

sno_df.to_csv(snakemake.output.df, sep='\t', index=False)

sp.call('rm temp_infernal.txt', shell=True)