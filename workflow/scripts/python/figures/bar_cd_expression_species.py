#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
import glob 

colors = {k:v for k,v in snakemake.params.colors.items() if k != 'negatives'}

# Load expressed and pseudogenes dfs
expressed_paths = glob.glob(
                f'{snakemake.params.tgirt_dir}/*expressed_snoRNAs*')
pseudo_paths = glob.glob(
                f'{snakemake.params.tgirt_dir}/*pseudogene_snoRNAs*')
species = ['homo_sapiens', 'mus_musculus', 'drosophila_melanogaster', 
            'saccharomyces_cerevisiae']

exp_dfs, pseudo_dfs = [], []
for s in species:
    exp_df = pd.read_csv([p for p in expressed_paths if s in p][0], sep='\t')
    exp_df['expression_status'] = 'expressed_CD_snoRNA'
    exp_dfs.append(exp_df[['gene_id', 'expression_status', 'species_name']])
    if s in "|".join(pseudo_paths):
        pseudo_df = pd.read_csv([p for p in pseudo_paths if s in p][0], 
                                                                    sep='\t')
        pseudo_df['expression_status'] = 'snoRNA_pseudogene'
        pseudo_dfs.append(pseudo_df[['gene_id', 'expression_status', 
                                    'species_name']])

all_cd = pd.concat(exp_dfs + pseudo_dfs)



# Nb of expressed snoRNA and snoRNA pseudogene per species
exp_status = list(pd.unique(all_cd.expression_status))
counts, sno_nb = [], []
for sp in species:
    c = []
    sp_df = all_cd[all_cd['species_name'] == sp]
    sno_nb.append(len(sp_df))
    for status in exp_status:
        c.append(len(sp_df[sp_df['expression_status'] == status]))
    counts.append(c)
percent_counts = ft.percent_count(counts)

# Change long species name for shorter one
sp_name_dict = snakemake.params.species_short_name
sp_name_dict = {v: k for k,v in sp_name_dict.items()}
species_short = [sp_name_dict[s] for s in species]

# Create bar chart
ft.stacked_bar(percent_counts, species_short,
                exp_status, '', 'Species name', 
                'Proportion of C/D box snoRNAs (%)', colors, 0, 105, 
                [f'({i})' for i in sno_nb], snakemake.output.bar) 

