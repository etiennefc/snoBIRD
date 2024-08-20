#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 

# Load dfs
pos = pd.read_csv(snakemake.input.positives, sep='\t', names=['gene_id', 'terminal_stem_mfe'])
neg = pd.read_csv(snakemake.input.negatives, sep='\t', names=['gene_id', 'terminal_stem_mfe'])
biotype_df = pd.concat([pd.read_csv(path, sep='\t') for path in snakemake.input.biotype_df])

# Define params and output
color_dict = snakemake.params.biotype_colors
color_dict = {k:v for k,v in color_dict.items() if k != 'expressed_CD_snoRNA'}
target_color_dict = snakemake.params.target_colors
density_all = snakemake.output.density_all
density_negatives = snakemake.output.density_negatives

# Add gene_biotype to negatives
neg = neg.merge(biotype_df[['gene_id', 'gene_biotype']], how='left', on='gene_id')

# Separate pseudosno from negatives
pseudo = neg[neg['gene_biotype'] == 'snoRNA_pseudogene'].gene_id
real_neg = neg[neg['gene_biotype'] != 'snoRNA_pseudogene']

# Create density for the negatives (hue: gene_biotype)
tempo = []
for biotype in color_dict.keys():
    tempo.append(real_neg[real_neg['gene_biotype'] == biotype].terminal_stem_mfe)
ft.density_x(tempo, 'Terminal stem stability (kcal/mol)', 'Density', 'linear', 'Terminal stem stability across all negative types', list(color_dict.values()), 
            list(color_dict.keys()), density_negatives)


# Create density plot for all examples (positives, negatives and pseudosno)
df_list = [pos.terminal_stem_mfe, neg[neg['gene_id'].isin(pseudo)].terminal_stem_mfe, real_neg.terminal_stem_mfe]
ft.density_x(df_list, 'Terminal stem stability (kcal/mol)', 'Density', 'linear', 'Terminal stem stability across all examples of the 3 classes', 
            list(target_color_dict.values()), list(target_color_dict.keys()), density_all)
