#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
import re

# Load dfs
pos = pd.read_csv(snakemake.input.positives, sep='\t', names=['gene_id', 'mfe'])
neg = pd.read_csv(snakemake.input.negatives, sep='\t', names=['gene_id', 'mfe'])
biotype_df = pd.concat([pd.read_csv(path, sep='\t') for path in snakemake.input.biotype_df])

# Define params and output
color_dict = snakemake.params.biotype_colors
color_dict = {k:v for k,v in color_dict.items() if k != 'expressed_CD_snoRNA'}
target_color_dict = snakemake.params.target_colors
density_all = snakemake.output.density_all
density_negatives = snakemake.output.density_negatives
density_all_normalized = snakemake.output.density_all_normalized
density_negatives_normalized = snakemake.output.density_negatives_normalized
density_all_length = snakemake.output.density_all_length
density_negatives_length = snakemake.output.density_negatives_length

# Get length of snoRNA
pos_len, neg_len = {}, {}
with open(snakemake.input.positives_fa, 'r') as f:
    for line in f:
        if line.startswith('>'):
            id_ = line.replace('>', '').replace('\n', '')
        elif re.match(r"^A|T|C|G|U", line):
            len_ = len(line.replace('\n', ''))
            pos_len[id_] = len_

with open(snakemake.input.negatives_fa, 'r') as f:
    for line in f:
        if line.startswith('>'):
            id_ = line.replace('>', '').replace('\n', '')
        elif re.match(r"^A|T|C|G|U", line):
            len_ = len(line.replace('\n', ''))
            neg_len[id_] = len_

# Add sno_len to df and add normalized mfe
pos['sno_len'] = pos['gene_id'].map(pos_len)
neg['sno_len'] = neg['gene_id'].map(neg_len)

pos['normalized_mfe'] = pos.mfe / pos.sno_len
neg['normalized_mfe'] = neg.mfe / neg.sno_len

# Add gene_biotype to negatives
neg = neg.merge(biotype_df[['gene_id', 'gene_biotype']], how='left', on='gene_id')

# Separate pseudosno from negatives
pseudo = neg[neg['gene_biotype'] == 'snoRNA_pseudogene'].gene_id
real_neg = neg[neg['gene_biotype'] != 'snoRNA_pseudogene']

# Create density for the negatives (hue: gene_biotype) (mfe or normalized_mfe)
tempo, tempo_normalized, tempo_len = [], [], []
for biotype in color_dict.keys():
    tempo.append(real_neg[real_neg['gene_biotype'] == biotype].mfe)
    tempo_normalized.append(real_neg[real_neg['gene_biotype'] == biotype].normalized_mfe)
    tempo_len.append(real_neg[real_neg['gene_biotype'] == biotype].sno_len)

ft.density_x(tempo, 'Structure stability (kcal/mol)', 'Density', 'linear', 'Structure stability across all negative types', list(color_dict.values()), 
            list(color_dict.keys()), density_negatives)
ft.density_x(tempo_normalized, 'Normalized structure stability (kcal/mol/nt)', 'Density', 'linear', 'Normalized structure stability across all negative types', list(color_dict.values()), 
            list(color_dict.keys()), density_negatives_normalized)
ft.density_x(tempo_len, 'Length (nt)', 'Density', 'linear', 'Length across all negative types', list(color_dict.values()), 
            list(color_dict.keys()), density_negatives_length)


# Create density plot for all examples (positives, negatives and pseudosno) (mfe or normalized_mfe)
df_list = [pos.mfe, neg[neg['gene_id'].isin(pseudo)].mfe, real_neg.mfe]
ft.density_x(df_list, 'Structure stability (kcal/mol)', 'Density', 'linear', 'Structure stability of all examples of the 3 classes', 
            list(target_color_dict.values()), list(target_color_dict.keys()), density_all)

df_list = [pos.normalized_mfe, neg[neg['gene_id'].isin(pseudo)].normalized_mfe, real_neg.normalized_mfe]
ft.density_x(df_list, 'Normalized structure stability (kcal/mol/nt)', 'Density', 'linear', 'Normalized structure stability of all examples of the 3 classes', 
            list(target_color_dict.values()), list(target_color_dict.keys()), density_all_normalized)

df_list = [pos.sno_len, neg[neg['gene_id'].isin(pseudo)].sno_len, real_neg.sno_len]
ft.density_x(df_list, 'Length (nt)', 'Density', 'linear', 'Length of all examples of the 3 classes', 
            list(target_color_dict.values()), list(target_color_dict.keys()), density_all_length)



