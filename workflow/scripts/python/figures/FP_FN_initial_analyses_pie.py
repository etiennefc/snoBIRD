#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 


snoreport = pd.read_csv(snakemake.input.snoreport, sep='\t')
snoscan = pd.read_csv(snakemake.input.snoscan, sep='\t')
infernal_rfam = pd.read_csv(snakemake.input.infernal_rfam, sep='\t')

# Merge biotype and species info to GRU_nn predictions df
# Also exclude snoRNA_pseudogene from this analysis
gru_nn = pd.read_csv(snakemake.input.gru_nn, sep='\t')
gru_nn = gru_nn.merge(snoreport[['gene_id', 'gene_biotype', 'species_name']], 
            how='left', on='gene_id').rename(columns={'y_true': 'target', 'y_pred':'GRU_NN_prediction'})
gru_nn = gru_nn[gru_nn.target != 1]

# Convert numerical predictions to string values
gru_nn[['target', 'GRU_NN_prediction']] = gru_nn[['target', 'GRU_NN_prediction']].replace({
                                                0: 'other', 2: 'expressed_CD_snoRNA'})
pie_species = snakemake.output.pie_species
pie_neg_type = snakemake.output.pie_neg_type
error = snakemake.wildcards.error
sp_colors = snakemake.params.species_colors
ty_colors = snakemake.params.biotype_colors

# Count the FP that are part of a given species/gene_biotype
mods = ['snoreport2', 'snoscan', 'infernal_rfam', 'GRU_NN']
species = sorted(list(pd.unique(snoreport.species_name)))
biotype = sorted(list(pd.unique(snoreport.gene_biotype)))
species_colors = [sp_colors[sp] for sp in species]
biotype_colors = [ty_colors[biot] for biot in biotype]
ax_titles = []
species_count, neg_type_count = [], []
for i, df in enumerate([snoreport, snoscan, infernal_rfam, gru_nn]):
    df = df[df['gene_biotype'] != 'snoRNA_pseudogene']  # don't count these in the predictions
    if error == 'FP':
        fp = df[(df['target'] == 'other') & (df[f'{mods[i]}_prediction'] == 'expressed_CD_snoRNA')]
        err = 'false positives'
    elif error == 'FN':
        fp = df[(df['target'] == 'expressed_CD_snoRNA') & (df[f'{mods[i]}_prediction'] == 'other')]
        err = 'false negatives'
    print(fp[['species_name', 'gene_biotype', 'gene_id']])
    ax_titles.append(f'{mods[i]} ({len(fp)})')
    temp_l, temp_l2 = [], []
    for sp in species:
        nb = len(fp[fp['species_name'] == sp])
        temp_l.append(nb)
        
    species_count.append(temp_l)
    for ty in biotype:
        nb2 = len(fp[fp['gene_biotype'] == ty])
        temp_l2.append(nb2)
    neg_type_count.append(temp_l2)

# Create pie chart for species distribution in FP
ft.pie_multiple(1, len(mods), species_count, species, species_colors, ax_titles, 
                f'Species distribution across the {err}', 
                '', pie_species)
# Create pie chart for gene_biotype distribution in FP
ft.pie_multiple(1, len(mods), neg_type_count, biotype, biotype_colors, ax_titles, 
                f'Negative type distribution across the {err}', 
                '', pie_neg_type)
