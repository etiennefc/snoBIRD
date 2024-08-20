#!/usr/bin/python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from upsetplot import UpSet

""" Generate an upset plot per confusion value (true positive, false positive,
    false negative or true negative) to see the intersection in the snoRNAs
    (mis)classified by all models (3 existing predictors and GRU NN)."""

sp_color_dict = snakemake.params.species_colors
biot_color_dict = snakemake.params.biotype_colors
outputs = snakemake.output


# Load predictions for each tool
snoreport = pd.read_csv(snakemake.input.snoreport, sep='\t')
snoreport = snoreport[snoreport.target != 1]
snoscan = pd.read_csv(snakemake.input.snoscan, sep='\t')
snoscan = snoscan[snoscan.target != 1]
infernal_rfam = pd.read_csv(snakemake.input.infernal_rfam, sep='\t')
infernal_rfam = infernal_rfam[infernal_rfam.target != 1]
gru_nn = pd.read_csv(snakemake.input.gru_nn, sep='\t').rename(columns={'y_true': 'target', 'y_pred':'GRU_NN_prediction'})
gru_nn[['target', 'GRU_NN_prediction']] = gru_nn[['target', 'GRU_NN_prediction']].replace({
                                                0: 'other', 2: 'expressed_CD_snoRNA'})
gru_nn = gru_nn[gru_nn.target != 1]


# Create conf_value column
mods = ['snoreport2', 'snoscan', 'infernal_rfam', 'GRU_NN']
for i, model in enumerate([snoreport, snoscan, infernal_rfam, gru_nn]):
    model.loc[(model['target'] == 'expressed_CD_snoRNA') & (model[f'{mods[i]}_prediction'] == 'expressed_CD_snoRNA'), f'conf_val_{mods[i]}'] = 'TP'
    model.loc[(model['target'] == 'expressed_CD_snoRNA') & (model[f'{mods[i]}_prediction'] != 'expressed_CD_snoRNA'), f'conf_val_{mods[i]}'] = 'FN'
    model.loc[(model['target'] == 'other') & (model[f'{mods[i]}_prediction'] != 'other'), f'conf_val_{mods[i]}'] = 'FP'
    model.loc[(model['target'] == 'other') & (model[f'{mods[i]}_prediction'] == 'other'), f'conf_val_{mods[i]}'] = 'TN'

# Merge 4 dfs
merged_df = snoreport[['gene_id', 'gene_biotype', 
                    'species_name', 'conf_val_snoreport2', 'target']].merge(snoscan[['gene_id', 'conf_val_snoscan']], 
                    how='left', on='gene_id')
merged_df = merged_df.merge(infernal_rfam[['gene_id', 'conf_val_infernal_rfam']], how='left', on='gene_id')
merged_df = merged_df.merge(gru_nn[['gene_id', 'conf_val_GRU_NN']], how='left', on='gene_id')



# Set index based on the wanted horizontal bar chart in the upset
conf_val = ['TP', 'FN', 'FP', 'TN']
for val in conf_val:
    merged_dff = merged_df.set_index(merged_df.conf_val_GRU_NN == val).set_index(
                merged_df.conf_val_snoreport2 == val, append=True).set_index(
                merged_df.conf_val_snoscan == val, append=True).set_index(
                merged_df.conf_val_infernal_rfam == val, append=True)
    # Select only examples where at least one predictor predicts it as val
    merged_dff = merged_dff[merged_dff.isin([val]).any(axis=1)]
    
    # Upset with hue of species
    upset = UpSet(merged_dff, show_counts=True, sort_by='cardinality', 
                intersection_plot_elements=0)  # disable default bar chart
    colors_species = list(sp_color_dict.values())
    plt.rcParams['svg.fonttype'] = 'none'
    upset.add_stacked_bars(by='species_name', colors=colors_species, title='Number of examples') # add stacked bar chart
    upset.plot()
    plt.suptitle(val)
    path = [p for p in outputs if val in p and 'species' in p][0]
    print(path)
    plt.savefig(path, dpi=600, bbox_inches='tight')
    #plt.show()

    # Upset with hue of gene_biotypes for negatives
    if val in ['TN', 'FP']:
        upset = UpSet(merged_dff, show_counts=True, sort_by='cardinality', 
                intersection_plot_elements=0)  # disable default bar chart
        colors_biot = list(biot_color_dict.values())
        plt.rcParams['svg.fonttype'] = 'none'
        upset.add_stacked_bars(by='gene_biotype', colors=colors_biot, title='Number of examples') # add stacked bar chart
        upset.plot()
        plt.suptitle(val)
        path = [p for p in outputs if val in p and 'biotype' in p][0]
        print(path)
        plt.savefig(path, dpi=600, bbox_inches='tight')




