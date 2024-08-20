#!/usr/bin/python3
import pandas as pd
import functions as ft
import matplotlib.pyplot as plt
import functions as ft
import seaborn as sns
import matplotlib.patches as mpatches

""" Create a density plot to compare all confusion values for the 4
    intrinsic features (one plot per predictor). We exclude snoRNA 
    pseudogenes here."""

color_dict = snakemake.params.conf_value_colors
output = snakemake.output.density
feature = snakemake.wildcards.intrinsic_feature
feature_df = pd.read_csv(snakemake.input.features, sep='\t')
feature_df = feature_df[['gene_id', feature]]

# Load predictions for each tool
snoreport = pd.read_csv(snakemake.input.snoreport, sep='\t')
snoscan = pd.read_csv(snakemake.input.snoscan, sep='\t')
infernal_rfam = pd.read_csv(snakemake.input.infernal_rfam, sep='\t')
gru_nn = pd.read_csv(snakemake.input.gru_nn, sep='\t').rename(columns={'y_true': 'target', 'y_pred':'GRU_NN_prediction'})
gru_nn[['target', 'GRU_NN_prediction']] = gru_nn[['target', 'GRU_NN_prediction']].replace({
                                                0: 'other', 2: 'expressed_CD_snoRNA'})
gru_nn = gru_nn[gru_nn.target != 1]

# Create conf_value column
mods = ['snoreport2', 'snoscan', 'infernal_rfam', 'GRU_NN']
example_per_conf_val = {}
for i, model in enumerate([snoreport, snoscan, infernal_rfam, gru_nn]):
    model.loc[(model['target'] == 'expressed_CD_snoRNA') & (model[f'{mods[i]}_prediction'] == 'expressed_CD_snoRNA'), 'conf_val'] = 'TP'
    model.loc[(model['target'] == 'expressed_CD_snoRNA') & (model[f'{mods[i]}_prediction'] == 'other'), 'conf_val'] = 'FN'
    model.loc[(model['target'] == 'other') & (model[f'{mods[i]}_prediction'] == 'expressed_CD_snoRNA'), 'conf_val'] = 'FP'
    model.loc[(model['target'] == 'other') & (model[f'{mods[i]}_prediction'] == 'other'), 'conf_val'] = 'TN'

    # Get df of all examples of a given confusion_value inside dict
    for conf_val in ['TN', 'TP', 'FN', 'FP']:
        df_temp = model[model['conf_val'] == conf_val]
        example_id_list = df_temp['gene_id'].to_list()
        df = feature_df[feature_df['gene_id'].isin(example_id_list)]
        if mods[i] not in example_per_conf_val.keys():
            example_per_conf_val[mods[i]] = {conf_val: df}
        else:
            example_per_conf_val[mods[i]][conf_val] = df



# Create density plot
colors = [color_dict['TN'], color_dict['TP'], color_dict['FN'], color_dict['FP']]
fig, axes = plt.subplots(1, 4, figsize=(45, 14))
ax = axes.flatten()

rc = {'ytick.labelsize': 60, 'xtick.labelsize': 60}
plt.rcParams.update(**rc)
plt.rcParams['svg.fonttype'] = 'none'
for i, model_name in enumerate(mods):
    dfs = example_per_conf_val[model_name]
    df_list = [dfs['TN'][feature], dfs['TP'][feature], dfs['FN'][feature], dfs['FP'][feature]]
    for j, df in enumerate(df_list):
        sns.kdeplot(df, fill=True, color=colors[j], ax=ax[i])
    ax[i].set_title(model_name, fontdict={'fontsize': 35})
    ax[i].set_xlabel(feature, fontdict={'fontsize': 35})
    ax[i].set_ylabel('Density', fontdict={'fontsize': 35})
    ax[i].set_xscale('linear')
    ax[i].spines['right'].set_linewidth(0)
    ax[i].spines['top'].set_linewidth(0)
legend_list = []
for i, crit in enumerate(['TN', 'TP', 'FN', 'FP']):
    legend_element = mpatches.Patch(color=colors[i], label=crit)
    legend_list.append(legend_element)
plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0,1.1),
            fontsize=20)
plt.savefig(output, bbox_inches='tight', dpi=500)
