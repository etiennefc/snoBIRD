#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
import collections as coll 
model = str(snakemake.wildcards.cd_predictors)
df_path = [path for path in snakemake.input if model in path][0]
df = pd.read_csv(df_path, sep='\t')
df = df[df['gene_biotype'] != 'snoRNA_pseudogene'] # don't count these as FP or FN
bar_all = snakemake.output.bar_all
bar_FN_FP = snakemake.output.bar_FN_FP
species_colors = snakemake.params.species_colors
conf_val_colors = snakemake.params.conf_value_colors

# Order species name by decreasing number of examples in the test set
d = coll.Counter(df.species_name)
species_ordered = sorted(d, key=d.get, reverse=True)

# Given a species name list, count the number of criteria in specific col of df
# that was previously filtered using species_name_list in global_col
def count_list_species(initial_df, species_name_list, global_col, criteria, specific_col):
    """
    Create a list of lists using initial_col to split the global list and
    specific_col to create the nested lists.
    """
    df_list = []

    #Sort in acending order the unique values in global_col and create a list of
    # df based on these values
    print(species_name_list)
    for val in species_name_list:
        temp_val = initial_df[initial_df[global_col] == val]
        df_list.append(temp_val)


    l = []
    for i, df in enumerate(df_list):
        temp = []
        for j, temp1 in enumerate(criteria):
            crit = df[df[specific_col] == temp1]
            crit = len(crit)
            temp.append(crit)
        l.append(temp)

    return l

# Count the FP/FN/TP/TN that are from of a given species in the test set
conf_val = ['FP', 'FN', 'TP', 'TN']

df.loc[(df['target'] == 'expressed_CD_snoRNA') & 
        (df[f'{model}_prediction'] == 'expressed_CD_snoRNA'), 
        f'{model}_confusion_value'] = 'TP'
df.loc[(df['target'] == 'other') & 
        (df[f'{model}_prediction'] == 'other'), 
        f'{model}_confusion_value'] = 'TN'
df.loc[(df['target'] == 'expressed_CD_snoRNA') & 
        (df[f'{model}_prediction'] == 'other'), 
        f'{model}_confusion_value'] = 'FN'
df.loc[(df['target'] == 'other') & 
        (df[f'{model}_prediction'] == 'expressed_CD_snoRNA'), 
        f'{model}_confusion_value'] = 'FP'
counts_per_feature = count_list_species(df, species_ordered, 'species_name',
                conf_val, f'{model}_confusion_value')

total_nb_examples = [str(sum(l)) for l in counts_per_feature] # per species in test set
xtick_labels = [label.capitalize().replace('_', ' ') for label in species_ordered]

# Convert to percent
percent = ft.percent_count(counts_per_feature)

# Create bar chart
ft.stacked_bar2(percent, xtick_labels,
                conf_val, '', 'Species name',
                f'{model}\n\nProportion of examples in test set (%)', 
                conf_val_colors, 0, 108, total_nb_examples, bar_all)

#Do the same but counting only the FP and FN
counts_per_error = count_list_species(df, species_ordered, 'species_name',
                ['FP', 'FN'], f'{model}_confusion_value')
ft.stacked_bar2(counts_per_error, xtick_labels,
                ['FP', 'FN'], '', 'Species name',
                f'{model}\n\nNumber of wrongly predicted\nexamples in test set', 
                conf_val_colors, 0, 33, total_nb_examples, bar_FN_FP)
