#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
import random

num_epoch = int(snakemake.params.num_epoch)

# F1 score before the first epoch of training
dfs_0 = [list(pd.read_csv(path, sep='\t').values[0])[0] for path in snakemake.input.f1_before_train]
dfs = []
# Iterate through all epochs across folds
for i, p in enumerate(snakemake.input.f1_score_tsv):
    df = pd.read_csv(p, sep='\t')
    new_row = pd.DataFrame({'f1_score_per_epoch': [dfs_0[i]]})  # add f1 before first epoch
    df = pd.concat([new_row, df])
    df['epoch'] = [i for i in range(0, num_epoch+1)]
    dfs.append(df)

f1_score_df = pd.concat(dfs)


# Grouby the df by epoch and compute mean and stdev across folds
mean = f1_score_df.groupby('epoch').mean().reset_index().rename(columns={'f1_score_per_epoch': 'avg_f1_score'})
std = f1_score_df.groupby('epoch').std().reset_index().rename(columns={'f1_score_per_epoch': 'std_f1_score'})
merged_df = mean.merge(std, how='left', on='epoch')
print(merged_df)

# Create lineplot of average f1_score across folds and epochs +/- stdev across folds
ft.lineplot_errorbars(merged_df, 'epoch', 'avg_f1_score', 'std_f1_score', None, 'Number of epoch', 'Average F1 score across folds', 
                'Average F1 score across folds and epochs', 'grey', snakemake.output.learning_curve)
