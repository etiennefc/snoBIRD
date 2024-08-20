#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 

f1_score_df = pd.read_csv(snakemake.input.f1_score_tsv, sep='\t')


# Grouby the df by epoch and compute mean and stdev
mean = f1_score_df.groupby('epoch').mean().reset_index().rename(columns={'f1_score': 'avg_f1_score'})
std = f1_score_df.groupby('epoch').std().reset_index().rename(columns={'f1_score': 'std_f1_score'})
merged_df = mean.merge(std, how='left', on='epoch')

# Create lineplot of average f1_score across folds and epochs +/- stdev across folds
ft.lineplot_errorbars(merged_df, 'epoch', 'avg_f1_score', 'std_f1_score', None, 'Number of epoch', 'Average F1 score across folds', 
                'Average F1 score across folds and epochs', 'grey', snakemake.output.learning_curve)