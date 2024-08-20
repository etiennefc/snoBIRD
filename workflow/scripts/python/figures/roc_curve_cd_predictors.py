#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
from sklearn import metrics


snoreport = pd.read_csv(snakemake.input.snoreport, sep='\t')
snoscan = pd.read_csv(snakemake.input.snoscan, sep='\t')
infernal_rfam = pd.read_csv(snakemake.input.infernal_rfam, sep='\t')
output = snakemake.output.roc
color_dict = snakemake.params.colors

# Merge dfs
df = snoreport.merge(snoscan[['gene_id', 'snoscan_prediction']], 
                        how='left', on='gene_id')
df = df.merge(infernal_rfam[['gene_id', 'infernal_rfam_prediction']], 
                        how='left', on='gene_id')

# Get the true/false positive rate (tpr/fpr) of all models
def get_metrics(df_, test_col, pred_col):
    df2 = df_.copy()
    df2[[test_col, pred_col]] = df2[[test_col, pred_col]].replace(
                                {'expressed_CD_snoRNA': 1, 'other': 0})
    fpr, tpr, _ = metrics.roc_curve(df2[test_col],  df2[pred_col])
    return fpr, tpr

f, t = get_metrics(df, 'target', 'snoscan_prediction')
plt.plot(f, t)
plt.show()
