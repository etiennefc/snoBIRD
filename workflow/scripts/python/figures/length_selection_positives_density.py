#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 


df = pd.read_csv(snakemake.input.positives, sep='\t')
output = snakemake.output.density
percentile_colors= snakemake.params.percent_colors

# Create length column
df['len'] = df.end - df.start + 1

# Get snoRNA length at given percentile
length_thresholds, percentile_labels = [], []
for i in range(75, 100, 5):
    percentile = i / 100
    length_at_percentile = round(df.len.quantile(q=percentile, 
                interpolation='nearest'))  # choose the nearest of the two points
                                           # if the quantile is between two points
    length_thresholds.append(length_at_percentile)
    percentile_labels.append(f'{i}th percentile ({length_at_percentile} nt)')


ft.density_percentile(df.len, 'Expressed C/D box snoRNA length (nt)', 
                    'Density', '', output, length_thresholds, 
                    percentile_colors, percentile_labels, color='grey')
