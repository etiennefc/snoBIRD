#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 



# Nb of prediction overlapping annotated C/D and total of predictions
infernal = [40, 0]
snoscan = [40, 912]
snoreport = [12, 1258]
snoBIRD = [32, 305]

counts = [snoreport, snoscan, infernal, snoBIRD]
xtick_labels = ['snoreport2', 'snoscan', 'infernal_rfam', 'snoBIRD']
labels = ['Annotated C/D', 'Unannotated prediction']
colors = ['black', 'lightgrey']

# Create bar chart
ft.stacked_bar2(counts, xtick_labels,
                labels, 'Predictions overlapping annotated or unannotated\npotential C/D box snoRNAs in S. cerevisiae', 'Predictor',
                'Number of predictions', 
                colors, 0, 1350, [f'({i[0]})' for i in counts], snakemake.output.bar) 

