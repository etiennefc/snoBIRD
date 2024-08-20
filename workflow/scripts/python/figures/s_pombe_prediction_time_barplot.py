#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, recall_score, matthews_corrcoef

# On S. pombe gneome in min; snoBIRD was used with step size=5
pred_time = [['snoreport2', 28.4], 
                ['snoscan', 166.0833333], 
                ['infernal_rfam', 16.9166667],
                ['snoBIRD', 52.883333]] 
df = pd.DataFrame(pred_time, columns=['predictor', 'prediction_time'])

print(df)

output = snakemake.output.barplot
color_dict = {'snoreport2': '#fda47a', 'snoscan': '#80b1d3', 'infernal_rfam': '#d73027', 'snoBIRD': '#66a61e'}
print(color_dict)

plt.rcParams['svg.fonttype'] = 'none'
rc = {'ytick.labelsize': 25, 'xtick.labelsize': 25}
plt.rcParams.update(**rc)
fig, ax = plt.subplots(1, 1, figsize=(12, 10))
sns.barplot(df, x='predictor', y='prediction_time', palette=color_dict)
plt.title('Prediction time on S. pombe genome', fontsize=30, y=1.01)
ax.set_xlabel('Predictor', fontdict={'fontsize': 30})
ax.set_ylabel('Prediction time (min)', fontdict={'fontsize': 30})
plt.savefig(snakemake.output.barplot, dpi=600, bbox_inches='tight')