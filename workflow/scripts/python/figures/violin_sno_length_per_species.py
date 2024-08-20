#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 

fixed_length = snakemake.wildcards.fixed_length
expressed_cd = pd.read_csv(snakemake.input.expressed_cd, sep='\t')
expressed_cd['gene_biotype'] = 'expressed_CD_snoRNA'
human_pseudo = pd.read_csv(snakemake.input.human_pseudo, sep='\t')
mouse_pseudo = pd.read_csv(snakemake.input.mouse_pseudo, sep='\t')
drosophila_pseudo = pd.read_csv(snakemake.input.drosophila_pseudo, sep='\t')
all_pseudo = pd.concat([human_pseudo, mouse_pseudo, drosophila_pseudo])
all_pseudo = all_pseudo[['gene_id', 'chr', 'strand', 'start', 'end', 'sequence', 'species_name',
       f'extended_{fixed_length}nt_sequence']]
all_pseudo['gene_biotype'] = 'snoRNA_pseudogene'

short_sp_names = snakemake.params.short_sp_names
short_sp_names = {v: k for k,v in short_sp_names.items()}

# Concat expressed and pseudo
df = pd.concat([expressed_cd, all_pseudo])
df['snoRNA_length'] = df.end.astype(int) - df.start.astype(int) + 1
df['species_name'] = df['species_name'].map(short_sp_names)



# Create the violin plot
colors = snakemake.params.target_colors
colors = {k:v for k,v in colors.items() if k != 'negatives'}
sp_order = ['H_sapiens', 'M_mulatta', 'M_musculus', 'O_anatinus', 
             'G_gallus', 'C_elegans', 'D_melanogaster', 
            'S_cerevisiae', 'A_fumigatus', 'C_albicans', 'N_crassa', 
            'A_thaliana', 'O_sativa', 'O_tauri', 
            'D_discoideum', 'G_lamblia', 'L_major', 'T_thermophila']

plt.rcParams['svg.fonttype'] = 'none'
plt.subplots(1, 1, figsize=(15, 10))
ax = sns.violinplot(data=df, x='species_name', y='snoRNA_length', hue='gene_biotype', palette=colors, order=sp_order, fill=None)
sns.stripplot(data=df, x='species_name', y='snoRNA_length', hue='gene_biotype', palette=colors, jitter=True, dodge=True, 
            order=sp_order, ax=ax, size=2, linewidth=0.5, edgecolor='black')
ax.hlines(y=int(fixed_length) - 30, xmin=0, xmax=len(pd.unique(df.species_name)) - 1,
            linestyles='dashed', colors='black')
plt.xticks(rotation=90)
plt.title('Length distribution of expressed C/D and snoRNA pseudogenes across eukaryotes', fontsize=20)
plt.xlabel('Species name', fontsize=20)
plt.ylabel('SnoRNA length (nt)', fontsize=20)
plt.savefig(snakemake.output.violin, bbox_inches='tight', dpi=600)
