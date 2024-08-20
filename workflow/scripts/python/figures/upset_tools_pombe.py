#!/usr/bin/python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from upsetplot import UpSet
from pybedtools import BedTool

""" Generate an upset plot to see the snoRNA prediction intersection between
    all models (3 existing predictors and SnoBIRD)."""

# Outputs
output_df = snakemake.output.df 
output_upset = snakemake.output.upset

# Load dfs and beds
bed_cols = ['chr', 'start', 'end', 'gene_name', 'score', 'strand', 'length']
df_infernal = pd.read_csv(snakemake.input.infernal_rfam, sep='\t', names=bed_cols)
df_snoscan = pd.read_csv(snakemake.input.snoscan, sep='\t', names=bed_cols)
df_snoreport = pd.read_csv(snakemake.input.snoreport, sep='\t', names=bed_cols)
df_snoBIRD = pd.read_csv(snakemake.input.snoBIRD, sep='\t', names=bed_cols)

bed_infernal = BedTool(snakemake.input.infernal_rfam)
bed_snoscan = BedTool(snakemake.input.snoscan)
bed_snoreport = BedTool(snakemake.input.snoreport)
bed_snoBIRD = BedTool(snakemake.input.snoBIRD)


# Intersect snoBIRD's prediction with that of the other tools
# Expect an overlap of at least 50 % across the other tool prediction
cols = ['chr', 'start', 'end', 'gene_name', 'score', 'strand', 'length', 
        'chr_tool', 'start_tool', 'end_tool', 'gene_name_tool', 'score_tool', 
        'strand_tool', 'length_tool']
final_cols = ['gene_id', 'snoBIRD_predicted', 'infernal_predicted', 'snoreport_predicted', 'snoscan_predicted']
snoBIRD_intersect = bed_snoBIRD.intersect(b=[bed_snoscan, bed_snoreport, bed_infernal], 
                    wa=True, wb=True, filenames=True, s=True, F=0.5).to_dataframe(names=cols)

print(snoBIRD_intersect)
ia = 0
j = 0
snoBIRD_intersect_dfs, snoBIRD = [], True
seen_snoscan, seen_snoreport, seen_infernal = [], [], []  # predictions that overlap with snoBIRD's predictions
for i, group in snoBIRD_intersect.groupby(['chr', 'start', 'end', 'strand']):
    snoBIRD_id = list(pd.unique(group['gene_name']))[0]
    tool_ids = list(pd.unique(group['gene_name_tool']))
    infernal, snoscan, snoreport = False, False, False 
    snorep, snosc = [], []
    for id_ in tool_ids:
        if id_.startswith('RF'):
            infernal = True
            seen_infernal.append(id_)
        elif id_.startswith('snoreport'):
            snoreport = True
            seen_snoreport.append(id_)
            
            snorep.append(id_)
        elif id_.startswith('snoscan'):
            snoscan = True
            seen_snoscan.append(id_)
            
            snosc.append(id_)
    if len(snorep) >1:
        j+= len(snorep) - 1
    if len(snosc) >1:
        ia+= len(snosc) -1
    snoBIRD_intersect_dfs.append(pd.DataFrame([[snoBIRD_id, snoBIRD, infernal, snoreport, snoscan]], columns=final_cols))

# Add snoBIRD-specific preds
snoBIRD_preds_overlap = pd.concat(snoBIRD_intersect_dfs)
snoBIRD_specific = []
for snobird_id in df_snoBIRD[~df_snoBIRD['gene_name'].isin(snoBIRD_preds_overlap.gene_id)]['gene_name']:
    snoBIRD_specific.append(pd.DataFrame([[snobird_id, True, False, False, False]], columns=final_cols))

snoBIRD_preds_overlap = pd.concat([snoBIRD_preds_overlap] + snoBIRD_specific)


print(snoBIRD_preds_overlap)
print(ia, j)


# Intersect snoscan's prediction with that of the other tools except snoBIRD (because already done above)
# Expect an overlap of at least 50 % across the other tool prediction
snoscan_intersect = bed_snoscan.intersect(b=[bed_snoreport, bed_infernal], 
                    wa=True, wb=True, filenames=True, s=True, F=0.5).to_dataframe(names=cols)

snoscan_intersect_dfs, snoscan = [], True
for i, group in snoscan_intersect.groupby(['chr', 'start', 'end', 'strand']):
    snoscan_id = list(pd.unique(group['gene_name']))[0]
    if snoscan_id not in list(seen_snoscan):  #not already counted in overlap with snoBIRD
        tool_ids = list(pd.unique(group['gene_name_tool']))
        infernal, snoreport = False, False
        for id_ in tool_ids:
            if id_.startswith('RF'):
                infernal = True
                seen_infernal.append(id_)
            elif id_.startswith('snoreport'):
                snoreport = True
                seen_snoreport.append(id_)
        snoscan_intersect_dfs.append(pd.DataFrame([[snoscan_id, False, infernal, snoreport, snoscan]], columns=final_cols))

# Add snoscan-specific preds
snoscan_preds_overlap = pd.concat(snoscan_intersect_dfs)
snoscan_specific = []
for snoscan_id in df_snoscan[~df_snoscan['gene_name'].isin(snoscan_preds_overlap.gene_id)]['gene_name']:
    if snoscan_id not in seen_snoscan:
        snoscan_specific.append(pd.DataFrame([[snoscan_id, False, False, False, True]], columns=final_cols))

snoscan_preds_overlap = pd.concat([snoscan_preds_overlap] + snoscan_specific)
print(snoscan_preds_overlap)


# Intersect snoreport's prediction with that of infernal only (because already done the other comparisons above)
# Expect an overlap of at least 50 % across the other tool prediction
snoreport_intersect = bed_snoreport.intersect(b=bed_infernal, 
                    wa=True, wb=True, filenames=True, s=True, F=0.5).to_dataframe(names=cols)

snoreport_intersect_dfs, snoreport = [], True
for i, group in snoreport_intersect.groupby(['chr', 'start', 'end', 'strand']):
    snoreport_id = list(pd.unique(group['gene_name']))[0]
    if snoreport_id not in list(seen_snoreport):  #not already counted in overlap with snoBIRD
        tool_ids = list(pd.unique(group['gene_name_tool']))
        infernal = False
        for id_ in tool_ids:
            if id_.startswith('RF'):
                infernal = True
                seen_infernal.append(id_)
        snoreport_intersect_dfs.append(pd.DataFrame([[snoreport_id, False, infernal, snoreport, False]], columns=final_cols))

# Add snoreport-specific preds
snoreport_preds_overlap = pd.concat(snoreport_intersect_dfs)
snoreport_specific = []
for snoreport_id in df_snoreport[~df_snoreport['gene_name'].isin(snoreport_preds_overlap.gene_id)]['gene_name']:
    if snoreport_id not in seen_snoreport:
        snoreport_specific.append(pd.DataFrame([[snoreport_id, False, False, True, False]], columns=final_cols))

snoreport_preds_overlap = pd.concat([snoreport_preds_overlap] + snoreport_specific)
print(snoreport_preds_overlap)


# Add infernal-specific preds
inf_specific = []
for id_ in df_infernal['gene_name']:
    if id_ not in seen_infernal:
        inf_specific.append(pd.DataFrame([[id_, False, True, False, False]], columns=final_cols))
inf_df = pd.concat(inf_specific)

# Concat all dfs
final_df = pd.concat([snoBIRD_preds_overlap, snoscan_preds_overlap, snoreport_preds_overlap, inf_df])
final_df['prediction_type'] = 'prediction'
final_df.to_csv('test_upset.tsv', sep='\t', index=False)
print(final_df)










# Set index based on the wanted horizontal bar chart in the upset
merged_dff = final_df.set_index(final_df.snoBIRD_predicted == True).set_index(
            final_df.infernal_predicted == True, append=True).set_index(
            final_df.snoscan_predicted == True, append=True).set_index(
            final_df.snoreport_predicted == True, append=True)
print(merged_dff)
# Upset with hue of species
upset = UpSet(merged_dff, show_counts=True, sort_by='cardinality', 
            intersection_plot_elements=0)  # disable default bar chart
#colors_species = list(sp_color_dict.values())
plt.rcParams['svg.fonttype'] = 'none'
upset.add_stacked_bars(by='prediction_type', colors=['black'], title='Number of examples') # add stacked bar chart
upset.plot()
#plt.suptitle(val)
#path = [p for p in outputs if val in p and 'species' in p][0]
#print(path)
plt.savefig('test_upset.svg', dpi=600, bbox_inches='tight')
#plt.savefig(path, dpi=600, bbox_inches='tight')


