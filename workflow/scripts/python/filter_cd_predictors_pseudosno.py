#!/usr/bin/python3
import pandas as pd
import subprocess as sp
import re

pseudosno_df = pd.read_csv(snakemake.input.all_pseudosno, sep='\t')
snoreport = snakemake.input.snoreport
snoscan = snakemake.input.snoscan
infernal = snakemake.input.infernal_rfam 
output_df = snakemake.output.pseudosno_preds

# Filter snoreport results to remove duplicates (i.e. id ending with _[2-9]*)
snoreport_preds_dict = {}
with open(snoreport, 'r') as f:
    for line in f:
        if '>' in line:
            example_id = line.split(' ')[0].replace('>', '')
            if example_id.endswith('_1'):  # select only 1 prediction per example (drop those ending with _2)
                id = re.sub('_1$', '', example_id)
                snoreport_preds_dict[id] = 'expressed_CD_snoRNA'

# Extract which snoRNA pseudogenes were predicted as snoRNAs by snoscan
predictions_id = {}
with open(snoscan, 'r') as f:
    for line in f:
        if line.startswith('>>'):
            gene_id = line.split(' ')[1]
            if gene_id not in predictions_id.keys():
                predictions_id[gene_id] = 'expressed_CD_snoRNA'

# Get the Rfam family id with awk (pandas does not work)
sp.call("awk 'NR>2 {print $2,$3,$15}' "+infernal+" | sed '/^[^RF]/d' > temp_infernal_pseudosno.txt", shell=True)
infernal_df = pd.read_csv('temp_infernal_pseudosno.txt', sep='\s', 
                        names=['rfam_family_id', 'gene_id', 'score'], 
                        engine='python')

# Keep only the Rfam family with the highest score if multiple families per snoRNA
infernal_df = infernal_df.loc[infernal_df.groupby('gene_id')['score'].idxmax()]


# Add snoreport and snoscan prediction columns
pseudosno_df['target'] = 'other'  # not considered as an expressed_CD_snoRNA
pseudosno_df['snoreport2_prediction'] = pseudosno_df['gene_id'].map(snoreport_preds_dict).fillna('other')
pseudosno_df['snoscan_prediction'] = pseudosno_df['gene_id'].map(predictions_id).fillna('other')


# Merge rfam family id to all C/D pseudogenes
preds = pseudosno_df.merge(infernal_df[['gene_id', 'rfam_family_id']], 
                        how='left', on='gene_id')

# Create infernal_rfam_prediction col
preds['infernal_rfam_prediction'] = preds.rfam_family_id
preds['infernal_rfam_prediction'] = preds.infernal_rfam_prediction.fillna('other')
preds.loc[preds.infernal_rfam_prediction != 'other', 'infernal_rfam_prediction'] = 'expressed_CD_snoRNA'
# Save df
preds[['gene_id', 'gene_biotype', 'species_name', 'extended_211nt_sequence', 'target', 
        'snoreport2_prediction', 'snoscan_prediction', 'infernal_rfam_prediction']].to_csv(output_df, 
        sep='\t', index=False)

sp.call('rm temp_infernal_pseudosno.txt', shell=True)