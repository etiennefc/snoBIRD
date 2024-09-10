#!/usr/bin/python3
import pandas as pd
import collections as coll
import re
import regex
import numpy as np
import utils as ut

# Load df and variables
len_c_box, len_d_box = 7, 4
fixed_length = int(snakemake.params.fixed_length)
output_type = str(snakemake.params.output_type)
preds = pd.read_csv(snakemake.input.preds, sep='\t', names=
        ['chr_window', 'start_window', 'end_window', 'gene_id', 'score', 
        'strand_window', 'block_name', f'extended_{fixed_length}nt_sequence'])
shap_df = pd.read_csv(snakemake.input.shap_df, sep='\t').rename(
                    columns={'probability': 'probability_CD'})
shap_df = shap_df[shap_df['predicted_label'] != 'Other']

half_window = int((fixed_length + 2) / 2)  # the +2 is to add for the [CLS] and 
                            # [SEP] tokens at the start and end of the sequence
shap_df = shap_df.merge(preds[['gene_id', 
                        f'extended_{fixed_length}nt_sequence']], how='left', 
                        on='gene_id')
shap_cols = [i for i in shap_df.columns if i.startswith('SHAP_')]

# The smallest reliably annotated C/D snoRNA is 50 nt long (NR_145814, 
# a C/D snoRNA in human). To not miss any snoRNA, we define the minimal length 
# to find the C or D box to be +-15 nt from the center of the predicted window 
# as it is in this snoRNA
min_box_dist = snakemake.params.min_box_dist  # default: 15
# Flanking nt extending after the snoRNA start/end are minimally of 15 nt, so 
# no C or D box should be found in the first and last 20 nt of the window 
# (15 nt + 5 nt within the snoRNA to reach C_start or D_end)
flanking_nt = snakemake.params.flanking_nt  # default: 20

# Thus, the ranges in which a C and D boxes can be searched are defined below
## + 1 to account for the [CLS] ##
C_range = range(flanking_nt + 1, half_window-min_box_dist)
D_range = range(half_window+min_box_dist, fixed_length - flanking_nt + 1)



# Based on SHAP values, find the C and D box (they're usually located at 
# peaks of SHAP values). Then find C' and D' boxes and the overall box score
df = ut.find_all_boxes(shap_df, fixed_length, shap_cols, C_range,
                            D_range, flanking_nt)

# Get genomic location of snoRNA based on C and D positions
df = df.merge(preds[['gene_id', f'extended_{fixed_length}nt_sequence',
                'chr_window', 'strand_window', 'start_window', 'end_window']], 
                how='left', on='gene_id')

df[['chr', 'start', 'end', 'strand']] = df.apply(ut.get_sno_location, axis=1)

# Get predicted sequence 
df['predicted_sequence'] = df.apply(ut.get_seq, axis=1)

# Correct box position based on the predicted sequence, not the extended window
df = df.apply(ut.correct_box_pos, axis=1, 
                motifs=['C', 'D', 'C_PRIME', 'D_PRIME'])


# Save df or fasta depending on the user-defined output type
final_cols = ['gene_id', 'chr', 'start', 'end', 'strand', 'probability_CD', 
            'box_score', 'C_MOTIF', 'C_START', 'C_END', 'D_MOTIF', 'D_START', 
            'D_END', 'C_PRIME_MOTIF', 'C_PRIME_START', 'C_PRIME_END', 
            'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END', 
            'predicted_sequence']
int_cols = ['start', 'end', 'C_START', 'C_END', 'D_START', 'D_END', 
            'C_PRIME_START', 'C_PRIME_END', 'D_PRIME_START', 'D_PRIME_END', 
            'box_score']
df_final = df[final_cols]
df_final[int_cols] = df_final[int_cols].astype(int)

if output_type == 'tsv':
    df_final.to_csv(snakemake.output.minimal_output, sep='\t', index=False)
if output_type == 'fa':
    with open(snakemake.output.minimal_output, 'w') as f:
        for row in df_final.iterrows():
            gene_id = row[1]['gene_id']
            chr_ = row[1]['chr']
            start = row[1]['start']
            end = row[1]['end']
            strand = row[1]['strand']
            proba = row[1]['probability_CD']
            box_score = row[1]['box_score']
            c_motif = row[1]['C_MOTIF']
            c_start = row[1]['C_START']
            c_end = row[1]['C_END']
            d_motif = row[1]['D_MOTIF']
            d_start = row[1]['D_START']
            d_end = row[1]['D_END']
            c_prime_motif = row[1]['C_PRIME_MOTIF']
            c_prime_start = row[1]['C_PRIME_START']
            c_prime_end = row[1]['C_PRIME_END']
            d_prime_motif = row[1]['D_PRIME_MOTIF']
            d_prime_start = row[1]['D_PRIME_START']
            d_prime_end = row[1]['D_PRIME_END']
            sequence = row[1]['predicted_sequence']
            f.write(f'>{gene_id} {chr_}:{start}-{end} ({strand}) '+
                    f'probability_CD={proba} C={c_motif}/{c_start}/{c_end} '+
                    f'D_PRIME={d_prime_motif}/{d_prime_start}/{d_prime_end} '+
                    f'C_PRIME={c_prime_motif}/{c_prime_start}/{c_prime_end} '+
                    f'D={d_motif}/{d_start}/{d_end} box_score={box_score}'+
                    f'\n{sequence}\n')
            



