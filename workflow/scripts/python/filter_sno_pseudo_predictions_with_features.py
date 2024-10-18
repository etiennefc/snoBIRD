#!/usr/bin/python3
import sys
import subprocess as sp
import pandas as pd
import os
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from math import ceil
import utils as ut

# Load dfs, params and define output
sno_df = pd.read_csv(snakemake.input.sno_limits, sep='\t')
sno_df = sno_df.drop_duplicates(subset=['chr', 'start', 'end', 'strand'])
second_model_preds = pd.read_csv(snakemake.input.sno_pseudo_preds, sep='\t', 
                    names=['chr_window', 'start_window', 'end_window', 
                    'gene_id', 'block_id', 'probability_expressed_pseudogene', 
                    'second_model_prediction'], header=0)
terminal_combined_thresh = float(snakemake.params.terminal_stem_score_thresh)
sno_mfe_thresh = float(snakemake.params.normalized_sno_stability_thresh)
box_score_thresh = int(snakemake.params.box_score_thresh)
score_c_thresh = int(snakemake.params.score_c_thresh)
score_d_thresh = int(snakemake.params.score_d_thresh)
prob_sno_pseudo_thresh = float(snakemake.params.prob_second_model)
output_ = snakemake.output.final_output
output_type = snakemake.params.output_type
fixed_length = int(snakemake.params.fixed_length)

# Compute snoRNA length, normalized structure stability and 
# terminal stem combined score
sno_df['sno_length'] = sno_df['predicted_sequence'].apply(
                                lambda x: len(x))

sno_df = ut.sno_mfe(sno_df, 'predicted_sequence',
                    'sno_length')


# Get terminal stem stability/score and combined
sno_df = ut.terminal_stem(sno_df, 
                f'extended_{fixed_length}nt_sequence', extended_col=True, 
                extended_col_name='predicted_extended_sequence')

# Merge with second model predictions and round values
sno_df = sno_df.merge(second_model_preds[
                    ['gene_id', 'probability_expressed_pseudogene', 
                    'second_model_prediction']], how='left', on='gene_id')

sno_df['normalized_sno_stability'] = sno_df[
                                        'normalized_sno_stability'].round(2)
sno_df['terminal_stem_score'] = sno_df['terminal_combined'].round(2)


# Filter predictions using the feature values
int_cols = ['score_c', 'score_d', 'box_score']
float_cols = ['probability_expressed_pseudogene', 'terminal_combined', 
            'normalized_sno_stability']
sno_df[int_cols] = sno_df[int_cols].astype(int)
sno_df[float_cols] = sno_df[float_cols].astype(float)
sno_df = ut.feature_filters(sno_df, 'probability_expressed_pseudogene', 
            'second_model_prediction', 'predicted_label', 
            terminal_combined=terminal_combined_thresh, sno_mfe=sno_mfe_thresh, 
            box_score=box_score_thresh, score_c=score_c_thresh, 
            score_d=score_d_thresh, prob_thresh=prob_sno_pseudo_thresh)


# Save final output
final_cols = ['gene_id', 'chr', 'start', 'end', 'strand', 
            'probability_CD', 'box_score', 'C_MOTIF', 'C_START', 'C_END',
            'D_MOTIF', 'D_START', 'D_END', 'C_PRIME_MOTIF', 'C_PRIME_START', 
            'C_PRIME_END', 'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END', 
            'terminal_stem_score', 'normalized_sno_stability', 
            'probability_expressed_pseudogene', 'predicted_label', 
            'predicted_sequence']

int_cols = ['start', 'end', 'C_START', 'C_END', 'D_START', 'D_END', 
            'C_PRIME_START', 'C_PRIME_END', 'D_PRIME_START', 'D_PRIME_END', 
            'box_score']
df_final = sno_df[final_cols]
df_final[int_cols] = df_final[int_cols].astype(int)
df_final = df_final.drop_duplicates(subset=['chr', 'start', 'end', 'strand'])


if output_type == 'tsv':
    df_final.to_csv(output_, sep='\t', index=False)

elif output_type == 'bed':
    df_final['empty_score'] = '.'
    bed_cols = ['chr', 'start', 'end', 'gene_id', 'empty_score', 'strand', 
                'attributes']
    att_cols = ['probability_CD', 'box_score', 'C_MOTIF', 'C_START', 'C_END',
            'D_MOTIF', 'D_START', 'D_END', 'C_PRIME_MOTIF', 'C_PRIME_START', 
            'C_PRIME_END', 'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END', 
            'terminal_stem_score', 'normalized_sno_stability', 
            'probability_expressed_pseudogene', 'predicted_label', 
            'predicted_sequence']
    def create_attributes(df, col_names):
        # create an attribute column as last column of bed file
        return df.apply(lambda row: '; '.join(
                        [f"{col}={row[col]}" for col in col_names]), axis=1)

    df_final['attributes'] = create_attributes(df_final, att_cols)
    df_final[bed_cols].to_csv(output_, sep='\t', index=False, header=False)

elif output_type == 'fa':
    with open(output_, 'w') as f:
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
            terminal_score = row[1]['terminal_stem_score']
            norm_mfe = row[1]['normalized_sno_stability']
            prob_sno_pseudo = row[1]['probability_expressed_pseudogene']
            pred_label = row[1]['predicted_label']
            sequence = row[1]['predicted_sequence']
            f.write(f'>{gene_id} {chr_}:{start}-{end} ({strand}) '+
                    f'probability_CD={proba} C={c_motif}/{c_start}/{c_end} '+
                    f'D_PRIME={d_prime_motif}/{d_prime_start}/{d_prime_end} '+
                    f'C_PRIME={c_prime_motif}/{c_prime_start}/{c_prime_end} '+
                    f'D={d_motif}/{d_start}/{d_end} box_score={box_score} '+
                    f'terminal_stem_score={terminal_score} '+
                    f'normalized_sno_stability={norm_mfe} '+
                    f'predicted_label={pred_label}\n{sequence}\n')


