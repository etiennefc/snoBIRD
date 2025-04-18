#!/usr/bin/python3
import pandas as pd
import collections as coll
import re
import regex
import numpy as np
import utils as ut
import sys
import warnings
warnings.filterwarnings("ignore")

# Load df and variables
len_c_box, len_d_box = 7, 4
fixed_length = int(sys.argv[4])
preds = pd.read_csv(str(sys.argv[2]), sep='\t', names=
        ['chr_window', 'start_window', 'end_window', 'gene_id', 'score', 
        'strand_window', 'block_name', f'extended_{fixed_length}nt_sequence'], 
        dtype={'chr_window': 'str'})
shap_df = pd.read_csv(str(sys.argv[1]), sep='\t').rename(
                    columns={'probability': 'probability_CD'})
shap_df = shap_df[shap_df['predicted_label'] != 'Other']

half_window = int((fixed_length + 2) / 2)  # the +2 is to add for the [CLS] and 
                            # [SEP] tokens at the start and end of the sequence

if (len(shap_df) == 0) | (len(preds) == 0):  # no snoRNA was predicted in input
	final_cols = ['gene_id', 'chr', 'start', 'end', 'strand', 'probability_CD', 
        'box_score', 'C_MOTIF', 'C_START', 'C_END', 'score_c', 'D_MOTIF', 
        'D_START', 'D_END', 'score_d', 'C_PRIME_MOTIF', 'C_PRIME_START', 
        'C_PRIME_END', 'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END', 
        'predicted_sequence', 'predicted_extended_sequence', 
        f'extended_{fixed_length}nt_sequence']
	df = pd.DataFrame(columns=final_cols)
else:
	shap_df = shap_df.merge(preds[['gene_id', 
	                        f'extended_{fixed_length}nt_sequence']], how='left', 
	                        on='gene_id')
	shap_cols = [i for i in shap_df.columns if i.startswith('SHAP_')]

	# The smallest reliably annotated C/D snoRNA is 50 nt long (NR_145814, 
	# a C/D snoRNA in human). To not miss any snoRNA, we define the minimal length 
	# to find the C or D box to be +-15 nt from the center of the predicted window 
	# as it is in this snoRNA
	min_box_dist = int(sys.argv[5])  # default: 15
	# Flanking nt extending after the snoRNA start/end are minimally of 15 nt, so 
	# no C or D box should be found in the first and last 20 nt of the window 
	# (15 nt + 5 nt within the snoRNA to reach C_start or D_end)
	flanking_nt = int(sys.argv[6])  # default: 20

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
	df['chr'] = df['chr'].astype(str)

	# Get predicted sequence as well as predicted extended sequence
	# These sequence will be used to determine the snoRNA structure stability 
	# and terminal stem stability
	df = df.merge(preds[['gene_id', f'extended_{fixed_length}nt_sequence']], 
	                how='left', on='gene_id')
	df['predicted_sequence'] = df.apply(ut.get_seq, axis=1)

	# find extended seq (15 nt flanking the snoRNA) to compute the 
	# terminal stem stability
	df['predicted_extended_sequence'] = df.apply(lambda row: ut.get_seq(
	                                                row, extended=True), axis=1)

	# Correct box position based on the predicted sequence, not the extended window
	df = df.apply(ut.correct_box_pos, axis=1, 
	                motifs=['C', 'D', 'C_PRIME', 'D_PRIME'])

	df = df.drop_duplicates(subset=['chr', 'start', 'end', 'strand'])

# Save df
final_cols = ['gene_id', 'chr', 'start', 'end', 'strand', 'probability_CD', 
        'box_score', 'C_MOTIF', 'C_START', 'C_END', 'score_c', 'D_MOTIF', 
        'D_START', 'D_END', 'score_d', 'C_PRIME_MOTIF', 'C_PRIME_START', 
        'C_PRIME_END', 'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END', 
        'predicted_sequence', 'predicted_extended_sequence', 
        f'extended_{fixed_length}nt_sequence']

df[final_cols].to_csv(str(sys.argv[3]), sep='\t', index=False)

