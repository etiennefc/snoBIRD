#!/usr/bin/python3
import pandas as pd
import os
import collections as coll
from pybedtools import BedTool
import subprocess as sp
from Bio import SeqIO
import utils as ut

# Define inputs, parameters and outputs
input_fasta = snakemake.input.input_fasta
input_fasta_dir = snakemake.input.input_fasta_dir
all_preds = snakemake.input.predictions
fixed_length = int(snakemake.params.fixed_length)
strand = snakemake.params.strand 
step_size = snakemake.params.step_size
chunk_size = snakemake.params.chunk_size
prob_threshold = float(snakemake.params.prob_threshold)
min_consecutive_windows = int(
                            snakemake.params.min_consecutive_windows_threshold)
output_df = snakemake.output.filtered_preds
output_df_center = snakemake.output.center_preds



# Concat all predictions of first SnoBIRD model into 1 df (preds from all 
# chr and/or chunks)
""" Apply 1st filter on probability."""
dfs = []
for p in all_preds:
    df = pd.read_csv(p, sep='\t')
    df = df.rename(columns={'probs': 'probability_first_model'})
    
    df = df[df['probability_first_model'] > prob_threshold]
    dfs.append(df)
df_preds = pd.concat(dfs)


# For chunks of chr, rectify genomic location of positive windows based on 
# the number of previous chunks
# Preds on - strand have already been corrected during genome_prediction 
chunk_fa = [path for path in os.listdir(
            input_fasta_dir+'/fasta/') if '_chunk_' in path]
rectify_dict = {}
if len(chunk_fa) > 0:
    for chunk in chunk_fa:
        with open(input_fasta_dir+'/fasta/'+chunk, 'r') as f:
            first_line = f.readline().replace('>', '').replace('\n', '')
            chunk_name = first_line.split('  ')[0]
            chunk_number = int(chunk_name.split('_')[-1])
            prev_size = int(first_line.split(' ')[2].replace('prev_size=', ''))
            total_size = int(first_line.split(' ')[3].replace(
                            'total_size=', ''))
            rectify_dict[chunk_name] = (prev_size, total_size)


def rectify_pos(row, chr_dict):
    # Rectify position based on the number of nt before a given chunk
    if row['chr'] in chr_dict.keys():
        nt_number = chr_dict[row['chr']][0]
        row['start'] = row['start'] + nt_number 
        row['end'] = row['end'] + nt_number
        row['chr'] = row['chr'].split('_chunk_')[0]
    return row

df_preds = df_preds.apply(rectify_pos, axis=1, chr_dict=rectify_dict)
df_preds['chr'] = df_preds['chr'].astype(str)
df_preds = df_preds.sort_values(by=['chr', 'start', 'strand'])




# Calculate difference in start positions between consecutive rows
# It creates a NA for each first row of a group in the groupby
# This is to further filter on the nb of consecutive windows
df_preds['diff'] = df_preds.groupby(['chr', 'strand'])['start'].diff()


# Create boolean mask to highlight where diff >1 (
# i.e. a new stretch of consecutive nt is created)
# Use cumulative sum to assign a different number to each stretch 
# of consecutive nts
mask = (df_preds['start'].shift(-1) == df_preds['start'] + 1)

# Manually change these NA to 2 because they should be counted in the 
# block length (first element of the block)
df_preds.loc[(df_preds['diff'].isna()) & (mask), 'diff'] = 2  

# Manually change the remaining NA that are not part of stretches so 
# that they are assigned a different stretch_id
df_preds.loc[df_preds['diff'].isna(), 'diff'] = 2  
df_preds['stretch_id'] = (df_preds['diff'] > 1).cumsum()
df_preds['stretch_id'] =  'block_' + df_preds['stretch_id'].astype(str)
df_preds = df_preds[['chr', 'start', 'end', 'stretch_id', 'strand', 
                    'probability_first_model']]
df_preds.to_csv('temp_preds.bed', sep='\t', index=False, header=False)
preds_bed = BedTool('temp_preds.bed')

# Merge consecutive positive preds into one entry and create bed out of 
# all blocks this is done via groupby over the stretch_id column and keeping 
# the lowest/highest start/end of all windows in that block
# and by averaging the probabilities over all windows in that block
merged_blocks = preds_bed.groupby(g=[4], c=[2, 3, 4, 6, 5, 1], 
                o=['min', 'max', 'first', 'mean', 'first', 'first'])
# reorder columns of bed in right order
merged_blocks = merged_blocks.cut([6, 1, 2, 3, 4, 5]).to_dataframe()  
# correct bed coordinates (since they are 1-based)
merged_blocks['start'] = merged_blocks['start'] + 1
merged_blocks = merged_blocks.sort_values(by =
                                ['chrom', 'strand', 'start', 'end'])

# Get length of block
merged_blocks['len'] = merged_blocks.end - merged_blocks.start + 1


# Merge blocks that overlap reciprocally to at least 50%
# Apply the merging function row-wise
merged_rows = []
current_row = None

for i in range(len(merged_blocks)):
    if i < len(merged_blocks) - 1:
        if current_row is None:  # first row or after previous block has ended
            current_row = merged_blocks.iloc[i].copy()
        next_row = merged_blocks.iloc[i + 1].copy()
        merged_row = ut.merge_rows(current_row, next_row)
        
        # Current_row (can be a block) is merged with next row/block
        if merged_row is not None:  
            current_row = pd.Series(merged_row)
        else:  # no merge
            merged_rows.append(current_row)
            current_row = None
    else:  # deal with last row if it's not included in last block
        last_row = merged_blocks.iloc[i]
        # Last block has ended and doesn't include last row
        if current_row is None:  
            merged_rows.append(last_row)
        # otherwise it has already been merged with before last row/block



# Create a DataFrame from the merged rows
result_df = pd.DataFrame(merged_rows)
""" Apply 2nd filter: length of merged blocks"""
#result_df = result_df[(result_df['len'] >= 204) & (result_df['len'] <= 264)]
result_df = result_df[result_df['len'] >= (
                                    fixed_length + min_consecutive_windows)]
result_df.reset_index(drop=True).reset_index(drop=True)[
    ['chrom', 'start', 'end', 'name', 'score', 'strand', 'len']
    ].to_csv(output_df, sep='\t', index=False, header=False)

# Deal with large blocks that can contain more than 1 snoRNA (ex: clustered 
# snoRNAs). They might not have a centered positively predicted window as there 
# are likely two or more in the same block, so not necessarily centered
# Include also step-size specific thresholds?





# Select the center window of fixed_length in the merged block
center_window = result_df.reset_index(drop=True).copy()
center_window[['centered_start', 'centered_end']] = center_window.apply(
                    lambda row: ut.centered_window(row, fixed_length), axis=1)

df1_windows = set(df_preds.apply(lambda row: (
                row['chr'], row['start'], row['end'], row['strand']), axis=1))

# Check if rows in df2 are present in df1
""" Apply 3rd filter: if center window is not predicted as CD in the merged 
    block, don't consider that prediction."""
center_window['is_present_in_df1'] = center_window.apply(lambda row: (
        row['chrom'], row['centered_start'], row['centered_end'], row['strand']
        ) in df1_windows, axis=1)
center_window = center_window[center_window['is_present_in_df1'] == True]

# Add prediction id
center_window['index_'] = center_window.index + 1
center_window['prediction_id'] = 'CD_' + center_window.index_.astype(str)
cols_bed = ['chrom', 'centered_start', 'centered_end', 'prediction_id', 
            'score', 'strand', 'name']
center_window[cols_bed].reset_index(drop=True).to_csv('center_window.bed', 
                        sep='\t', index=False, header=False)
bed_center_window = BedTool('center_window.bed')

# Get sequence of the centered predicted CD window
fasta = bed_center_window.sequence(fi=input_fasta, nameOnly=True, s=True)
d = {}
with open(fasta.seqfn, 'r') as f:
    for line in f:
        if '>' in line:
            sno_id = line.replace('>', '').replace('\n', '')
            sno_id = sno_id.replace('(+)', '').replace('(-)', '')
        else:
            seq = line.replace('\n', '')
            d[sno_id] = seq

center_window[f'extended_{fixed_length}nt_sequence'] = center_window[
                                                        'prediction_id'].map(d)
center_window[cols_bed + [f'extended_{fixed_length}nt_sequence']].to_csv(
                        output_df_center, sep='\t', index=False, header=False)




sp.call('rm temp_preds.bed center_window.bed', shell=True)
