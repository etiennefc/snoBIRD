#!/usr/bin/python3
import pandas as pd
import os
import collections as coll
from pybedtools import BedTool
import subprocess as sp
from Bio import SeqIO
import utils as ut
import numpy as np
import torch
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification

# Define inputs, outputs and parameters
pretrained_model = str(snakemake.input.pretrained_model)
tokenizer_path = str(snakemake.input.tokenizer)
all_preds = list(snakemake.input.predictions)
model_path = str(snakemake.input.snoBIRD)
input_fasta_dir = str(snakemake.input.input_fasta_dir)
input_fasta = str(snakemake.input.input_fasta)

output_df = str(snakemake.output.filtered_preds)
output_df_center = str(snakemake.output.center_preds)

fixed_length = int(snakemake.params.fixed_length)
step_size = int(snakemake.params.step_size)
chunk_size = int(snakemake.params.chunk_size)
batch_size = int(snakemake.params.batch_size)
num_labels = int(snakemake.params.num_labels)
prob_threshold = float(snakemake.params.prob_threshold)
min_consecutive_windows = int(snakemake.params.min_consecutive_windows_threshold)

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
merged_blocks['start'] = merged_blocks['start'] 
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


# Deal with large blocks that can contain more than 1 snoRNA (ex: clustered 
# snoRNAs). They might not have a centered positively predicted window as there 
# are likely two or more in the same block, so not necessarily centered
# Include also step-size specific thresholds?





# Identify the center window of fixed_length in the merged block
center_window = result_df.reset_index(drop=True).copy()
center_window[['centered_start', 'centered_end']] = center_window.apply(
                    lambda row: ut.centered_window(row, fixed_length), axis=1)

df1_windows = set(df_preds.apply(lambda row: (
                row['chr'], row['start'], row['end'], row['strand']), axis=1))

# Check if rows in df2 are present in df1
""" Apply 3rd filter: if center window is not predicted as CD in the merged 
    block, don't consider that prediction (if step_size=1, thus this centered 
    window was predicted by SnoBIRD, but not as a CD) or predict on that center 
    window (if step_size>1, thus this centered window was never predicted on by 
    SnoBIRD because the larger step size passed over it)."""
center_window['is_present_in_df1'] = center_window.apply(lambda row: (
        row['chrom'], row['centered_start'], row['centered_end'], row['strand']
        ) in df1_windows, axis=1)

if step_size == 1:
    # keep only the centered window predicted as CD in those blocks
    center_window = center_window[center_window['is_present_in_df1'] == True]
elif step_size > 1:
    center_window_pos = center_window[
                                    center_window['is_present_in_df1'] == True]
    non_predicted_windows = center_window[
                                center_window['is_present_in_df1'] == False]
    if len(non_predicted_windows) > 0: # some center windows were not predicted
        # correct for 0-based bed so that we get the right sequence
        non_predicted_windows['centered_start'] = (
                                non_predicted_windows['centered_start'] - 1)
        non_predicted_windows[['chrom', 'centered_start', 'centered_end', 
                            'name', 'score', 'strand', 'len']].reset_index(
                            drop=True).to_csv('center_window_false.bed', 
                            sep='\t', index=False, header=False)
        bed_center_non_window = BedTool('center_window_false.bed')
        fasta = bed_center_non_window.sequence(fi=input_fasta, nameOnly=True, 
                                                s=True)
        preds_seqs = []
        with open(fasta.seqfn, 'r') as f:
            for line in f:
                if '>' not in line:
                    preds_seqs.append(line.strip('\n'))

        # Limit the number of threads that torch can spawn with (to avoid core 
        # oversubscription) i.e. set the num of threads to the number of CPUs 
        # requested (not all CPUs physically installed)
        N_CPUS = os.environ.get("SLURM_CPUS_PER_TASK")
        torch.set_num_threads(int(N_CPUS))

        # Force to parallelize tokenizing to increase SnoBIRD speed
        os.environ["TOKENIZERS_PARALLELISM"] = "true"

        # Allow TF32 on matrix multiplication to speed up computations
        torch.backends.cuda.matmul.allow_tf32 = True

        # Allow TF32 when using cuDNN library (GPU-related library usually 
        # automatically installed on the cluster)
        torch.backends.cudnn.allow_tf32 = True


        # Show packages versions
        sp.call(f'echo PANDAS VERSION: {pd.__version__}', shell=True)
        sp.call(f'echo TORCH VERSION: {torch.__version__}', shell=True)
        sp.call(f'echo NUMPY VERSION: {np.__version__}', shell=True)
        sp.call(f'echo TRANSFORMERS VERSION: {transformers.__version__}', 
                shell=True)
        sp.call(f'echo IS CUDA AVAILABLE?: {torch.cuda.is_available()}', 
                shell=True)



        # Define if we use GPU (CUDA) or CPU
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


        # Kmerize the sequences and predict per batch if these centered_windows 
        # are predicted as CD. If so, they will be in the final predictions
        kmer_seqs = [ut.seq2kmer(s, 6) for s in preds_seqs]

        # Tokenize data in right format and create dataloader
        tokenizer = AutoTokenizer.from_pretrained(tokenizer_path)  
        eval_dataset = tokenizer(kmer_seqs, return_tensors='pt', padding=True)
        eval_labels = torch.tensor([1] * len(kmer_seqs)).to(device, 
                                    non_blocking=True)  # fake labels
        eval_dataset = TensorDataset(eval_dataset.input_ids, 
                                    eval_dataset.attention_mask, eval_labels)
        eval_dataloader = DataLoader(eval_dataset, batch_size=batch_size)

        # Load model
        model = BertForSequenceClassification.from_pretrained(pretrained_model, 
                                                        num_labels=num_labels)
        model.load_state_dict(torch.load(model_path)) 
        model.to(device, non_blocking=True)
        model.classifier.to(device, non_blocking=True)

        # Run SnoBIRD model to predict if centered window is a C/D
        model.eval()  # put model in evaluation mode
        ev_preds, probs = [], []
        for i, ev_batch in enumerate(eval_dataloader):
            sp.call(f'echo EVAL BATCH {i+1}', shell=True)
            ev_input_ids, ev_attention_mask, ev_batch_labels = ev_batch
            ev_input_ids = ev_input_ids.to(device, non_blocking=True)
            ev_attention_mask = ev_attention_mask.to(device, non_blocking=True)
            ev_batch_labels = ev_batch_labels.to(device, non_blocking=True)
            with torch.no_grad():  # nor gradient computation
                # Predict the labels of eval_dataset (returns logits here)
                # where logits are the model's prediction without applying any
                # activation function (>0: more probable; <0: less probable)
                output = model(ev_input_ids, attention_mask=ev_attention_mask, 
                                labels=ev_batch_labels)

                # Convert logits to probabilities via the softmax activation 
                # function (in 2nd dimension of output) 
                probabilities = torch.softmax(output.logits, dim=1)
                highest_prob, _ = torch.max(probabilities, dim=1)

                # Get the predicted labels from the probabilities of each class
                pred_labels = torch.argmax(probabilities, dim=1)  # get index 
                                                # (0 or 1) of the highest prob
                ev_preds += pred_labels.tolist()
                probs += highest_prob.tolist()

        # Create df out of prediction of center windows
        df_center = pd.DataFrame({'name': list(non_predicted_windows['name']), 
                    'probability_CD': probs, 
                    'first_model_prediction': ev_preds})

        # Filter these positive predictions (1) based on probability 
        df_center = df_center[(df_center['first_model_prediction'] == 1) & (
                                df_center['probability_CD'] > prob_threshold)]
    
        # Concat the newly predicted center windows with those that were 
        # already predicted as CD in their center window
        center_window = pd.concat([
            non_predicted_windows[non_predicted_windows['name'].isin(
            df_center['name'])].reset_index(drop=True), 
            center_window_pos.reset_index(drop=True)])    
    
    else:  # all center windows were already predicted before
        center_window = center_window[
                                    center_window['is_present_in_df1'] == True]


# Save the large merged blocks    
result_df[result_df['name'].isin(center_window['name'])][
    ['chrom', 'start', 'end', 'name', 'score', 'strand', 'len']].reset_index(
    drop=True).to_csv(output_df, sep='\t', index=False, header=False)


# Add prediction id to each center window that are predicted as CD
center_window = center_window.reset_index(drop=True)
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




sp.call('rm temp_preds.bed center_window.bed center_window_false.bed', 
        shell=True)
