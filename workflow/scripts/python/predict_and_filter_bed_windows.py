#!/usr/bin/python3
import sys
import pandas as pd
import subprocess as sp
from Bio import SeqIO
import os
import numpy as np
from math import ceil
import time 
import torch
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification
from utils import seq2kmer
import warnings
from pybedtools import BedTool
import collections as coll
warnings.filterwarnings("ignore")

genome = str(sys.argv[1])
bed_path = str(sys.argv[2])
pretrained_model = str(sys.argv[3])
tokenizer_path = str(sys.argv[4])
model_path = str(sys.argv[5])
output_path = str(sys.argv[6])
window_size = int(sys.argv[7])
batch_size = int(sys.argv[8])
num_labels = int(sys.argv[9])
prob_thresh = float(sys.argv[10])
profile = str(sys.argv[11])
gpu = str(sys.argv[12])
final_output = str(sys.argv[13])

# Create chr_size file
chr_len_dict = {record.id: len(str(record.seq)) 
                for record in SeqIO.parse(genome, "fasta")}
with open('chr_size_'+final_output, 'w') as f:
    for k, v in chr_len_dict.items():
        f.write(f'{k}\t{v}\n')


# Load bed of windows
bed_cols = ['chr', 'start', 'end', 'gene_id', 'score', 'strand']
df = pd.read_csv(bed_path, sep='\t', header=None).reset_index(drop=True)
len_cols = len(df.columns)
new_cols = [f'OTHERCol_{i}' for i in range(1, len_cols-len(bed_cols)+1)]
df.columns = bed_cols + new_cols

# Remove first line (i.e. a header) if a header was provided unnecessarily 
# by the user
try:  # if first value of 2nd column (start) can be converted to int
    _ = int(df.iloc[0, 1])
except ValueError:
    df = df.iloc[1:].reset_index(drop=True)

# Keep only relevant columns and get the right data_type
df = df[[c for c in df.columns if not c.startswith('OTHERCol_')]]
df[['chr', 'gene_id', 'score', 'strand']] = df[
                            ['chr', 'gene_id', 'score', 'strand']].astype(str)
df[['start', 'end']] = df[['start', 'end']].astype(int)

# Drop if length of window is > 194
big_preds = list(df[df['end'] - df['start'] + 1 > 194]['gene_id'])
df = df[df['end'] - df['start'] + 1 <= 194]
if len(big_preds) > 0:
    sp.call('echo "WARNING: The following bed entries were filtered out since'+ 
        ' their length exceeds the maximal 194 nt window size of SnoBIRD:\n'+ 
        f'{big_preds}\n"', 
        shell=True)

# Drop duplicate rows based on gene_id
dup_id_dict = dict(coll.Counter(df['gene_id']))
dup_id = [dup for dup, num in dup_id_dict.items() if num > 1]
df = df.drop_duplicates('gene_id')
df = df.drop_duplicates(['chr', 'start','end', 'strand'])
if len(dup_id) > 0:
    sp.call('echo "WARNING: The following bed entries were found to be '+
            'duplicated one or multiple times based on their gene_or_locus_id'+
            ' and/or genomic location. All duplicates were filtered out and '+
            'only one entry was kept per duplicated group:\n' + f'{dup_id}\n"', 
            shell=True)

# If strand info is not present, duplicate rows and predict on both strands
no_strand = list(df[~df['strand'].isin(['+', '-'])]['gene_id'])
no_strand_rows = df[~df['strand'].isin(['+', '-'])].copy()
no_strand_rows['strand'] = '-'
no_strand_rows['gene_id'] = no_strand_rows['gene_id'].astype(str) + '_-'

df.loc[~df['strand'].isin(['+', '-']), 'gene_id'] = df['gene_id'].astype(
                                                                    str) + '_+'
df.loc[~df['strand'].isin(['+', '-']), 'strand'] = '+'
df_bed = pd.concat([df, no_strand_rows]).sort_values(
                                    ['chr', 'start']).reset_index(drop=True)

if len(no_strand) > 0:
    sp.call('echo "WARNING: The following bed entries did not have an ' +
            'appropriate strand value (neither + nor -). Therefore, '+ 
            'SnoBIRD will predict on both the + and - strand by default at '+ 
            'their respective genomic location and will create two new '+
            'gene_or_location_id for these entries (i.e. gene_id_+ and '+
            'gene_id_-): \n'+ f'{no_strand}\n"',
            shell=True)

# Create bed to extend windows to 194nt if necessary
bed = BedTool.from_dataframe(df_bed)

# Extend window and get 194 nt sequence
ext_seq_dict = {}
for i_row in bed:
    length_window = int(i_row[2]) - int(i_row[1]) 
    difference = window_size - length_window
    if difference >= 0:  # bed window is smaller than 194 nt
        # even number: split the remaining nt equally each side of the window
        if difference % 2 == 0: 
            l_extension = int(difference / 2)
            r_extension = int(l_extension)
        # odd number: split the remaining nt almost equally each side of the 
        # window (1 nt more on the 3'end)
        else: 
            l_extension = int((difference - 1) / 2)
            r_extension = int(l_extension + 1)
        extended_window = BedTool(str(i_row), from_string=True).slop(
                                    g='chr_size_'+final_output, 
                                    r=r_extension, l=l_extension)
    else:
        extended_window = BedTool(str(i_row), from_string=True)
    ext_seq = extended_window.sequence(fi=genome, nameOnly=True, s=True)
    with open(ext_seq.seqfn, 'r') as fa:
        for line in fa:
            if '>' in line:
                window_name = line.strip('\n>').replace('(+)', '').replace(
                                                                    '(-)', '')
            else:
                seq = line.strip('\n')
                ext_seq_dict[window_name] = seq

sp.call(f'echo TOTAL WINDOWS {len(ext_seq_dict.values())}', shell=True)

# Get sequences to predict on
pred_df = pd.DataFrame(list(ext_seq_dict.items()), 
            columns=['gene_id', f'extended_{window_size}nt_sequence'])
preds_seqs = list(pred_df[f'extended_{window_size}nt_sequence'])
kmer_seqs = [seq2kmer(s, 6) for s in preds_seqs]


# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Force to not parallelize tokenizing before dataloader (causes forking 
# errors otherwise)
os.environ["TOKENIZERS_PARALLELISM"] = "true"

# Limit the number of threads that torch can spawn with (to avoid core 
# oversubscription) i.e. set the number of threads to the number of CPUs 
# requested (not all CPUs physically installed)
# Do only on cluster, not in local
if profile != 'local':
    N_CPUS = os.environ.get("SLURM_CPUS_PER_TASK")
    torch.set_num_threads(int(N_CPUS))
    
    # Allow TF32 on matrix multiplication to speed up computations
    torch.backends.cuda.matmul.allow_tf32 = True
    
    # Allow TF32 when using cuDNN library (GPU-related library usually 
    # automatically installed on the cluster)
    torch.backends.cudnn.allow_tf32 = True

# Show packages versions
sp.call(f'echo PANDAS VERSION: {pd.__version__}', shell=True)
sp.call(f'echo TORCH VERSION: {torch.__version__}', shell=True)
sp.call(f'echo NUMPY VERSION: {np.__version__}', shell=True)
sp.call(f'echo TRANSFORMERS VERSION: {transformers.__version__}', shell=True)
sp.call(f'echo IS CUDA AVAILABLE?: {torch.cuda.is_available()}', shell=True)

# Load model and tokenizer 
start_time = time.time()
tokenizer = AutoTokenizer.from_pretrained(tokenizer_path)
model = BertForSequenceClassification.from_pretrained(pretrained_model, 
            num_labels=num_labels)
model.load_state_dict(torch.load(model_path)) 
model.to(device, non_blocking=True)
model.classifier.to(device, non_blocking=True)
model.eval()
end_time = time.time()
sp.call(f'echo LOAD INITIAL MODEL: {end_time-start_time}s', shell=True)

# Tokenize data in right format and create dataloader
tokenizer = AutoTokenizer.from_pretrained(tokenizer_path)  # BertTokenizerFast
eval_dataset = tokenizer(kmer_seqs, return_tensors='pt', padding=True)
eval_labels = torch.tensor([1] * len(pred_df)).to(device, 
                            non_blocking=True) # fake labels
eval_dataset = TensorDataset(eval_dataset.input_ids, 
                            eval_dataset.attention_mask, eval_labels)
eval_dataloader = DataLoader(eval_dataset, batch_size=batch_size)

# Load model
model = BertForSequenceClassification.from_pretrained(pretrained_model, 
                                                        num_labels=num_labels)
model.load_state_dict(torch.load(model_path)) 
model.to(device, non_blocking=True)
model.classifier.to(device, non_blocking=True)

# Run 1st SnoBIRD model to predict if there is a C/D snoRNA in the sequence
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
        
        # Get the predicted labels from these probabilities of each class
        pred_labels = torch.argmax(probabilities, dim=1)  # get the index 
                                                # (0 or 1) of the highest prob
        ev_preds += pred_labels.tolist()
        probs += highest_prob.tolist()

# Create df out of predictions
final_df = pd.DataFrame({'probability_CD': probs, 
                    'first_model_prediction': ev_preds})
final_df = pd.concat([df_bed[
                ['chr', 'start', 'end', 'gene_id', 'strand']
                ].reset_index(drop=True), final_df], axis=1)

# Filter based on probability
final_df = final_df[(final_df['first_model_prediction'] == 1) & (
                                final_df['probability_CD'] > prob_thresh)]

# Add relevant columns
final_df['gene_name'] = final_df['gene_id']
final_df = final_df[['chr', 'start', 'end', 'gene_id', 'probability_CD', 
                    'strand', 'gene_name']]
final_df[f'extended_{window_size}nt_sequence'] = final_df[
                                                'gene_id'].map(ext_seq_dict)
final_df.to_csv(output_path, sep='\t', index=False, header=False)


# If job is run on a HPC cluster, update the time limit for the next job 
# (shap_snoBIRD) based on the number of predictions and GPU type
if profile != 'local':
    # cannot overrule 'Unknown' gpu type, which allows user to choose a fixed 
    # time value for both genome_prediction and shap_snoBIRD rules
    if gpu != 'Unknown':
        pred_nb = len(final_df)
        # number of predictions for which SHAP values are computed per hour 
        # depending on the GPU
        gpu_rate = {'P100': 4500, 'V100': 11200, 'A100': 26000, 'H100': 62000}
        rate = gpu_rate[gpu]
        time_l = ceil(pred_nb/rate)
        # Optimize the default time limit based on the number of predicted C/D
        sp.call("scontrol update jobid=$(squeue -u $USER "+
                "--format='%.20i %.112j %.8u %.10M %.6D %R' | "+
                f"grep -w 'shap_snoBIRD.{final_output}' | "+
                "awk '{print $1}')"+f" TimeLimit={time_l}:00:00", shell=True)


sp.call('rm chr_size_'+final_output, shell=True)
