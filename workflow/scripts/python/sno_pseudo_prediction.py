#!/usr/bin/python3
import sys
import os
import pandas as pd
import torch
import torch.nn as nn
import random
import numpy as np
import subprocess as sp
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification, logging
from utils import seq2kmer
#logging.set_verbosity_error()

# Load inputs, params and output
pretrained_model = str(sys.argv[1])
tokenizer_path = str(sys.argv[2])  
model_path = str(sys.argv[4])
fixed_length = int(sys.argv[5])
df_preds = sys.argv[6]
batch_size = int(sys.argv[7])
num_labels = int(sys.argv[8])

# Limit the number of threads that torch can spawn with (to avoid core 
# oversubscription) i.e. set the number of threads to the number of CPUs 
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
sp.call(f'echo TRANSFORMERS VERSION: {transformers.__version__}', shell=True)
sp.call(f'echo IS CUDA AVAILABLE?: {torch.cuda.is_available()}', shell=True)



# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# Load predicted C/D by first SnoBIRD model and transform in kmers (6-mers)
first_preds_df = pd.read_csv(sys.argv[3], sep='\t', names=['chr', 'start', 
                    'end', 'gene_id', 'probability_CD', 'strand', 'block_id', 
                    f'extended_{fixed_length}nt_sequence'])
preds_seqs = list(first_preds_df[f'extended_{fixed_length}nt_sequence'])
kmer_seqs = [seq2kmer(s, 6) for s in preds_seqs]

# Tokenize data in right format and create dataloader
tokenizer = AutoTokenizer.from_pretrained(tokenizer_path)  # BertTokenizerFast
eval_dataset = tokenizer(kmer_seqs, return_tensors='pt', padding=True)
eval_labels = torch.tensor([1] * len(first_preds_df)).to(device, 
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

# Run 2nd SnoBIRD model to predict if snoRNA is expressed or a pseudogene
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

# Save predictions
df = pd.DataFrame({'probability_expressed_pseudogene': probs, 
                    'second_model_prediction': ev_preds})
df = pd.concat([first_preds_df[
                ['chr', 'start', 'end', 'gene_id', 'block_id']
                ].reset_index(drop=True), df], axis=1)
df.to_csv(df_preds, sep='\t', index=False)

