#!/usr/bin/python3
import sys
import os
import pandas as pd
import torch
import torch.nn as nn
import random
import numpy as np
import subprocess as sp
import sklearn
from sklearn.metrics import f1_score, accuracy_score
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification, logging
#logging.set_verbosity_error()

# Load inputs and params
pretrained_model = str(sys.argv[1])  # pretrained DNABert6 model
best_hyperparams = pd.read_csv(sys.argv[4], sep='\t')
batch_size = int(best_hyperparams.batch_size.values[0])  # nb of example per batch
num_labels = 2
model_path = str(sys.argv[5])
fixed_length = sys.argv[2].split('nt.ts')[0].split('_')[-1]

# Load test data and create tensor matrix (format used by pytorch)
### Remove pseudosno from test set for now
X_test = pd.read_csv(sys.argv[2], sep='\t')
#X_test = X_test[X_test['target'] != 'CD_snoRNA_pseudogene']
y_test = pd.read_csv(sys.argv[3], sep='\t')
y_simple = y_test[y_test['gene_id'].isin(X_test.gene_id)]
y_simple = y_simple.drop(columns=['gene_id'])
y_simple['target'] = y_simple['target'].replace(2, 1)

# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Load outputs
df_metrics = sys.argv[6]
df_preds = sys.argv[7]


# Limit the number of threads that torch can spawn with (to avoid core oversubscription)
# i.e. set the number of threads to the number of CPUs requested (not all CPUs physically installed)
N_CPUS = os.environ.get("SLURM_CPUS_PER_TASK")
torch.set_num_threads(int(N_CPUS))

# Force to not parallelize tokenizing before dataloader (causes forking errors otherwise)
os.environ["TOKENIZERS_PARALLELISM"] = "false"

# Allow TF32 on matrix multiplication to speed up computations
#torch.backends.cuda.matmul.allow_tf32 = True

# Allow TF32 when using cuDNN library (GPU-related library usually automatically installed on the cluster)
#torch.backends.cudnn.allow_tf32 = True


# Transform sequence of examples in test set into kmers (6-mers)
def seq2kmer(seq, k):
    """
    Convert original str sequence to kmers (str of all kmers)
    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers

test_seqs = list(X_test[f'extended_{fixed_length}nt_sequence'])
kmer_seqs = [seq2kmer(s, 6) for s in test_seqs]

# Tokenize test data in right format and create dataloader
tokenizer = AutoTokenizer.from_pretrained(pretrained_model)  # BertTokenizerFast
eval_dataset = tokenizer(kmer_seqs, return_tensors='pt', padding=True).to(device)
eval_labels = torch.tensor(list(y_simple.target)).to(device)
eval_dataset = TensorDataset(eval_dataset.input_ids, eval_dataset.attention_mask, eval_labels)
eval_dataloader = DataLoader(eval_dataset, batch_size=batch_size)

# Load model
model = BertForSequenceClassification.from_pretrained(pretrained_model, num_labels=num_labels)
model.load_state_dict(torch.load(model_path)) 
model.to(device)
model.classifier.to(device)

# Test the model
model.eval()  # no more dropout
ev_preds, ev_labels = [], []
for i, ev_batch in enumerate(eval_dataloader):
    sp.call(f'echo EVAL BATCH {i+1}', shell=True)
    ev_input_ids, ev_attention_mask, ev_batch_labels = ev_batch
    ev_input_ids = ev_input_ids.to(device)
    ev_attention_mask = ev_attention_mask.to(device)
    ev_batch_labels = ev_batch_labels.to(device)
    with torch.no_grad():  # nor gradient computation
        # Predict the labels of eval_dataset (returns logits here)
        # where logits are the model's prediction without applying any
        # activation function (positive: more probable; negative: less probable)
        output = model(ev_input_ids, attention_mask=ev_attention_mask, labels=ev_batch_labels)
        
        # Convert logits to probabilities via the softmax activation function (in 2nd dimension of output 
        probabilities = torch.softmax(output.logits, dim=1)
        
        # Get the predicted labels from these probabilities of each class
        pred_labels = torch.argmax(probabilities, dim=1)  # get the index (0 or 1) of the highest prob
        ev_preds += pred_labels.tolist()
        ev_labels += ev_batch_labels.tolist()


# Save predictions
df = pd.DataFrame({'y_true': ev_labels, 'y_pred': ev_preds})
df = pd.concat([X_test[['gene_id']].reset_index(drop=True), df], axis=1)
df.to_csv(df_preds, sep='\t', index=False)

# Compute test prediction metrics
fscore = f1_score(ev_labels, ev_preds, average='macro')  # macro avg across 3 classes
sp.call(f'echo FSCORE {fscore}', shell=True)
sno_df = df[(df.y_true == 1) | (df.y_pred == 1)]
accuracy = accuracy_score(df.y_true, df.y_pred)
sp.call(f'echo Accuracy {accuracy}', shell=True)
TP_sno = len(sno_df[(sno_df.y_true == 1) & (sno_df.y_pred == 1)])
FP_sno = len(sno_df[(sno_df.y_true != 1) & (sno_df.y_pred == 1)])
FN_sno = len(sno_df[(sno_df.y_true == 1) & (sno_df.y_pred != 1)])
if TP_sno + FP_sno == 0:
    precision_sno = 0
else:
    precision_sno = TP_sno/(TP_sno + FP_sno)
if TP_sno + FN_sno == 0:
    recall_sno = 0
else:
    recall_sno = TP_sno/(TP_sno + FN_sno)

f1_sno = 2 * precision_sno * recall_sno / (precision_sno + recall_sno)
sp.call(f'echo PRECISION sno: {precision_sno}', shell=True)
sp.call(f'echo RECALL sno: {recall_sno}', shell=True)
sp.call(f'echo FSCORE sno: {f1_sno}', shell=True)

metrics_df = pd.DataFrame([[accuracy, fscore, precision_sno,
                            recall_sno, f1_sno]],
                            columns=['accuracy_2_classes', 'f1_score_2_classes',
                            'precision_sno', 'recall_sno', 'f1_score_sno'])
metrics_df.to_csv(df_metrics, sep='\t', index=False)

