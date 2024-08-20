#!/usr/bin/python3
import sys
import pandas as pd
import torch
import torch.nn as nn
import random
import numpy as np
import subprocess as sp
import sklearn
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, accuracy_score, precision_recall_fscore_support
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification, logging
#logging.set_verbosity_error()

# Load params
pretrained_model = sys.argv[1]  # pretrained DNABert6 model
fold_num = str(sys.argv[2])  # training fold number
rs = int(sys.argv[3])  # random_state
fixed_length = sys.argv[4].split('nt.ts')[0].split('_')[-1]
best_hyperparams = pd.read_csv(sys.argv[6], sep='\t')
batch_size = int(best_hyperparams.batch_size.values[0])  # nb of example per batch

# Load inputs
X_train = pd.read_csv(sys.argv[4], sep='\t')
y_train = pd.read_csv(sys.argv[5], sep='\t')

# Select only the C/D pseudogenes (0) and expressed C/D (1) for the training
X_train = X_train[X_train['target'] != 'other'].reset_index(drop=True)
y_train = y_train[y_train['gene_id'].isin(X_train.gene_id)]
y_simple = y_train.drop(columns=['gene_id']).reset_index(drop=True)
y_simple['target'] = y_simple['target'].replace(1, 0)
y_simple['target'] = y_simple['target'].replace(2, 1)


# Get path of outputs
output_loss = sys.argv[7]
output_f1 = sys.argv[8]

# Show packages versions
sp.call(f'echo PANDAS VERSION: {pd.__version__}', shell=True)
sp.call(f'echo TORCH VERSION: {torch.__version__}', shell=True)
sp.call(f'echo NUMPY VERSION: {np.__version__}', shell=True)
sp.call(f'echo SKLEARN VERSION: {sklearn.__version__}', shell=True)
sp.call(f'echo TRANSFORMERS VERSION: {transformers.__version__}', shell=True)
sp.call(f'echo IS CUDA AVAILABLE?: {torch.cuda.is_available()}', shell=True)

# Set reproducible randomness
torch.backends.cudnn.deterministic = True
random.seed(rs)
torch.manual_seed(rs)
torch.cuda.manual_seed(rs)
np.random.seed(rs)

# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Transform sequence of examples in training set into kmers (6-mers)
def seq2kmer(seq, k):
    """
    Convert original str sequence to kmers (str of all kmers)
    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers

train_seqs = list(X_train[f'extended_{fixed_length}nt_sequence'])
kmer_seqs = [seq2kmer(s, 6) for s in train_seqs]

# Load pre-trained DNABERT model
tokenizer = AutoTokenizer.from_pretrained(pretrained_model)  # BertTokenizerFast
model = BertForSequenceClassification.from_pretrained(pretrained_model, num_labels=2)  # BertModel
print(model)
model.to(device)
model.classifier.to(device)

# Define loss function
loss_fn = torch.nn.CrossEntropyLoss(weight=torch.tensor([1/20, 1]).to(device))

# Train over given fold (fold_num) in stratified 10-fold CV
skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=rs)
fold_dict = {str(fold_index+1): [train_index, test_index]
            for fold_index, (train_index, test_index) in
            enumerate(skf.split(kmer_seqs, y_simple))}

train_index = fold_dict[fold_num][0]
test_index = fold_dict[fold_num][1]


# Load eval (test) dataset
x_test = [k for k in kmer_seqs if kmer_seqs.index(k) in test_index]
Y_test = [y_simple.loc[i, 'target'] for i in test_index]



# Get prediction performance without any training
eval_dataset_0 = tokenizer(x_test, return_tensors='pt').to(device)
eval_labels_0 = torch.tensor(Y_test).to(device)
eval_dataset_0 = TensorDataset(eval_dataset_0.input_ids, eval_dataset_0.attention_mask, eval_labels_0)
eval_dataloader_0 = DataLoader(eval_dataset_0, batch_size=batch_size)
ev_preds_0, ev_labels_0 = [], []
total_loss = 0.0
epoch_f_scores, epoch_loss = [], []
for i, ev_batch_0 in enumerate(eval_dataloader_0):
    sp.call(f'echo EVAL BATCH BEFORE TRAINING {i+1}', shell=True)
    ev_input_ids_0, ev_attention_mask_0, ev_batch_labels_0 = ev_batch_0
    ev_input_ids_0 = ev_input_ids_0.to(device)
    ev_attention_mask_0 = ev_attention_mask_0.to(device)
    ev_batch_labels_0 = ev_batch_labels_0.to(device)
    with torch.no_grad():  # nor gradient computation
        # Predict the labels of eval_dataset (returns logits here)
        # where logits are the model's prediction without applying any
        # activation function (positive: more probable; negative: less probable)
        output_0 = model(ev_input_ids_0, attention_mask=ev_attention_mask_0, labels=ev_batch_labels_0)
        loss = loss_fn(output_0.logits.to(device), ev_batch_labels_0)
        total_loss += loss.item()


        # Convert logits to probabilities via the softmax activation function (in 2nd dimension of output)
        probabilities_0 = torch.softmax(output_0.logits, dim=1)

        # Get the predicted labels from these probabilities of each class
        pred_labels_0 = torch.argmax(probabilities_0, dim=1)  # get the index (0 or 1) of the highest probability
        ev_preds_0 += pred_labels_0.tolist()
        ev_labels_0 += ev_batch_labels_0.tolist()

# Compute F1 score
# The 'macro' makes it that each class (0 or 1) are taken into account equally (which is what we want)
fscore_0 = f1_score(ev_labels_0, ev_preds_0, average='macro')
sp.call(f'echo Fscore: {fscore_0}', shell=True)

# Save that f1_score before training
epoch_f_scores.append(fscore_0)

# Save loss before training
avg_loss = total_loss / len(eval_dataloader_0)  # across batches
epoch_loss.append(avg_loss)
sp.call(f"echo Epoch 0, Loss: {avg_loss}", shell=True)


# Save those epoch f1_scores and loss for each fold
f_df = pd.DataFrame(epoch_f_scores, columns=['f1_score_per_epoch'])
f_df.to_csv(output_f1, sep='\t', index=False)
l_df = pd.DataFrame(epoch_loss, columns=['loss_per_epoch'])
l_df.to_csv(output_loss, sep='\t', index=False)


