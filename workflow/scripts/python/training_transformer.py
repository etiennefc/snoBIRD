#!/usr/bin/python3
print('start')
import pandas as pd
import torch
import torch.nn as nn
import random
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, accuracy_score, precision_recall_fscore_support
from torch.utils.data import TensorDataset, DataLoader
from transformers import AutoTokenizer, BertForSequenceClassification,logging
#logging.set_verbosity_error()

print(torch.cuda.is_available())
pretrained_model = snakemake.params.pretrained_model
fold_num = snakemake.wildcards.fold_num
output_model = snakemake.output.model
output_loss = snakemake.output.fold_loss
output_f1 = snakemake.output.fold_f1_score

# Set reproducible randomness
rs = int(snakemake.params.random_state)
torch.backends.cudnn.deterministic = True
random.seed(rs)
torch.manual_seed(rs)
torch.cuda.manual_seed(rs)
np.random.seed(rs)
print('before cpu-gpu')
# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if use_cuda else "cpu")
print('after cpu-gpu')
# Load dfs
X_train = pd.read_csv(snakemake.input.X_train, sep='\t')  # to modify .[0]
y_train = pd.read_csv(snakemake.input.y_train, sep='\t')  # to modify .[0]
x_simple = X_train['extended_211nt_sequence']
y_simple = y_train.drop(columns=['gene_id'])


# Transform sequence of examples in training set into kmers (6-mers)
def seq2kmer(seq, k):
    """
    Convert original str sequence to kmers (str of all kmers)
    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers

train_seqs = list(X_train['extended_211nt_sequence'])
kmer_seqs = [seq2kmer(s, 6) for s in train_seqs]

print('Before model load')
# Load pre-trained DNABERT model
tokenizer = AutoTokenizer.from_pretrained(pretrained_model)  # BertTokenizerFast
model = BertForSequenceClassification.from_pretrained(pretrained_model, num_labels=3)  # BertModel
model = model.to(device)
print('PRE-TRAINED MODEL LOAD FINISHED')

# Set number of batches (per epoch) and epochs
num_epochs = 10
batch_size = 100  # nb of example per batch


# Define optimizer and loss function
optimizer = torch.optim.AdamW(model.parameters(), lr=2e-5)
loss_fn = torch.nn.CrossEntropyLoss(weight=torch.tensor([1/20, 1, 1]))


# Train over given fold (fold_num) in stratified 10-fold CV
skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=rs)
fold_dict = {str(fold_index+1): [train_index, test_index]
            for fold_index, (train_index, test_index) in
            enumerate(skf.split(kmer_seqs, y_simple))}

train_index = fold_dict[fold_num][0]
test_index = fold_dict[fold_num][1]


print(f'FOLD {fold_num}')

# Load train and eval (test) datasets
x_train = [k for k in kmer_seqs if kmer_seqs.index(k) in train_index]
x_test = [k for k in kmer_seqs if kmer_seqs.index(k) in test_index]
Y_train = [y_simple.loc[i, 'target'] for i in train_index]
Y_test = [y_simple.loc[i, 'target'] for i in test_index]


# Load input sequences in right format (tokenize it for BERT)
inputs = tokenizer(x_train, return_tensors='pt').to(device)
labels = torch.tensor(Y_train).to(device)
dataset = TensorDataset(inputs.input_ids, inputs.attention_mask, labels)
kwargs = {'num_workers': 1, 'pin_memory': True} if use_cuda else {}  # to speed up data load on GPU
train_dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True, **kwargs)

# Iterate over epochs and batches per epoch
epoch_f_scores, epoch_loss = [], []
for epoch in range(num_epochs):
    p = f'EPOCH {epoch}'
    print(p)
    total_loss = 0.0

    for i, batch in enumerate(train_dataloader):
        print(f'BATCH Train {i}')
        input_ids, attention_mask, batch_labels = batch
        model.train()
        outputs = model(input_ids, attention_mask=attention_mask, labels=batch_labels)
        loss = loss_fn(outputs.logits, batch_labels)
        total_loss += loss.item()

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    avg_loss = total_loss / len(train_dataloader)  # across batches
    epoch_loss.append(avg_loss)
    print(f"Epoch {epoch + 1}/{num_epochs}, Loss: {avg_loss}")


    ## Compute f1_score on held-out fold
    model.eval()  # no more dropout

    # First, load input sequences in right format (tokenize it for BERT)
    eval_dataset = tokenizer(x_test, return_tensors='pt')
    eval_labels = torch.tensor(Y_test)
    eval_dataset = TensorDataset(eval_dataset.input_ids, eval_dataset.attention_mask, eval_labels)
    eval_dataloader = DataLoader(eval_dataset, batch_size=500)

    # Evaluate model on eval dataset per batch for 1 epoch
    ev_preds, ev_labels = [], []
    for i, ev_batch in enumerate(eval_dataloader):
        print(f'EVAL BATCH {i+1}')
        ev_input_ids, ev_attention_mask, ev_batch_labels = ev_batch
        with torch.no_grad():  # nor gradient computation
            # Predict the labels of eval_dataset (returns logits here)
            # where logits are the model's prediction without applying any
            # activation function (positive: more probable; negative: less probable)
            output = model(ev_input_ids, attention_mask=ev_attention_mask, labels=ev_batch_labels)

            # Convert logits to probabilities via the softmax activation function (in 2nd dimension of output)
            probabilities = torch.softmax(output.logits, dim=1)

            # Get the predicted labels from these probabilities of each class
            pred_labels = torch.argmax(probabilities, dim=1)  # get the index (0, 1 or 2) of the highest probability
            ev_preds += pred_labels.tolist()
            ev_labels += ev_batch_labels.tolist()

    # Compute F1 score
    # The 'macro' makes it that each class (0, 1 or 2) are taken into account equally (which is what we want)
    fscore = f1_score(ev_labels, ev_preds, average='macro')
    print(f'Fscore: {fscore}')

    # Save that f1_score for each epoch
    epoch_f_scores.append(fscore)

# Save those epoch f1_scores and loss for each fold
f_df = pd.DataFrame(epoch_f_scores, columns=['f1_score_per_epoch'])
f_df.to_csv(output_f1, sep='\t', index=False)
l_df = pd.DataFrame(epoch_loss, columns=['loss_per_epoch'])
l_df.to_csv(output_loss, sep='\t', index=False)




# Save model for that given fold (only save weights and parameters as it is lighter than saving the whole model)
model.to('cpu')
torch.save(model.state_dict(), output_model)

