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
import optuna
from optuna.trial import TrialState
from optuna import pruners
#logging.set_verbosity_error()

# Load params
rs = int(float(sys.argv[1]))  # random_state
pretrained_model = sys.argv[2]  # pretrained DNABert6 model
fixed_length = sys.argv[3].split('nt.ts')[0].split('_')[-1]

# Load input tuning set
X_tuning = pd.read_csv(sys.argv[3], sep='\t')
y_tuning = pd.read_csv(sys.argv[4], sep='\t')

# Select only the C/D pseudogenes (0) and expressed C/D (1) for the tuning
X_tuning = X_tuning[X_tuning['target'] != 'other'].reset_index(drop=True)
y_tuning = y_tuning[y_tuning['gene_id'].isin(X_tuning.gene_id)]
y_simple = y_tuning.drop(columns=['gene_id']).reset_index(drop=True)
y_simple['target'] = y_simple['target'].replace(1, 0)
y_simple['target'] = y_simple['target'].replace(2, 1)

# Get path of outputs
output_best_hyperparams = sys.argv[5]

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

# Define hyperparams space
search_space = {"batch_size": [16, 32, 100], "learning_rate": [0.00002, 0.00003, 0.00004, 0.00005]}

# Set number of epochs
num_epochs = 30


# Transform sequence of examples in training set into kmers (6-mers)
def seq2kmer(seq, k):
    """
    Convert original str sequence to kmers (str of all kmers)
    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers

tuning_seqs = list(X_tuning[f'extended_{fixed_length}nt_sequence'])
kmer_seqs = [seq2kmer(s, 6) for s in tuning_seqs]

# Load tokenizer
tokenizer = AutoTokenizer.from_pretrained(pretrained_model)  # BertTokenizerFast


# Define optimizer loss function (weight equal for both expressed and pseudogene CD)
loss_fn = torch.nn.CrossEntropyLoss().to(device)


# Define objective function
def objective(trial):
    """ Objective function we want to maximize (f1-score) to get the best lr and batch size"""
    # Learning rate (used by the optimizer)
    learning_rate = trial.suggest_float("learning_rate", min(search_space['learning_rate']), max(search_space['learning_rate']))
    # Batch size
    batch_size = trial.suggest_int("batch_size", min(search_space["batch_size"]), max(search_space['batch_size']))


    # Create stratified CV strategy
    skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=rs)
    

    # Iterate over fold in stratified 3-fold CV
    f1_scores = []
    for fold_index, (train_index, test_index) in enumerate(skf.split(kmer_seqs, y_simple)):
        x_train = [k for k in kmer_seqs if kmer_seqs.index(k) in train_index]
        x_test = [k for k in kmer_seqs if kmer_seqs.index(k) in test_index]
        y_train = [y_simple.loc[i, 'target'] for i in train_index]
        y_test = [y_simple.loc[i, 'target'] for i in test_index]
        sp.call(f'echo FOLD {fold_index + 1}', shell=True)
        
        # Instantiate model for given fold (and also optimizer)
        model = BertForSequenceClassification.from_pretrained(pretrained_model, num_labels=2)  # BertModel
        model.to(device)
        model.classifier.to(device)
        optimizer = torch.optim.AdamW(model.parameters(), lr=learning_rate)

        # Iterate over epochs
        for epoch in range(num_epochs):
            sp.call(f'echo EPOCH {epoch + 1}', shell=True)
            # Load input sequences in right format (tokenize it for BERT)
            inputs = tokenizer(x_train, return_tensors='pt').to(device)
            labels = torch.tensor(y_train).to(device)
            dataset = TensorDataset(inputs.input_ids, inputs.attention_mask, labels)
            train_dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)


            #Iterate over batches comprised in 1 epoch
            for i, batch in enumerate(train_dataloader):
                sp.call(f'echo BATCH Train {i}', shell=True)
                input_ids, attention_mask, batch_labels = batch
                model.train()
                input_ids, batch_labels = input_ids.to(device), batch_labels.to(device)
                attention_mask = attention_mask.to(device)
                outputs = model(input_ids, attention_mask=attention_mask, labels=batch_labels)
                loss = loss_fn(outputs.logits.to(device), batch_labels)

                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

            # First, load input sequences in right format (tokenize it for BERT)
            eval_dataset = tokenizer(x_test, return_tensors='pt').to(device)
            eval_labels = torch.tensor(y_test).to(device)
            eval_dataset = TensorDataset(eval_dataset.input_ids, eval_dataset.attention_mask, eval_labels)
            eval_dataloader = DataLoader(eval_dataset, batch_size=batch_size)


            ## Compute f1_score on held-out fold
            model.eval()  # no more dropout

            # Evaluate model on eval dataset per batch for 1 epoch
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
        
                    # Convert logits to probabilities via the softmax activation function (in 2nd dimension of output)
                    probabilities = torch.softmax(output.logits, dim=1)
        
                    # Get the predicted labels from these probabilities of each class
                    pred_labels = torch.argmax(probabilities, dim=1)  # get the index (0 or 1) of the highest probability
                    ev_preds += pred_labels.tolist()
                    ev_labels += ev_batch_labels.tolist()
        
            # Compute F1 score
            # The 'macro' makes it that each class (0 or 1) are taken into account equally (which is what we want)
            fscore = f1_score(ev_labels, ev_preds, average='macro')
            sp.call(f'echo Fscore: {fscore}', shell=True)

            # Save fscore for last epoch of training of given fold
            if epoch + 1 == num_epochs:
                f1_scores.append(fscore)
        
            # Save that f1_score for each epoch
            if epoch == num_epochs - 1:
                sp.call('nvidia-smi', shell=True)


    # Compute the average f_score across the 3-fold CV
    avg_fscore = sum(f1_scores) / len(f1_scores)

    return avg_fscore



# Initialize study (to maximize f1_score using Grid Search 
study = optuna.create_study(direction="maximize", 
                            sampler=optuna.samplers.GridSampler(search_space))
study.optimize(objective)


# Return best trial only (best combination of hyperparams which maximizes f1_score)
best_trial = study.best_trial
sp.call(f"echo Best trial: {best_trial.params}", shell=True)
sp.call(f"echo Best f1_score: {best_trial.value}", shell=True)

# Create df out of the best hyperparams with the highest f1_score
best_hyperparams_dict = best_trial.params
best_hyperparams_dict['avg_f1_score_3fold_tuning'] = best_trial.value
params_df = pd.DataFrame(best_hyperparams_dict, index=[0])
params_df.to_csv(output_best_hyperparams, sep='\t', index=False)


