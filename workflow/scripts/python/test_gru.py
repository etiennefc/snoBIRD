#!/usr/bin/python3
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import random
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, accuracy_score, precision_recall_fscore_support
from figures import functions as ft
import matplotlib.pyplot as plt 
import seaborn as sns

# Outputs
df_metrics = snakemake.output.df_metrics_on_test
df_preds = snakemake.output.test_predictions

# Load training data and create tensor matrix (format used by pytorch)
X = pd.read_csv(snakemake.input.X_test, sep='\t')
x_tensor = torch.tensor(X.drop(columns=['gene_id']).values)
y = pd.read_csv(snakemake.input.y_test, sep='\t')
y_tensor = torch.tensor(y.drop(columns=['gene_id']).values)

# Define custom GRU_nn class
# Define a class to be able to separate dataset in batch
class CustomDataset(torch.utils.data.Dataset):
    def __init__(self, data, labels):
        self.data = data
        self.labels = labels

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        x = self.data[index]
        y = self.labels[index]
        return x, y


# Define custom GRU class
class GRU_nn(nn.Module):
    # Define a bidirectional GRU (processes the input sequence in both direction, which gives better context)
    # The activation function is tanh (hyperbolic tangent function, which returns value between -1 and 1)
    def __init__(self, input_size, hidden_sizes, num_layers, output_size, dropout_rate, 
                bidirectional=False, activation=torch.tanh):
        super(GRU_nn, self).__init__()
        self.hidden_sizes = hidden_sizes # number of units/nodes in hidden layers
        self.num_layers = num_layers  # number of layers
        self.gru_layers = nn.ModuleList()  # this is a list where we will append the hidden layers
        self.dropout = nn.Dropout(dropout_rate)  # save as a member variable so 
                                                 # that it will be used later in the forward pass

        # Append the different gru layers
        for i in range(num_layers):
            if i == 0:  # first hidden layer (its input size is the number of input feature)
                self.gru_layers.append(nn.GRU(input_size, hidden_sizes[i], batch_first=True))  # batch_first allows for variable input_size
            else:  # all subsequent hidden layers (if num_layers >1)
                self.gru_layers.append(nn.GRU(hidden_sizes[i-1], hidden_sizes[i], batch_first=True))
        # Connects the last hidden layer to the output layer
        self.last_layer = nn.Linear(hidden_sizes[-1], output_size)  # this returns logits (will need to be converted to probabilities)

    # Forward pass (computing the output of a layer given a input tensor x)
    def forward(self, x):
        hiddens = []
        for i in range(self.num_layers):
            # By default, the initial hidden states (values in each node of each layer) is set to 0
            out, _ = self.gru_layers[i](x)
            x = self.dropout(out)
            hiddens.append(out)

        # Return the output of the last layer (forward pass) only
        out = self.last_layer(hiddens[-1])  # no dropout in that last layer
        return out

# Load best model
input_size = len([col for col in X.columns 
            if 'gene_id' not in col])  # number of input features (5 (ATCGN) nt * 211 of length + 4 intrinsic features)
output_size = len(pd.unique(y.target))  # number of class to predict
best_hyperparams = pd.read_csv(snakemake.input.best_hyperparams, sep='\t')
best_hyperparams_dict = {k: best_hyperparams[k].values[0] for k in best_hyperparams.columns}
hidden_sizes = [v for k,v in best_hyperparams_dict.items() if 'hidden_size' in k]
num_layers = [v for k,v in best_hyperparams_dict.items() if 'num_layers' in k][0]
dropout_rate = [v for k,v in best_hyperparams_dict.items() if 'dropout_rate' in k][0]
model_path = [p for p in snakemake.input.all_models if 'fold_7.pt' in p][0]
model_gru = GRU_nn(input_size=input_size, hidden_sizes=hidden_sizes, num_layers=num_layers, 
                output_size=output_size, dropout_rate=dropout_rate)

model_gru.load_state_dict(torch.load(model_path)) 


# Test the model
eval_dataset = x_tensor.reshape(1, len(x_tensor), input_size)
eval_labels = y_tensor.reshape(len(y_tensor))  # reshape to 1d 
model_gru.eval()  # no more dropout 
with torch.no_grad():  # nor gradient computation
    # Predict the labels of eval_dataset (returns logits here)
    # where logits are the model's prediction without applying any 
    # activation function (positive: more probable; negative: less probable)
    output = model_gru(eval_dataset.float())
    # Convert logits to probabilities via the softmax activation function (in 2nd dimension of output)
    probabilities = torch.softmax(output, dim=2)
    # Get the predicted labels from these probabilities of each class  
    pred_labels = torch.argmax(probabilities, dim=2)  # get the index (0, 1 or 2) of the highest probability
    pred_labels = pred_labels.reshape(len(eval_labels))  # reshape to a 1d tensor of same length as eval_labels
    fscore = f1_score(eval_labels.numpy(), pred_labels.numpy(), average='macro')  # macro avg across 3 classes
    
    df = pd.DataFrame({'y_true': eval_labels, 'y_pred': pred_labels})
    # Merge with gene_id and save that prediction df
    df = pd.concat([X[['gene_id']], df], axis=1)
    df.to_csv(df_preds, sep='\t', index=False)

    # Compute accuracy and metrics 
    t_df = df[df.y_true != 1]
    sno_df = df[(df.y_true == 2) | (df.y_pred == 2)]
    sno_df = sno_df[sno_df['y_pred'] != 1]
    pseudosno_df = df[(df.y_true == 1) | (df.y_pred == 1)]
    #accuracy = accuracy_score(eval_labels.numpy(), pred_labels.numpy())  # all 3 classes combined
    accuracy = accuracy_score(t_df.y_true, t_df.y_pred)  # just sno and negatives
    # For class expressed_CD_snoRNA (2)
    TP_sno = len(sno_df[(sno_df.y_true == 2) & (sno_df.y_pred == 2)]) 
    FP_sno = len(sno_df[(sno_df.y_true != 2) & (sno_df.y_pred == 2)])
    FN_sno = len(sno_df[(sno_df.y_true == 2) & (sno_df.y_pred != 2)])
    if TP_sno + FP_sno == 0:
        precision_sno = 0
    else:
        precision_sno = TP_sno/(TP_sno + FP_sno)
    if TP_sno + FN_sno == 0:
        recall_sno = 0
    else:
        recall_sno = TP_sno/(TP_sno + FN_sno)
    fscore = 2*precision_sno*recall_sno /(precision_sno + recall_sno) 
    
    # For class snoRNA_pseudogene (1)
    TP_pseudosno = len(pseudosno_df[(pseudosno_df.y_true == 1) & (pseudosno_df.y_pred == 1)]) 
    FP_pseudosno = len(pseudosno_df[(pseudosno_df.y_true != 1) & (pseudosno_df.y_pred == 1)])
    FN_pseudosno = len(pseudosno_df[(pseudosno_df.y_true == 1) & (pseudosno_df.y_pred != 1)])
    if TP_pseudosno + FP_pseudosno == 0:
        precision_pseudosno = 0
    else:
        precision_pseudosno = TP_pseudosno/(TP_pseudosno + FP_pseudosno)
    if TP_pseudosno + FN_pseudosno == 0:
        recall_pseudosno = 0
    else:
        recall_pseudosno = TP_pseudosno/(TP_pseudosno + FN_pseudosno) 
    metrics_df = pd.DataFrame([[accuracy, fscore, precision_sno, 
                            recall_sno, precision_pseudosno, recall_pseudosno]], 
                            columns=['accuracy_3_classes', 'f1_score_3_classes',
                            'precision_sno', 'recall_sno', 'precision_pseudosno', 'recall_pseudosno'])
    metrics_df.to_csv(df_metrics, sep='\t', index=False)


