
#!/usr/bin/python3
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import optuna
from optuna import pruners
from optuna.trial import TrialState
import math
import random
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score

# Load tuning data and create tensor matrix (format used by pytorch)
#X = pd.read_csv(snakemake.input.X_tuning, sep='\t')
#x_tensor = torch.tensor(X.drop(columns=['gene_id']).values)
#y = pd.read_csv(snakemake.input.y_tuning, sep='\t')
#y_tensor = torch.tensor(y.drop(columns=['gene_id']).values)

"""TEMPPP TO DELETE"""
X = pd.read_csv(snakemake.input.X_tuning, sep='\t')
y = pd.read_csv(snakemake.input.y_tuning, sep='\t')
print(X)
print(y)
n, p, s = y[y['target'] == 0], y[y['target'] == 1], y[y['target'] == 2]
print(n,p,s)
n = n.sample(n=886, random_state=42)
p = p.sample(n=33, random_state=42)
s = s.sample(n=44, random_state=42)
y2 = pd.concat([n,p,s]).sample(frac=1, random_state=42)
ddd = y2.merge(X, on='gene_id', how='left')
X2 = ddd.drop(columns=['target'])
print(X2)
print(y2)
x_tensor = torch.tensor(X2.drop(columns=['gene_id']).values)
y_tensor = torch.tensor(y2.drop(columns=['gene_id']).values)

hyperparam_space = snakemake.params.hyperparams_search_space
output = snakemake.output.best_hyperparams


# Set reproducible randomness 
rs = int(snakemake.params.random_state)
torch.backends.cudnn.deterministic = True
random.seed(rs)
torch.manual_seed(rs)
torch.cuda.manual_seed(rs)
np.random.seed(rs)

# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")


# Define constants
input_size = len([col for col in X2.columns 
            if 'gene_id' not in col])  # number of input features (5 (ATCGN) nt * 211 of length + 4 intrinsic features)
output_size = len(pd.unique(y.target))  # number of class to predicts
total_length = len(X)  # i.e. nb of examples in input dataset
# 1 epoch corresponds to 1 forward pass and 1 backpropagation 
# on all examples (so there there is generally more than 1 epoch)
num_epochs = 5 
batch_size = 107  # nb of example per batch (this is an intermediate batch size)
num_batches = int(total_length / batch_size)  # the number of batches
num_trials = 300  # the number of tested hyperparameter combinations

# Define hyperparams space
def split_dict(d, key):
    # Return first and second value of list associated to key in dict d
    if len(d[key]) == 1:
        return d[key][0], d[key][0]
    else:
        return d[key][0], d[key][1]

min_layer, max_layer = split_dict(hyperparam_space, 'n_layers')
min_node, max_node = split_dict(hyperparam_space, 'n_nodes')
min_dropout, max_dropout = split_dict(hyperparam_space, 'dropout_rate')
min_l_rate, max_l_rate = split_dict(hyperparam_space, 'learning_rate')




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


# Define custom LSTM class
class LSTM_nn(nn.Module):
    # Define a bidirectional LSTM (processes the input sequence in both direction, which gives better context)
    # The activation function is tanh (hyperbolic tangent function, which returns value between -1 and 1)
    def __init__(self, input_size, hidden_sizes, num_layers, output_size, dropout_rate, 
                bidirectional=True):
        super(LSTM_nn, self).__init__()
        self.hidden_sizes = hidden_sizes # number of units/nodes in hidden layers
        self.num_layers = num_layers  # number of layers
        self.lstm_layers = nn.ModuleList()  # this is a list where we will append the hidden layers
        self.dropout = nn.Dropout(dropout_rate)  # save as a member variable so 
                                                 # that it will be used later in the forward pass

        # Append the different lstm layers
        for i in range(num_layers):
            if i == 0:  # first hidden layer (its input size is the number of input feature)
                self.lstm_layers.append(nn.LSTM(input_size, hidden_sizes[i], batch_first=True))  # batch_first allows for variable input_size
            else:  # all subsequent hidden layers (if num_layers >1)
                self.lstm_layers.append(nn.LSTM(hidden_sizes[i-1], hidden_sizes[i], batch_first=True))
        # Connects the last hidden layer to the output layer
        self.last_layer = nn.Linear(hidden_sizes[-1], output_size)  # this returns logits (will need to be converted to probabilities)

    # Forward pass (computing the output of a layer given a input tensor x)
    def forward(self, x):
        hiddens = []
        for i in range(self.num_layers):
            # By default, the initial hidden states (values in each node of each layer) is set to 0
            out, _ = self.lstm_layers[i](x)
            x = self.dropout(out)
            hiddens.append(out)

        # Return the output of the last layer (forward pass) only
        out = self.last_layer(hiddens[-1])  # no dropout in that last layer
        return out



def objective(trial):
    """ Objective function we want to mazimize in the case of """
    # Number of layers
    if max_layer == 1:
        num_layers = 1
    else:
        num_layers = trial.suggest_int("num_layers", min_layer, max_layer, step=1)
    # Number of nodes per layer (variable across layers)
    hidden_sizes = [trial.suggest_int(f"hidden_size_{i+1}", min_node, max_node) for i in range(num_layers)]
    # Dropout rate (constant across layers): probability of ignoring a neuron in a forward pass
    # This helps to prevent overfitting (the higher, the less chance of overfitting, but the more underfitting)
    dropout_rate = trial.suggest_float("dropout_rate", min_dropout, max_dropout)
    # Learning rate (used by the optimizer)
    learning_rate = trial.suggest_float("learning_rate", min_l_rate, max_l_rate, log=True)
    # Optimizer
    optimizer_name = trial.suggest_categorical("optimizer", hyperparam_space["optimizer"])
    
    # Instantiate model and loss function (compute this loss function equally with regards to the 3 classes)
    model = LSTM_nn(input_size, hidden_sizes, num_layers, output_size, dropout_rate).to(device)
    loss_fn = nn.CrossEntropyLoss(weight=torch.tensor([1/3, 1/3, 1/3])) 

    # Instantiate the optimizer
    if optimizer_name == "Adam":
        optimizer = optim.Adam(model.parameters(), lr=learning_rate, weight_decay=0.7)
    elif optimizer_name == "AdamW":
        optimizer = optim.AdamW(model.parameters(), lr=learning_rate, weight_decay=0.7)
    else:
        optimizer = optim.SGD(model.parameters(), lr=learning_rate, weight_decay=0.7)


    # Iterate over fold in stratified 3-fold CV
    skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=rs)
    f1_scores = []
    for fold_index, (train_index, test_index) in enumerate(skf.split(x_tensor.numpy(), y_tensor.numpy())):
        x_train, x_test = x_tensor[train_index], x_tensor[test_index]
        y_train, y_test = y_tensor[train_index], y_tensor[test_index]
        print(f'FOLD {fold_index + 1}')

        # Iterate over epochs 
        for epoch in range(num_epochs):  
            print(f'EPOCH {epoch + 1}')
            train_dataset = CustomDataset(x_train, y_train)
            train_dataset = torch.utils.data.DataLoader(train_dataset, 
                                batch_size=batch_size, shuffle=True)  # shuffle train_dataset between epochs
            #Iterate over batches comprised in 1 epoch
            for i, (input_batch, label_batch) in enumerate(train_dataset):
                # where input_batch: input samples with features in that batch
                # where label_batch: target labels of that batch to predict
                print(f'BATCH {i + 1}')
                if len(label_batch) != batch_size:  # to account for the last batch which will be <107
                    label_batch = label_batch.reshape(len(label_batch))  # reshape to 1d tensor
                else:
                    label_batch = label_batch.reshape(batch_size)  # reshape to 1d tensor

                # Enable training mode (activate dropout and gradient computation (for backprop))
                model.train()
                optimizer.zero_grad()  # set gradient values to 0 (gradient are 
                                    # computed in the backprop to update the model's params)
                output = model(input_batch.float())  # output is computed through forward pass
                loss = loss_fn(output, label_batch)  # loss is computed (comparison between predictions and true labels)
                loss.backward()  # backpropagation computation to update the model's params
                optimizer.step()  # update the model's params using the computed gradients (via optimizer)

                # Prune/stop trial if early stopping condition is met
                # i.e. if the trial's best intermediate result is worse 
                # than median of intermediate results of previous trials at the same step
                if trial.should_prune():
                    raise optuna.exceptions.TrialPruned()

        # Go into evaluation mode
        print(x_test.shape)
        eval_dataset = x_test.reshape(int(len(x_test)/batch_size), batch_size, input_size)
        eval_labels = y_test.reshape(len(y_test))  # reshape to 1d 
        model.eval()  # no more dropout 
        with torch.no_grad():  # nor gradient computation
            # Predict the labels of eval_dataset (returns logits here)
            # where logits are the model's prediction without applying any 
            # activation function (positive: more probable; negative: less probable)
            output = model(eval_dataset.float())
            # Convert logits to probabilities via the softmax activation function (in 2nd dimension of output)
            probabilities = torch.softmax(output, dim=2)
            # Get the predicted labels from these probabilities of each class  
            pred_labels = torch.argmax(probabilities, dim=2)  # get the index (0, 1 or 2) of the highest probability
            pred_labels = pred_labels.reshape(len(eval_labels))  # reshape to a 1d tensor of same length as eval_labels

            # Optimize for F1 score (so equally for both precision and recall)
            # The 'macro' makes it that each class (0, 1 or 2) are taken into account equally (which is what we want)
            fscore = f1_score(eval_labels.numpy(), pred_labels.numpy(), average='macro')
            f1_scores.append(fscore)
    
    # Compute the average f_score across the 3-fold CV
    avg_fscore = sum(f1_scores) / len(f1_scores)

    return avg_fscore





# Initialize study (to maximize f1_score using the TPE (Tree-structured Parzen Estimator: a 
# sequential-based optimization (SMBO) algorithm that uses bayesian optimization and 
# tree-structured search spaces) algorithm to sample hyperparameters)
study = optuna.create_study(direction="maximize", pruner=pruners.MedianPruner(), 
                            sampler=optuna.samplers.TPESampler(seed=rs)) 
study.optimize(objective, n_trials=num_trials)


# Return best trial only (best combination of hyperparams which maximizes f1_score)
best_trial = study.best_trial
print(f"Best trial: {best_trial.params}")
print(f"Best f1_score: {best_trial.value}")

# Create df out of the best hyperparams with the highest f1_score
best_hyperparams_dict = best_trial.params
best_hyperparams_dict['avg_f1_score_3fold_tuning'] = best_trial.value
params_df = pd.DataFrame(best_hyperparams_dict, index=[0])
params_df.to_csv(output, sep='\t', index=False)




