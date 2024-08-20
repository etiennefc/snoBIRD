#!/usr/bin/python3
import pandas as pd
import random
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, accuracy_score, precision_recall_fscore_support
from figures import functions as ft
import matplotlib.pyplot as plt 
import seaborn as sns
from pomegranate import *
import transformers
print(transformers.__version__)

# Load dfs
#pos = pd.read_csv(snakemake.input.positives, sep='\t')
#neg = pd.concat([pd.read_csv(p, sep='\t') for p in snakemake.input.negatives])
#pseudo = pd.read_csv(snakemake.input.pseudo, sep='\t')
X_train = pd.read_csv(snakemake.input.X_train[0], sep='\t')
X_test = pd.read_csv(snakemake.input.X_test[0], sep='\t')
y_train = pd.read_csv(snakemake.input.y_train[0], sep='\t')
#y_train['target'] = y_train['target'].replace({0:'other', 1:'CD_snoRNA_pseudogene', 2:'expressed_CD_snoRNA'})
y_test = pd.read_csv(snakemake.input.y_test[0], sep='\t')

# Select sequence of examples in training set
train_seqs = list(X_train['extended_211nt_sequence'])
encoded_seqs = [[['A', 'T',  'C', 'G', 'N'].index(nt) for nt in seq] for seq in train_seqs]

#dists = [NormalDistribution(5, 1), NormalDistribution(1, 7), NormalDistribution(8,2)]
#trans_mat = numpy.array([[0.7, 0.3, 0.0],
#[0.0, 0.8, 0.2],
#[0.0, 0.0, 0.9]])
#starts = numpy.array([1.0, 0.0, 0.0])
#ends = numpy.array([0.0, 0.0, 0.1])
#model = HiddenMarkovModel.from_matrix(trans_mat, dists, starts, ends)
#model.bake()
#print(np.array([[0.25, 0.25, 0.25, 0.25, 0], [0.25, 0.25, 0.25, 0.25, 0], [0.25, 0.25, 0.25, 0.25, 0]]))
#model.fit(encoded_seqs, algorithm='labeled', labels=[[i] for i in list(y_train.target)], verbose=True)
#model.fit(encoded_seqs[0],  algorithm='labeled', labels=(y_train.target[0]))



# HMM parameters
class_labels = ['other', 'CD_snoRNA_pseudogene', 'expressed_CD_snoRNA']
n_classes = 3  # nb of hidden_state (classes to predict) # 0, 1, 2 (other, pseudo and expressed CD)
transition_matrix = np.array([[0.8, 0.1, 0.1], [0.4, 0.3, 0.3], [0.4, 0.3, 0.3]])
emission_prob = {0: 0.25, 1: 0.25, 2: 0.25, 3: 0.25, 4: 0}  # A,T,C,G,N

#model = HiddenMarkovModel.from_samples(distribution=MultivariateGaussianDistribution, n_components=n_classes,
 #       X=encoded_seqs, labels=y_train, algorithm='labeled')

d1 = DiscreteDistribution({0 : 0.35, 1 : 0.20, 2 : 0.05, 3 : 0.40, 4: 0})
d2 = DiscreteDistribution({0 : 0.25, 1 : 0.25, 2 : 0.25, 3 : 0.25, 4: 0})
d3 = DiscreteDistribution({0 : 0.10, 1 : 0.40, 2 : 0.40, 3 : 0.10, 4: 0})
dd = [d1, d2, d3]






#
#model = HiddenMarkovModel(name="CD_predictor")
#
##print(sorted(model.__dir__()))
#
## Add states and their emission probabilities to model
#states = []
#for i in range(n_classes):
#    #dist = DiscreteDistribution(emission_prob)
#    #print(dist)
#    state = State(dd[i], name=class_labels[i])
#    model.add_state(state)
#    states.append(state)
#
## Add transition probabilities between states
#for i in range(n_classes):
#    for j in range(n_classes):
#        model.add_transition(states[i], states[j], transition_matrix[i][j])
#
## Create training dataset
#model.bake()
#print(model)
#print(len([[i] for i in list(y_train.target)]))
#print([len(i) for i in encoded_seqs[0:2]])
#print([[label]*211 for i, label in enumerate(list(y_train.target)) if i < 2])
#print(y_train.target[0:2])
#model.fit(encoded_seqs, algorithm='labeled', labels=[[label]*211 for i, label in enumerate(list(y_train.target))], verbose=True)
##model.fit(encoded_seqs[0:2], algorithm='labeled', verbose=True)
### weights= array of weights equal in length to the nb of training examples
#n = [1, 1, 1, 3, 3, 3, 0, 0, 1, 1, 3, 1, 0, 1, 3, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 3, 3, 0, 0, 3, 0, 1, 3, 3, 0, 2, 2, 0, 3, 3, 2, 0, 3, 1, 1, 3, 1, 3, 1, 1, 0, 3, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 2, 0, 2, 1, 1, 2, 0, 3, 0, 3, 3, 1, 1, 3, 1, 3, 0, 1, 3, 0, 0, 1, 1, 3, 2, 1, 1, 1, 0, 0, 1, 1, 1, 2, 1, 3, 0, 2, 0, 2, 2, 1, 2, 0, 1, 0, 1, 3, 0, 0, 0, 3, 2, 1, 3, 2, 0, 2, 1, 1, 3, 2, 0, 3, 1, 2, 1, 1, 3, 2, 0, 0, 3, 0, 2, 0, 3, 3, 2, 1, 1, 3, 3, 3, 3, 3, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 3, 0, 0, 3, 0, 0, 0, 0, 0, 0, 3, 0, 3, 0, 3, 0, 1, 0, 2, 1, 1, 3, 3, 3, 3, 1, 3, 3, 3, 0, 3, 3, 2, 0, 0, 3, 3, 2, 0, 3, 1, 1, 3, 3, 0, 0, 3, 3, 1, 1]
#
#log_probs = model.log_probability(n)
#predicted_class_index = np.argmax(log_probs)
#predicted_class_label = class_labels[predicted_class_index]
#print(predicted_class_label)
#
#
#
#

import torch
import torch.nn as nn

class SupervisedHMM(nn.Module):
    def __init__(self, n_hidden_states, n_observation_states):
        super(SupervisedHMM, self).__init__()
        self.n_hidden_states = n_hidden_states
        self.n_observation_states = n_observation_states

        # Define model parameters
        self.transition_probs = nn.Parameter(torch.randn(n_hidden_states, n_hidden_states))
        self.emission_probs = nn.Parameter(torch.randn(n_hidden_states, n_observation_states))
        self.initial_probs = nn.Parameter(torch.randn(n_hidden_states))

    def forward(self, observations, labels, sequence_length):
        # Implement the forward algorithm to compute the log-likelihood of the observations
        # based on the current model parameters
        #batch_size, sequence_length = observations.shape
        #batch_size = 233  # 13514 * 58
        #sequence_length = 211

        # Compute the log-likelihood of the observations given the model parameters
        log_likelihood = torch.zeros(batch_size)

        # Iterate over each sequence in the batch
        for i in range(batch_size):
            sequence = observations[i]
            label = labels[i]

            # Compute the initial probabilities based on the label
            initial_prob = self.initial_probs[label]
            print(sequence, label)
            print(initial_prob)

            # Compute the emission probabilities based on the label
            emission_prob = self.emission_probs[label]
            print(emission_prob)
            # Compute the transition probabilities based on the label
            transition_prob = self.transition_probs[label]
            print(transition_prob)
            # Forward Algorithm
            alpha = torch.zeros(211, self.n_hidden_states)
            alpha[0] = initial_prob * emission_prob[sequence[0]]

            for t in range(1, 211):
                alpha[t] = (alpha[t - 1] @ transition_prob) * emission_prob[sequence[t]]

            # Compute the log-likelihood of the sequence and add it to the batch log-likelihood
            log_likelihood[i] = alpha[-1].sum().log()

        return log_likelihood

    #def fit(self, observations, labels, n_iterations, learning_rate):
    #    optimizer = torch.optim.Adam(self.parameters(), lr=learning_rate)
    #    criterion = nn.NLLLoss()  # Negative Log-Likelihood Loss
#
    #    for iteration in range(n_iterations):
    #        print(str(iteration)+'ITTT')
    #        optimizer.zero_grad()
    #        log_likelihood = self.forward(observations, labels)
    #        loss = -log_likelihood.mean()
    #        print(str(loss)+'  LOSS')
    #        loss.backward()
    #        optimizer.step()

# Load your labeled sequences and class labels
# ...
# Encode the DNA sequences and class labels

# Create the Supervised HMM model
n_hidden_states = 3  # Number of hidden states (classes)
n_observation_states = 5  # Number of observation states (A, T, C, G, N)
num_batches = 58
batch_size = 233
hmm_model = SupervisedHMM(n_hidden_states, n_observation_states)

# Convert the sequences and labels to PyTorch tensors
observations = torch.tensor(encoded_seqs, dtype=torch.long)
print(observations)
print(len(observations[0]))
print(observations.shape)
print(y_train.target)
labels = torch.tensor(y_train.target, dtype=torch.long)

# Train the HMM on the labeled sequences
n_iterations = 5
learning_rate = 0.001
for epoch in range(n_iterations):
    total_loss = 0.0
    for batch_idx in range(num_batches):
        print(str(batch_idx)+'\n\n\n')
        start_idx = batch_idx * batch_size
        end_idx = start_idx + batch_size

        batch_observations = observations[start_idx:end_idx]
        batch_labels = labels[start_idx:end_idx]

        # Train the batch
        log_likelihood = hmm_model.forward(batch_observations, batch_labels, len(observations[0]))
        loss = -log_likelihood.mean()

        hmm_model.zero_grad()
        print(loss)
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    print(f"Epoch {epoch + 1}/{n_iterations}, Loss: {total_loss / num_batches}")



'''
# Load training data and create tensor matrix (format used by pytorch)
X = pd.read_csv(snakemake.input.X_train, sep='\t')
x_tensor = torch.tensor(X.drop(columns=['gene_id']).values)
y = pd.read_csv(snakemake.input.y_train, sep='\t')
y_tensor = torch.tensor(y.drop(columns=['gene_id']).values)

# Load best hyperparams
best_hyperparams_df = pd.read_csv(snakemake.input.best_hyperparams, sep='\t')
best_hyperparams_dict = {k: v[0] for k,v in 
                        best_hyperparams_df.to_dict(orient='list').items() 
                        if k != 'avg_f1_score_3fold_tuning'}

# Outputs
output_model = snakemake.output.trained_model
output_metrics = snakemake.output.training_metrics_per_fold
output_figure = snakemake.output.learning_curves

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
input_size = len([col for col in X.columns 
            if 'gene_id' not in col])  # number of input features (5 (ATCGN) nt * 211 of length + 4 intrinsic features)
output_size = len(pd.unique(y.target))  # number of class to predicts
total_length = len(X)  # i.e. nb of examples in input dataset
print(total_length)

num_epochs = 50 
batch_size = 107  # nb of example per batch (this is an intermediate batch size)
num_batches = int(total_length / batch_size)  # the number of batches
print(num_batches)


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



# Iterate over fold in stratified 10-fold CV
skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=rs)
fold_f1_scores, last_epoch_metrics_df, all_epochs_df = [], [], []
for fold_index, (train_index, test_index) in enumerate(skf.split(x_tensor.numpy(), y_tensor.numpy())):
    fold_i = str(fold_index + 1)
    print(f'FOLD {fold_i}')

    # Get best hyperparams separated
    hidden_sizes = [best_hyperparams_dict[k] for k in 
                    sorted(best_hyperparams_dict.keys()) if 'hidden_size_' in k]
    if 'num_layers' not in best_hyperparams_dict.keys():
        num_layers = len(hidden_sizes)
    else:
        num_layers = best_hyperparams_dict['num_layers']
    dropout_rate = best_hyperparams_dict['dropout_rate']
    learning_rate = best_hyperparams_dict['learning_rate']
    optimizer_name = best_hyperparams_dict['optimizer']

    # Initialize GRU and loss function
    model = GRU_nn(input_size, hidden_sizes, num_layers, output_size, dropout_rate)
    loss_fn = nn.CrossEntropyLoss(weight=torch.tensor([1/3, 1/3, 1/3]))

    # Initialize optimizer
    if optimizer_name == "Adam":
        optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    else:
        optimizer = optim.SGD(model.parameters(), lr=learning_rate)
    
    # Load dataset
    x_train, x_test = x_tensor[train_index], x_tensor[test_index]
    y_train, y_test = y_tensor[train_index], y_tensor[test_index]

    # Initialize early stopping

    ##implement EarlyStopping class!!!!!!!!!!!
    #early_stopping = EarlyStopping(patience=3, min_delta=0.01)  # early stops after 3 epochs where less than 0.01 f1 score improvement

    # Iterate over epoch
    epoch_f_scores = []
    for epoch in range(num_epochs):
        print(f'Epoch {epoch}')
        train_dataset = CustomDataset(x_train, y_train)
        train_dataset = torch.utils.data.DataLoader(train_dataset, 
                            batch_size=batch_size, shuffle=True)  # shuffle train_dataset between epochs
        
        #Iterate over batches comprised in 1 epoch
        for i, (input_batch, label_batch) in enumerate(train_dataset):
            # where input_batch: input samples with features in that batch
            # where label_batch: target labels of that batch to predict
            print(f'BATCH {i + 1}')
            print(label_batch.size())
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
        
        # Compute f1_score on held-out fold
        eval_dataset = x_test.reshape(1, len(x_test), input_size)
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

            # Save that f1_score for each epoch
            epoch_f_scores.append(fscore)

            # Save the metrics of the last epoch
            if epoch == (num_epochs - 1):
                df = pd.DataFrame({'y_true': eval_labels, 'y_pred': pred_labels})
                sno_df = df[(df.y_true == 2) | (df.y_pred == 2)]
                pseudosno_df = df[(df.y_true == 1) | (df.y_pred == 1)]
                accuracy = accuracy_score(eval_labels.numpy(), pred_labels.numpy())  # all 3 classes combined
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
                print([fold_i, accuracy, fscore, precision_sno, 
                                        recall_sno, precision_pseudosno, recall_pseudosno])
                metrics_df = pd.DataFrame([[fold_i, accuracy, fscore, precision_sno, 
                                        recall_sno, precision_pseudosno, recall_pseudosno]], 
                                        columns=['fold', 'accuracy_3_classes', 'f1_score_3_classes',
                                        'precision_sno', 'recall_sno', 'precision_pseudosno', 'recall_pseudosno'])
                last_epoch_metrics_df.append(metrics_df)
        
        # Implement early stopping
        #if early_stopping(fscore):  # if f1 score improved
            # Update the best model
         #   early_stopping.best_model_state = model.state_dict()
        #else:
        #    if early_stopping.early_stop:
        #        print("EARLY STOPPING!!")
        #        break



    # Save those epoch f1_scores for each fold
    fold_f1_scores.append(epoch_f_scores)

    # Plot the f1 score as a function of the number of epoch and save that graph for each fold
    output_path = [path for path in output_figure if f'_fold_{fold_i}.svg' in path][0]
    temp_df = pd.DataFrame({'epoch': [i for i in range(1, num_epochs+1)], 
                            'f1_score': epoch_f_scores})
    all_epochs_df.append(temp_df)
    ft.lineplot(temp_df, 'epoch', 'f1_score', None, 'Number of epoch', 'F1 score', 
                f'F1 score across epochs in fold #{fold_i}', 'grey', output_path)
    


    
    # Save model for that given fold (only save weights and parameters as it is lighter than saving the whole model)
    torch.save(model.state_dict(), [p for p in output_model if f'fold_{fold_i}.pt' in p][0])
    
# Concat f1_score across folds for all epochs
all_fold_epochs_df = pd.concat(all_epochs_df)
all_fold_epochs_df.to_csv(snakemake.output.all_fold_epochs_df, index=False, sep='\t')

# Save metrics df (metrics for the last epoch of each fold)
final_metrics_df = pd.concat(last_epoch_metrics_df)
avg = ['average_fold']
avg = avg + list(final_metrics_df[['accuracy_3_classes', 'f1_score_3_classes',
                            'precision_sno', 'recall_sno', 'precision_pseudosno', 
                            'recall_pseudosno']].mean())
final_metrics_df = pd.concat([final_metrics_df, pd.DataFrame([avg], 
                    columns=['fold', 'accuracy_3_classes', 'f1_score_3_classes',
                            'precision_sno', 'recall_sno', 'precision_pseudosno', 
                            'recall_pseudosno'])])    
final_metrics_df.to_csv(output_metrics, sep='\t', index=False)
'''

