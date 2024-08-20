#!/usr/bin/python3
import sys
import GPUtil
import pandas as pd
import torch
import torch.nn as nn
import random
import numpy as np
import subprocess as sp
import sklearn
import time
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, accuracy_score, precision_recall_fscore_support
from torch.utils.data import TensorDataset, DataLoader
from torch.nn import BCEWithLogitsLoss, CrossEntropyLoss, MSELoss
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification, logging
from transformers.modeling_outputs import SequenceClassifierOutput
from typing import Optional, Union, Tuple
#logging.set_verbosity_error()

# Load params
pretrained_model = sys.argv[1]  # pretrained DNABert6 model
fold_num = str(sys.argv[2])  # training fold number
rs = int(sys.argv[3])  # random_state

# Load inputs
X_train = pd.read_csv(sys.argv[4], sep='\t')
y_train = pd.read_csv(sys.argv[5], sep='\t')
y_simple = y_train.drop(columns=['gene_id'])

# Convert sno labels so that expressed and pseudogene 
# snoRNAs are considered the same label (i.e. 1)
y_simple = y_simple.replace(2, 1)

# Get path of outputs
output_model = sys.argv[6]
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

# Define custom BertForSequenceClassification2 child class from the 
# parent BertForSequenceClassification class
class BertForSequenceClassification2(BertForSequenceClassification):
    def __init__(self, config, num_new_features):
        super().__init__(config)  # call parent class (BertForSequenceClassification)'s constructor
        # Add the number of new numerical features as input features of the classifier layer
        self.classifier = nn.Linear(config.hidden_size+num_new_features, config.num_labels)

    def forward(
        self,
        input_ids: Optional[torch.Tensor] = None,
        attention_mask: Optional[torch.Tensor] = None,
        token_type_ids: Optional[torch.Tensor] = None,
        position_ids: Optional[torch.Tensor] = None,
        head_mask: Optional[torch.Tensor] = None,
        inputs_embeds: Optional[torch.Tensor] = None,
        labels: Optional[torch.Tensor] = None,
        numerical_features: Optional[torch.Tensor] = None,
        output_attentions: Optional[bool] = None,
        output_hidden_states: Optional[bool] = None,
        return_dict: Optional[bool] = None,
    ) -> Union[Tuple[torch.Tensor], SequenceClassifierOutput]:
        r"""
        labels (`torch.LongTensor` of shape `(batch_size,)`, *optional*):
            Labels for computing the sequence classification/regression loss. Indices should be in `[0, ...,
            config.num_labels - 1]`. If `config.num_labels == 1` a regression loss is computed (Mean-Square loss), If
            `config.num_labels > 1` a classification loss is computed (Cross-Entropy).
        """
        return_dict = return_dict if return_dict is not None else self.config.use_return_dict

        outputs = self.bert(
            input_ids,
            attention_mask=attention_mask,
            token_type_ids=token_type_ids,
            position_ids=position_ids,
            head_mask=head_mask,
            inputs_embeds=inputs_embeds,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict,
        )

        # Output of pooling layer which comes from the last hidden transformer layer
        pooled_output = outputs[1]

        # Numerical features
        new_features = numerical_features

        # Concat the sequence features to the numerical features
        pooled_output = torch.cat((pooled_output, new_features), dim=1)

        # Add dropout to that layer
        pooled_output = self.dropout(pooled_output)
        logits = self.classifier(pooled_output)

        loss = None
        if labels is not None:
            if self.config.problem_type is None:
                if self.num_labels == 1:
                    self.config.problem_type = "regression"
                elif self.num_labels > 1 and (labels.dtype == torch.long or labels.dtype == torch.int):
                    self.config.problem_type = "single_label_classification"
                else:
                    self.config.problem_type = "multi_label_classification"

            if self.config.problem_type == "regression":
                loss_fct = MSELoss()
                if self.num_labels == 1:
                    loss = loss_fct(logits.squeeze(), labels.squeeze())
                else:
                    loss = loss_fct(logits, labels)
            elif self.config.problem_type == "single_label_classification":
                loss_fct = CrossEntropyLoss()
                loss = loss_fct(logits.view(-1, self.num_labels), labels.view(-1))
            elif self.config.problem_type == "multi_label_classification":
                loss_fct = BCEWithLogitsLoss()
                loss = loss_fct(logits, labels)
        if not return_dict:
            output = (logits,) + outputs[2:]
            return ((loss,) + output) if loss is not None else output

        return SequenceClassifierOutput(
            loss=loss,
            logits=logits,
            hidden_states=outputs.hidden_states,
            attentions=outputs.attentions,
        ) 


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

# Load pre-trained DNABERT model
tokenizer = AutoTokenizer.from_pretrained(pretrained_model)  # BertTokenizerFast
model = BertForSequenceClassification2.from_pretrained(pretrained_model, num_labels=2, 
        output_hidden_states=True, num_new_features=4)  # BertModel
model.to(device)
model.classifier.to(device)

# Set number of batches (per epoch) and epochs
num_epochs = 200
batch_size = 100  # nb of example per batch


# Define optimizer and loss function
optimizer = torch.optim.AdamW(model.parameters(), lr=2e-5)
loss_fn = torch.nn.CrossEntropyLoss(weight=torch.tensor([1/20, 1]).to(device))


# Train over given fold (fold_num) in stratified 10-fold CV
skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=rs)
fold_dict = {str(fold_index+1): [train_index, test_index]
            for fold_index, (train_index, test_index) in
            enumerate(skf.split(kmer_seqs, y_simple))}

train_index = fold_dict[fold_num][0]
test_index = fold_dict[fold_num][1]


# Standardize numerical features (separately to avoid data leakage between train and test sets)
scaler = StandardScaler()
num_cols = ['box_score', 'structure_mfe', 'terminal_stem_mfe', 'length']
train_num = X_train.loc[train_index][num_cols]
train_num = scaler.fit_transform(train_num)
train_num = torch.tensor(train_num).to(device)

test_num = X_train.loc[test_index][num_cols]
test_num = scaler.transform(test_num)
test_num = torch.tensor(test_num).to(device)


# Load train and eval (test) datasets (kmers of sequences)
x_train = [k for k in kmer_seqs if kmer_seqs.index(k) in train_index]
x_test = [k for k in kmer_seqs if kmer_seqs.index(k) in test_index]
Y_train = [y_simple.loc[i, 'target'] for i in train_index]
Y_test = [y_simple.loc[i, 'target'] for i in test_index]


# Load input sequences in right format (tokenize it for BERT)
inputs = tokenizer(x_train, return_tensors='pt').to(device)
tokens_dict = {i: example_tensor.to(torch.device("cpu")) for i, example_tensor in enumerate(inputs.input_ids.to(torch.device("cpu")))}
labels = torch.tensor(Y_train).to(device)
dataset = TensorDataset(inputs.input_ids, inputs.attention_mask, labels)
train_dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=False)

bt_time = time.time()
sp.call(f'echo Before training time: {bt_time}', shell=True)
# Iterate over epochs and batches per epoch
epoch_f_scores, epoch_loss = [], []
for epoch in range(num_epochs):
    p = f'EPOCH {epoch}'
    sp.call("echo " + p, shell=True)
    total_loss = 0.0

    for i, batch in enumerate(train_dataloader):
        sp.call(f'echo BATCH Train {i}', shell=True)
        input_ids, attention_mask, batch_labels = batch
        model.train()
        sp.call(f'echo After batch definition time: {time.time() - bt_time}', shell=True)
        t1 = time.time()
        input_ids, batch_labels = input_ids.to(torch.device("cpu")), batch_labels.to(device)
        attention_mask = attention_mask.to(device)
        sp.call(f'echo After to.(device) time: {time.time() - t1}', shell=True)
        t2 = time.time()

        # Get the index (in inputs) of chosen examples in batch i
        # if exact sequences are present multiple times, it's handled 
        # by the set (conserves only once the duplicated index)
        batch_indexes = list(set([i for i, tensor_ in tokens_dict.items() 
                            if any(torch.equal(tensor_.to(torch.device("cpu")), input_ids[j]) 
                            for j in range(input_ids.shape[0]))]))
        sp.call(f'echo After batch index time: {time.time() - t2}', shell=True)
        t3 = time.time()
        # Get the actual numerical features for the given indexes
        num_features_batch = train_num[batch_indexes, :].to(torch.float32).to(device)
        sp.call(f'echo After getting num features time: {time.time() - t3}', shell=True)
        t4 = time.time()

        input_ids = input_ids.to(device)
        
        # Forward pass
        outputs = model(input_ids, attention_mask=attention_mask, labels=batch_labels, 
                    numerical_features=num_features_batch)
        sp.call(f'echo After forward pass time: {time.time() - t4}', shell=True)
        t5 = time.time()
        #sp.call('nvidia-smi', shell=True)
        #sp.call(f"{GPUtil.showUtilization()}", shell=True)
        # Backpropagation
        loss = loss_fn(outputs.logits.to(device), batch_labels)
        total_loss += loss.item()

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        sp.call(f'echo After loss computation time: {time.time() - t5}', shell=True)

    avg_loss = total_loss / len(train_dataloader)  # across batches
    epoch_loss.append(avg_loss)
    sp.call(f"echo Epoch {epoch + 1}/{num_epochs}, Loss: {avg_loss}", shell=True)


    ## Compute f1_score on held-out fold
    model.eval()  # no more dropout

    # First, load input sequences in right format (tokenize it for BERT)
    eval_dataset = tokenizer(x_test, return_tensors='pt').to(device)
    eval_labels = torch.tensor(Y_test).to(device)
    eval_dataset = TensorDataset(eval_dataset.input_ids, eval_dataset.attention_mask, eval_labels)
    eval_dataloader = DataLoader(eval_dataset, batch_size=batch_size)
    n_eval_batch = len([i for i, b in enumerate(eval_dataloader)])

    # Evaluate model on eval dataset per batch for 1 epoch
    ev_preds, ev_labels = [], []
    for i, ev_batch in enumerate(eval_dataloader):
        sp.call(f'echo EVAL BATCH {i+1}', shell=True)
        ev_input_ids, ev_attention_mask, ev_batch_labels = ev_batch
        ev_input_ids = ev_input_ids.to(device)
        ev_attention_mask = ev_attention_mask.to(device)
        ev_batch_labels = ev_batch_labels.to(device)

        if i < n_eval_batch - 1:
            ev_new_features = test_num[i*batch_size:(i*batch_size+batch_size), :].to(torch.float32).to(device)
        else:  # for last batch (which can be smaller than batch_size)
            ev_new_features = test_num[i*batch_size:, :].to(torch.float32).to(device)

        with torch.no_grad():  # no gradient computation
            # Predict the labels of eval_dataset (returns logits here)
            # where logits are the model's prediction without applying any
            # activation function (positive: more probable; negative: less probable)
            output = model(ev_input_ids, attention_mask=ev_attention_mask, labels=ev_batch_labels, 
                    numerical_features=ev_new_features)

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

    # Save that f1_score for each epoch
    epoch_f_scores.append(fscore)
    if epoch == num_epochs - 1:
        sp.call('nvidia-smi', shell=True)


# Save those epoch f1_scores and loss for each fold
f_df = pd.DataFrame(epoch_f_scores, columns=['f1_score_per_epoch'])
f_df.to_csv(output_f1, sep='\t', index=False)
l_df = pd.DataFrame(epoch_loss, columns=['loss_per_epoch'])
l_df.to_csv(output_loss, sep='\t', index=False)




# Save model for that given fold (only save weights and parameters as it is lighter than saving the whole model)
model.to('cpu')
torch.save(model.state_dict(), output_model)


