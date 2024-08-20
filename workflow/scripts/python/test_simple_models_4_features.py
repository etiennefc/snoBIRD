#!/usr/bin/python3
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import metrics
import pickle

""" Test each model performance on unseen test data and report their accuracy."""

X_test = pd.read_csv(snakemake.input.X_test, sep='\t', index_col='gene_id')
X_test = X_test[['box_score_norm', 'structure_mfe_norm', 'terminal_stem_mfe_norm', 'length_norm']]
y_test = pd.read_csv(snakemake.input.y_test, sep='\t')
y_test_copy = y_test.copy()
y_test = y_test.drop(columns=['gene_id'])

# Unpickle and thus instantiate the trained model defined by the 'models' wildcard
model = pickle.load(open(snakemake.input.pickled_trained_model, 'rb'))

# Predict label (other (0), snoRNA_pseudogene (1) or expressed_CD_snoRNA (2)) on test data and compare to y_test
y_pred = model.predict(X_test)
print(snakemake.wildcards.simple_models)
print(metrics.accuracy_score(y_test, y_pred))
acc = {}
acc[snakemake.wildcards.simple_models+'_test_accuracy'] = metrics.accuracy_score(y_test, y_pred)
acc_df = pd.DataFrame(acc, index=[0])
acc_df.to_csv(snakemake.output.test_accuracy, sep='\t', index=False)

# Save predictions to file
target_pred_col = str(snakemake.wildcards.simple_models)+'_prediction'
y_pred_df = pd.DataFrame(y_pred, columns=[target_pred_col])
preds = pd.concat([y_test_copy, y_pred_df], axis=1)
preds.to_csv(snakemake.output.y_preds, sep='\t', index=False)

# Calculate the performance on the snoRNA pseudogene class only
pseudosno_preds = preds[(preds['target'] == 1) | (preds[target_pred_col] == 1)]
tp = len(pseudosno_preds[(pseudosno_preds['target'] == 1) & (pseudosno_preds[target_pred_col] == 1)])
fp = len(pseudosno_preds[(pseudosno_preds['target'] != 1) & (pseudosno_preds[target_pred_col] == 1)])
fn = len(pseudosno_preds[(pseudosno_preds['target'] == 1) & (pseudosno_preds[target_pred_col] != 1)])
if (tp + fp) > 0:
    precision = tp / (tp + fp)
else:
    precision = 0
if (tp + fn) > 0:
    recall = tp / (tp + fn)
else:
    recall = 0
if (precision + recall) > 0:
    f1_score = (2 * precision * recall) / (precision + recall)
else:
    f1_score = 0

# Create df of these performance
score_df = pd.DataFrame([[precision, recall, f1_score]], 
                columns=['precision', 'recall', 'f1_score'])
score_df.to_csv(snakemake.output.pseudosno_performance, sep='\t', index=False)

#pseudosno_only = pseudosno_preds[pseudosno_preds['target'] == 1]
pseudosno_preds.to_csv(snakemake.output.pseudosno_preds, sep='\t', index=False)

