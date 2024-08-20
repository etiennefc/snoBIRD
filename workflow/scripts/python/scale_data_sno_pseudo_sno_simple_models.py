#!/usr/bin/python3
import pandas as pd
from sklearn.preprocessing import StandardScaler
tune = pd.read_csv(snakemake.input.X_tuning, sep='\t')
train = pd.read_csv(snakemake.input.X_train, sep='\t')
test = pd.read_csv(snakemake.input.X_test, sep='\t')
cols = ['box_score', 'structure_mfe', 'terminal_stem_mfe', 'length']

# Keep only the expressed C/D and C/D pseudogenes
tune = tune[tune['target'] != 'other']
train = train[train['target'] != 'other']
test = test[test['target'] != 'other']

# Create and fit the scaler on the training data only
scaler = StandardScaler()
scaler.fit(train[cols])

# Apply scaler to all the sets separately
tune_scaled = pd.concat([tune[['gene_id']].reset_index(drop=True), pd.DataFrame(scaler.transform(tune[cols]), columns=cols).reset_index(drop=True)], axis=1)
train_scaled = pd.concat([train[['gene_id']].reset_index(drop=True), pd.DataFrame(scaler.transform(train[cols]), columns=cols).reset_index(drop=True)], axis=1)
test_scaled = pd.concat([test[['gene_id']].reset_index(drop=True), pd.DataFrame(scaler.transform(test[cols]), columns=cols).reset_index(drop=True)], axis=1)

# Save dfs
tune_scaled.to_csv(snakemake.output.X_tuning_scaled, index=False, sep='\t')
train_scaled.to_csv(snakemake.output.X_train_scaled, index=False, sep='\t')
test_scaled.to_csv(snakemake.output.X_test_scaled, index=False, sep='\t')
