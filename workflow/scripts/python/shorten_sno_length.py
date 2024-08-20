#!/usr/bin/python3
import pandas as pd

# Load dfs
fixed_length = snakemake.wildcards.fixed_length
positives = pd.read_csv(snakemake.input.positives, sep='\t')
negatives = pd.read_csv(snakemake.input.negatives, sep='\t')
all_df = pd.concat([positives, negatives])
X_tun = pd.read_csv([p for p in snakemake.input.initial_sets if 'tuning_set' in p][0], sep='\t')
y_tun = pd.read_csv([p for p in snakemake.input.initial_sets if 'tuning_target' in p][0], sep='\t')
X_train = pd.read_csv([p for p in snakemake.input.initial_sets if 'training_set' in p][0], sep='\t')
y_train = pd.read_csv([p for p in snakemake.input.initial_sets if 'training_target' in p][0], sep='\t')
X_test = pd.read_csv([p for p in snakemake.input.initial_sets if 'test_set' in p][0], sep='\t')
y_test = pd.read_csv([p for p in snakemake.input.initial_sets if 'test_target' in p][0], sep='\t')

# Get predicted sno start and end in dict
predicted_start = dict(zip(all_df.gene_id, all_df.predicted_start))
predicted_end = dict(zip(all_df.gene_id, all_df.predicted_end))

def extract_substring(row):
    # Extract sno sequence based on predicted start/end in the window of fixed_length
    start, end = None, None
    gene_id = row.gene_id
    if gene_id in predicted_start.keys() and gene_id in predicted_end.keys():
        start = predicted_start[gene_id]
        end = predicted_end[gene_id]
    if start is not None and end is not None:
        return row[f'extended_{fixed_length}nt_sequence'][start-1:end]
    else:
        return row[f'extended_{fixed_length}nt_sequence']

# Extract sno sequence for X_tun, X_train and X_test
X_tun['seq2'] = X_tun.apply(extract_substring, axis=1)
X_tun = X_tun.rename(columns={f'extended_{fixed_length}nt_sequence': f'whole_extended_{fixed_length}nt_sequence', 
                    'seq2': f'extended_{fixed_length}nt_sequence'})

X_train['seq2'] = X_train.apply(extract_substring, axis=1)
X_train = X_train.rename(columns={f'extended_{fixed_length}nt_sequence': f'whole_extended_{fixed_length}nt_sequence', 
                    'seq2': f'extended_{fixed_length}nt_sequence'})

X_test['seq2'] = X_test.apply(extract_substring, axis=1)
X_test = X_test.rename(columns={f'extended_{fixed_length}nt_sequence': f'whole_extended_{fixed_length}nt_sequence', 
                    'seq2': f'extended_{fixed_length}nt_sequence'})

# Save dfs
X_tun.to_csv(snakemake.output.X_tuning, sep='\t', index=False)
y_tun.to_csv(snakemake.output.y_tuning, sep='\t', index=False)
X_train.to_csv(snakemake.output.X_train, sep='\t', index=False)
y_train.to_csv(snakemake.output.y_train, sep='\t', index=False)
X_test.to_csv(snakemake.output.X_test, sep='\t', index=False)
y_test.to_csv(snakemake.output.y_test, sep='\t', index=False)
