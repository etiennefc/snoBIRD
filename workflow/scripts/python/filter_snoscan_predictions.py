#!/usr/bin/python3
import pandas as pd

fa = snakemake.input.predictions_fa
output = snakemake.output.predictions_tsv
test_set = pd.read_csv(snakemake.input.test_set, sep='\t')

# Extract which examples were predicted as snoRNAs by snoscan
predictions_id = {}
with open(fa, 'r') as f:
    for line in f:
        if line.startswith('>>'):
            gene_id = line.split(' ')[1]
            if gene_id not in predictions_id.keys():
                predictions_id[gene_id] = 'expressed_CD_snoRNA'

# Drop snoRNA pseudogene predictions
#test_set = test_set[test_set['gene_biotype'] != 'snoRNA_pseudogene']

# Create the snoscan_prediction column
test_set['snoscan_prediction'] = test_set['gene_id'].map(predictions_id)
test_set['snoscan_prediction'] = test_set['snoscan_prediction'].fillna('other')


test_set.to_csv(output, sep='\t', index=False)
