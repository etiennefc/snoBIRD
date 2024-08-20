#!/usr/bin/python3
import pandas as pd
import re

fa = snakemake.input.predictions_fa
output = snakemake.output.predictions_tsv
test_set = pd.read_csv(snakemake.input.test_set, sep='\t')

# Get the positive examples predicted by snoreport in the test set
# Keep in mind that some of these examples were used in snoreport2 training set
# so it might be an overly optimistic overview of snoreport performance
positives_dict = {}
with open(fa, 'r') as f:
    for line in f:
        if '>' in line:
            example_id = line.split(' ')[0].replace('>', '')
            if example_id.endswith('_1'):  # select only 1 prediction per example (drop those ending with _2)
                id = re.sub('_1$', '', example_id)
                positives_dict[id] = 'expressed_CD_snoRNA'

# Drop snoRNA pseudogene predictions
#test_set = test_set[test_set['gene_biotype'] != 'snoRNA_pseudogene']

# Create the snoreport2_prediction column
test_set['snoreport2_prediction'] = test_set['gene_id'].map(positives_dict)
test_set['snoreport2_prediction'] = test_set['snoreport2_prediction'].fillna('other')

test_set.to_csv(output, sep='\t', index=False)
