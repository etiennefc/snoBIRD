#!/usr/bin/python3
import pandas as pd
import re

preds = snakemake.input.predictions_table
output = snakemake.output.predictions_tsv
test_set = pd.read_csv(snakemake.input.test_set, sep='\t')

positives_dict = {}
with open(preds, 'r') as f:
    i = 0
    for line in f:
        # Select only lines with C/D snoRNAs for which a rfam family was identified
        if ('#' and 'SNORA' not in line) & ('Small nucleolar RNA' in line):
            gene_id = re.split('RF[0-9]{5}   ', line)[1].split('-')[0].strip(' ')
            if gene_id not in positives_dict.keys():
                positives_dict[gene_id] = 'expressed_CD_snoRNA'

# Drop snoRNA pseudogene predictions
#test_set = test_set[test_set['gene_biotype'] != 'snoRNA_pseudogene']

# Create the infernal_rfam_prediction column
test_set['infernal_rfam_prediction'] = test_set['gene_id'].map(positives_dict)
test_set['infernal_rfam_prediction'] = test_set['infernal_rfam_prediction'].fillna('other')

test_set.to_csv(output, sep='\t', index=False)
