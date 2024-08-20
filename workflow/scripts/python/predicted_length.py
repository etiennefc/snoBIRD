#!/usr/bin/python3
import pandas as pd
import re 

output_pos = snakemake.output.positives_df
output_neg = snakemake.output.negatives_df

# Get the length of the predicted snoRNA
pos_len, neg_len = {}, {}
with open(snakemake.input.positives_fa, 'r') as f:
    for line in f:
        if line.startswith('>'):
            id_ = line.replace('>', '').replace('\n', '')
        elif re.match(r"^A|T|C|G|U|N", line):
            len_ = len(line.replace('\n', ''))
            pos_len[id_] = len_

with open(snakemake.input.negatives_fa, 'r') as f:
    for line in f:
        if line.startswith('>'):
            id_ = line.replace('>', '').replace('\n', '')
        elif re.match(r"^A|T|C|G|U|N", line):
            len_ = len(line.replace('\n', ''))
            neg_len[id_] = len_

pos_df = pd.DataFrame(list(pos_len.items()), columns=['gene_id', 'predicted_length'])
neg_df = pd.DataFrame(list(neg_len.items()), columns=['gene_id', 'predicted_length'])

pos_df.to_csv(output_pos, sep='\t', index=False)
neg_df.to_csv(output_neg, sep='\t', index=False)
