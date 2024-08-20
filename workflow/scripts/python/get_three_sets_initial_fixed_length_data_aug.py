#!/usr/bin/python3
import pandas as pd
from sklearn.utils import shuffle
import collections as coll

positives_paths = snakemake.input.positives
negatives_paths = snakemake.input.negatives
rs = snakemake.params.random_state
fixed_length = snakemake.wildcards.fixed_length
outputs = [snakemake.output.tuning, snakemake.output.training, 
            snakemake.output.test]
target_outputs = [snakemake.output.tuning_target, snakemake.output.training_target, 
                snakemake.output.test_target]
short_name_dict = snakemake.params.short_name_dict

# Load dfs
tuning_pos = pd.read_csv([p for p in positives_paths 
                        if 'tuning' in p][0], sep='\t')
training_pos = pd.read_csv([p for p in positives_paths 
                        if 'training' in p][0], sep='\t')
test_pos = pd.read_csv([p for p in positives_paths 
                        if 'test' in p][0], sep='\t')
tuning_neg = pd.read_csv([p for p in negatives_paths 
                        if 'tuning' in p][0], sep='\t')
training_neg = pd.read_csv([p for p in negatives_paths 
                        if 'training' in p][0], sep='\t')
test_neg = pd.read_csv([p for p in negatives_paths 
                        if 'test' in p][0], sep='\t')
neg = [tuning_neg, training_neg, test_neg] 


# Add gene_biotype column and keep only relevant columns for positives
pos = []
for df in [tuning_pos, training_pos, test_pos]:
    df = df[['gene_id', 'gene_biotype', 'species_name', 
            f'extended_{fixed_length}nt_sequence']]
    df['species_name'] = df['species_name'].map(short_name_dict)
    # For droso
    df.loc[df['gene_id'].str.startswith('FBg'), 'species_name'] = 'drosophila_melanogaster'
    pos.append(df)

# Concat and shuffle respective sets; add also a target column 
# (what will be predicted)
for i, df in enumerate(pos):
    concat_df = pd.concat([df, neg[i].drop(columns=['dataset'])])
    concat_df['target'] = 'other'
    concat_df.loc[concat_df['gene_biotype'] == 'expressed_CD_snoRNA', 
                            'target'] = 'expressed_CD_snoRNA'
    concat_df.loc[concat_df['gene_biotype'] == 'snoRNA_pseudogene', 'target'] = 'snoRNA_pseudogene'
    # Remove duplicates if any
    concat_df = concat_df.drop_duplicates(subset=f'extended_{fixed_length}nt_sequence').reset_index(drop=True)
    concat_df = shuffle(shuffle(concat_df, random_state=rs), 
                        random_state=rs)
    concat_df.to_csv(outputs[i], sep='\t', index=False)
    print(coll.Counter(concat_df.species_name))

    # Create and save target df
    target_df = concat_df[['gene_id', 'target']]
    target_df['target'] = target_df['target'].replace('other', 0)
    target_df['target'] = target_df['target'].replace('snoRNA_pseudogene', 1)
    target_df['target'] = target_df['target'].replace('expressed_CD_snoRNA', 2)
    target_df.to_csv(target_outputs[i], sep='\t', index=False)
