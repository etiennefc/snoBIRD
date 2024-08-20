#!/usr/bin/python3
import pandas as pd
from sklearn.utils import shuffle

positives_paths = snakemake.input.positives
negatives_paths = snakemake.input.negatives
rs = snakemake.params.random_state
outputs = snakemake.output
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
for df in neg:
    import collections as coll
    print(coll.Counter(df.gene_biotype))
    print(coll.Counter(df.species_name))

# Add gene_biotype column and keep only relevant columns for positives
pos = []
for df in [tuning_pos, training_pos, test_pos]:
    df['gene_biotype'] = 'expressed_CD_snoRNA'
    df = df[['gene_id', 'gene_biotype', 'species_name', 
            'extended_sequence']]
    df['species_name'] = df['species_name'].map(short_name_dict)
    pos.append(df)

# Concat and shuffle respective sets; add also a target column 
# (what will be predicted)
for i, df in enumerate(pos):
    concat_df = pd.concat([df, neg[i]])
    concat_df['target'] = 'other'
    concat_df.loc[concat_df['gene_biotype'] == 'expressed_CD_snoRNA', 
                            'target'] = 'expressed_CD_snoRNA'
    concat_df = shuffle(shuffle(concat_df, random_state=rs), 
                        random_state=rs)
    concat_df.to_csv(outputs[i], sep='\t', index=False)
    import collections as coll 
    print(coll.Counter(concat_df.species_name))
