#!/usr/bin/python3
import pandas as pd
import collections as coll 
from sklearn.utils import shuffle

positives_paths = snakemake.input.positives
negatives_paths = snakemake.input.negatives
rs = snakemake.params.random_state
fixed_length = snakemake.wildcards.fixed_length
outputs = snakemake.output
short_name_dict = snakemake.params.short_name_dict

# Load dfs
box_score = pd.concat([pd.read_csv(p, sep='\t') 
            for p in snakemake.input.box_score]).reset_index(drop=True).filter(regex=r'gene_id|box_score')
structure_stability = pd.concat([pd.read_csv(p, sep='\t', names=['gene_id', 'structure_mfe']) 
            for p in snakemake.input.structure_stability 
            if '.tsv' in p])
terminal_stem = pd.concat([pd.read_csv(p, sep='\t', names=['gene_id', 'terminal_stem_mfe']) 
            for p in snakemake.input.terminal_stem_stability 
            if '.tsv' in p])
length = pd.concat([pd.read_csv(p, sep='\t', names=['gene_id', 'length'], header=0) 
            for p in snakemake.input.length])

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
    print(coll.Counter(df.gene_biotype))
    print(coll.Counter(df.species_name))

# Add gene_biotype column and keep only relevant columns for positives
pos = []
for df in [tuning_pos, training_pos, test_pos]:
    df['gene_biotype'] = 'expressed_CD_snoRNA'
    df = df[['gene_id', 'gene_biotype', 'species_name', 
            f'extended_{fixed_length}nt_sequence']]
    df['species_name'] = df['species_name'].map(short_name_dict)
    pos.append(df)

# Concat and shuffle respective sets; add also a target column 
# (what will be predicted)
# Merge intrinsic features 
added_features = box_score.merge(structure_stability, how='left', on='gene_id')
added_features = added_features.merge(terminal_stem, how='left', on='gene_id')
added_features = added_features.merge(length, how='left', on='gene_id')


for i, df in enumerate(pos):
    concat_df = pd.concat([df, neg[i]])
    concat_df = concat_df.merge(added_features, how='left', on='gene_id')
    concat_df['target'] = 'other'
    concat_df.loc[concat_df['gene_biotype'] == 'expressed_CD_snoRNA', 
                            'target'] = 'expressed_CD_snoRNA'
    concat_df.loc[concat_df['gene_biotype'] == 'snoRNA_pseudogene', 
                            'target'] = 'CD_snoRNA_pseudogene'
    concat_df = shuffle(shuffle(concat_df, random_state=rs), 
                        random_state=rs)
    concat_df.to_csv(outputs[i], sep='\t', index=False)
    print(coll.Counter(concat_df.species_name))
    print(coll.Counter(concat_df['target']))
