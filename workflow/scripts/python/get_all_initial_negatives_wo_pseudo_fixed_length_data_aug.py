#!/usr/bin/python3
import pandas as pd
import collections as coll
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
le = preprocessing.LabelEncoder()
from math import floor

short_dict = snakemake.params.short_name_dict
rs = snakemake.params.random_state
positives = pd.read_csv(snakemake.input.positives, sep='\t')
fixed_length = snakemake.wildcards.fixed_length

tuning_output = snakemake.output.tuning
training_output = snakemake.output.training
test_output = snakemake.output.test

def concat_dfs(dfs_path_list):
    """ Concat vertically a list of df paths into one df"""
    dfs = []
    for path in dfs_path_list:
        df = pd.read_csv(path, sep='\t')
        if len(df) > 0:
            dfs.append(df)
    return pd.concat(dfs)

# Load all negative dfs
ncRNA = pd.read_csv(snakemake.input.ncRNA, sep='\t')
haca = concat_dfs(snakemake.input.haca)
shuffle_sno = pd.read_csv(snakemake.input.shuffle_sno, sep='\t')
intronic = concat_dfs(snakemake.input.random_intronic_regions)
intergenic = concat_dfs(snakemake.input.random_intergenic_regions)
exonic = concat_dfs(snakemake.input.random_exonic_regions)


# Filter out shuffled C/D sequences that are not present in the positives 
# (because we limited the number of sno per RFAM family to max 10)
shuffle_sno['id_modified'] = shuffle_sno['gene_id'].str.replace('_shuffle$', '', regex=True)
shuffle_sno = shuffle_sno[shuffle_sno.id_modified.isin(positives.gene_id)]
shuffle_sno = shuffle_sno.drop(columns=['id_modified'])

# Format species column
haca['species'] = haca['species_name'].map(short_dict)
shuffle_sno['species'] = shuffle_sno['species_name']
intronic['species'] = intronic['gene_id'].replace(r'_intronic_region_[0-9]*', '', regex=True)
intergenic['species'] = intergenic['gene_id'].replace(r'_intergenic_region_[0-9]*', '', regex=True)
exonic['species'] = exonic['gene_id'].replace(r'_exonic_region_[0-9]*', '', regex=True)


# Format biotype column
haca['gene_biotype'] = 'HACA_snoRNA'
shuffle_sno['gene_biotype'] = 'shuffled_expressed_CD_snoRNA'

# Format gene_id column
ncRNA['gene_id'] = ncRNA['rnacentral_id']
haca['gene_id'] = haca['gene_name']

# Format extended sequence column
for i, df in enumerate([shuffle_sno, intronic, intergenic, exonic]):
    if i == 0:
        df[f'extended_{fixed_length}nt_sequence'] = df[f'shuffled_extended_{fixed_length}nt_sequence']
    else:
        df[f'extended_{fixed_length}nt_sequence'] = df['sequence']  # this is also an extended sequence


# Concat all negatives
all_dfs = []
for neg_df in [ncRNA, haca, shuffle_sno, intronic, intergenic, exonic]:
    df_ = neg_df[['gene_id', 'gene_biotype', 'species', f'extended_{fixed_length}nt_sequence']]
    all_dfs.append(df_)

all_negatives = pd.concat(all_dfs)

# Drop duplicate negatives (only in intergenic regions that have multiple NNNNNN)
all_negatives = all_negatives.drop_duplicates(subset=[f'extended_{fixed_length}nt_sequence'])

# Keep all shuffled snoRNA sequence (1:1 ratio compared to positives)
final_negatives = [all_negatives[all_negatives['gene_biotype'].isin(
                        ['shuffled_expressed_CD_snoRNA'])]]

# Keep a ratio of ~1:1 for H/ACA, tRNA, snRNA and premiRNA compared to positives (total 4:1)
for mncRNA in ['HACA_snoRNA', 'pre_miRNA', 'tRNA', 'snRNA']:
    mncRNA_df = all_negatives[all_negatives.gene_biotype == mncRNA]
    mncRNA_df = shuffle(shuffle(mncRNA_df, random_state=rs), random_state=rs)  # double shuffle
    if len(mncRNA_df) > len(positives):
        selected_mncRNA_df = mncRNA_df.sample(n=len(positives), random_state=rs)
    else: # if not enough negatives to sample, keep all of them
        selected_mncRNA_df = mncRNA_df.copy()
    final_negatives.append(selected_mncRNA_df)

# We tried 3:1 exon, 6:1 for intron and intergenic before so 15:1 (and 150:1 total and 75:1 and 30:1 and 15:1 (5:1/5:1/5:1)), now try to keep
# a ratio of 1:1 for exonic, intronic and intergenic regions (total 3:1)
for region in ['random_intronic_region', 'random_intergenic_region', 'random_exonic_region']:
    region_df = all_negatives[all_negatives.gene_biotype == region]
    n_total = 1 * len(positives)
    # To ensure all species are at least well represented in the random regions negatives:
    # Take the species with the minimal nb of regions
    # Select that number for each species and then randomly select in the remaining regions to get to ~50:1 ratio
    all_spp = list(pd.unique(region_df.species))
    min_n = min([len(region_df[region_df['species'] == sp]) for sp in all_spp])
    selected_ids = []
    if min_n * len(all_spp) <= n_total:
        for sp in all_spp:
            sp_df = region_df[region_df['species'] == sp]
            sp_df = shuffle(shuffle(sp_df, random_state=rs), random_state=rs)  # double shuffle
            minimal_chosen_df = sp_df.sample(n = min_n, random_state=rs)
            final_negatives.append(minimal_chosen_df)
            selected_ids.extend(list(minimal_chosen_df.gene_id))
        remaining_region_df = region_df[~region_df['gene_id'].isin(selected_ids)]
        remaining_region_df = shuffle(shuffle(remaining_region_df, random_state=rs), random_state=rs)  # double shuffle
        selected_region_df = remaining_region_df.sample(n = n_total - len(selected_ids), random_state=rs)  # total of 50:1 ratio
        final_negatives.append(selected_region_df)
    else:  
    # if min_n is too big (too many regions (when multiplied by the number of all the other 
    # species) in the species with less region), limit to the max number of region per species
        max_n = n_total // len(all_spp)
        for sp in all_spp:
            sp_df = region_df[region_df['species'] == sp]
            sp_df = shuffle(shuffle(sp_df, random_state=rs), random_state=rs)  # double shuffle
            minimal_chosen_df = sp_df.sample(n = max_n, random_state=rs)
            final_negatives.append(minimal_chosen_df)
            selected_ids.extend(list(minimal_chosen_df.gene_id))
        remaining_region_df = region_df[~region_df['gene_id'].isin(selected_ids)]
        remaining_region_df = shuffle(shuffle(remaining_region_df, random_state=rs), random_state=rs)  # double shuffle
        selected_region_df = remaining_region_df.sample(n = n_total - len(selected_ids), random_state=rs)  # total of 50:1 ratio
        final_negatives.append(selected_region_df)




# Concat all selected negatives (except the snoRNA pseudogenes)
final_negatives_df = pd.concat(final_negatives)
final_negatives_df = final_negatives_df.rename(columns={
                            'species': 'species_name'})


# Shuffle one last time and split the negatives into the tuning, training and test set
final_negatives_df = shuffle(shuffle(final_negatives_df, random_state=rs), random_state=rs)

le.fit(list(pd.unique(final_negatives_df.gene_biotype)))
y = le.transform(list(final_negatives_df.gene_biotype))  # this is to stratify according to the negative types

rest, tuning_neg = train_test_split(final_negatives_df, train_size=0.9,  
                                                test_size=0.1, random_state=rs, stratify=y)
y2 = le.transform(list(final_negatives_df[final_negatives_df.gene_id.isin(rest.gene_id)].gene_biotype))
train_neg, test_neg = train_test_split(rest, train_size=0.78, # 78% of the remaining 90% gives ~ 70%
                            test_size=0.22, random_state=rs, stratify=y2) # 22% of the remaining 90% gives ~ 20%




# Finalize dfs
final_tuning = tuning_neg.copy()
final_train = train_neg.copy()
final_test = test_neg.copy()
final_tuning['dataset'] = 'Tuning'
final_train['dataset'] = 'Training'
final_test['dataset'] = 'Test'


final_tuning.to_csv(snakemake.output.tuning, sep='\t', index=False)
final_train.to_csv(snakemake.output.training, sep='\t', index=False)
final_test.to_csv(snakemake.output.test, sep='\t', index=False)
final_all = pd.concat([tuning_neg, train_neg, test_neg])
final_all.to_csv(snakemake.output.all_negatives, sep='\t', index=False)
print(f'Final length of negative tuning, training and test: {len(final_tuning)}, {len(final_train)}, {len(final_test)}\n')



print('\nComparison of gene_biotype proportion of the negatives in Tuning, training and test sets:')
print(coll.Counter(final_tuning.gene_biotype))
print(coll.Counter(final_train.gene_biotype))
print(coll.Counter(final_test.gene_biotype), '\n')

print('Comparison of species proportion of the negatives in Tuning, training and test sets:')
print(coll.Counter(final_tuning.species_name))
print(coll.Counter(final_train.species_name))
print(coll.Counter(final_test.species_name), '\n')

