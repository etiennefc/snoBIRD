#!/usr/bin/python3
import pandas as pd
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
le = preprocessing.LabelEncoder()

short_dict = snakemake.params.short_name_dict
rs = snakemake.params.random_state
positives = pd.read_csv(snakemake.input.positives, sep='\t')

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
shuffle_sno['species'] = shuffle_sno['species_name'].str.replace('homo_sapiens', 'H_sapiens')
shuffle_sno['species'] = shuffle_sno['species'].str.replace('saccharomyces_cerevisiae', 'S_cerevisiae')
shuffle_sno['species'] = shuffle_sno['species'].str.replace('mus_musculus', 'M_musculus')
shuffle_sno['species'] = shuffle_sno['species'].map(short_dict)
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
for i, df in enumerate([ncRNA, haca, shuffle_sno, intronic, intergenic, exonic]):
    if i <= 1:
        df['extended_seq'] = df['extended_sequence']
    elif i == 2:
        df['extended_seq'] = df['shuffled_extended_seq']
    else:
        df['extended_seq'] = df['sequence']  # this is also an extended sequence


# Concat all negatives
all_dfs = []
for neg_df in [ncRNA, haca, shuffle_sno, intronic, intergenic, exonic]:
    df_ = neg_df[['gene_id', 'gene_biotype', 'species', 'extended_seq']]
    all_dfs.append(df_)

all_negatives = pd.concat(all_dfs)

# Drop duplicate negatives (only in intergenic regions that have multiple NNNNNN)
all_negatives = all_negatives.drop_duplicates(subset=['extended_seq'])

# Keep all shuffled snoRNA sequence (1:1 ratio compared to positives)
final_negatives = [all_negatives[all_negatives['gene_biotype'] == 'shuffled_expressed_CD_snoRNA']]

# Keep a ratio of 1:1 for H/ACA, tRNA, snRNA and premiRNA compared to positives (total 4:1)
for mncRNA in ['HACA_snoRNA', 'pre_miRNA', 'tRNA', 'snRNA']:
    mncRNA_df = all_negatives[all_negatives.gene_biotype == mncRNA]
    mncRNA_df = shuffle(shuffle(mncRNA_df, random_state=rs), random_state=rs)  # double shuffle
    selected_mncRNA_df = mncRNA_df.sample(n=len(positives), random_state=rs)
    final_negatives.append(selected_mncRNA_df)

# Keep a ratio of 3:1 for exonic regions and 6:1 for intronic and intergenic regions (total 15: 1)
for region in ['random_intronic_region', 'random_intergenic_region', 'random_exonic_region']:
    if region == 'random_exonic_region':  # less represented in the genome than introns and intergenic regions
        region_df = all_negatives[all_negatives.gene_biotype == region]
        region_df = shuffle(shuffle(region_df, random_state=rs), random_state=rs)  # double shuffle
        selected_region_df = region_df.sample(n = 3 * len(positives), random_state=rs)  # 3:1 ratio
        final_negatives.append(selected_region_df)
    else:
        region_df = all_negatives[all_negatives.gene_biotype == region]
        region_df = shuffle(shuffle(region_df, random_state=rs), random_state=rs)  # double shuffle
        selected_region_df = region_df.sample(n = 6 * len(positives), random_state=rs)  # 6:1 ratio
        final_negatives.append(selected_region_df)

# Concat all selected negatives
final_negatives_df = pd.concat(final_negatives)
final_negatives_df = final_negatives_df.rename(columns={
                            'extended_seq': 'extended_sequence',
                            'species': 'species_name'})


# Shuffle on last time and split the negatives into the tuning, training and test set
final_negatives_df = shuffle(shuffle(final_negatives_df, random_state=rs), random_state=rs)
le.fit(list(pd.unique(final_negatives_df.gene_biotype)))
y = le.transform(list(final_negatives_df.gene_biotype))  # this is to stratify according to the negative types

rest, tuning = train_test_split(final_negatives_df, train_size=0.9,  
                                                test_size=0.1, random_state=rs, stratify=y)
y2 = le.transform(list(final_negatives_df[final_negatives_df.gene_id.isin(rest.gene_id)].gene_biotype))
train, test = train_test_split(rest, train_size=0.78, # 78% of the remaining 90% gives ~ 70%
                            test_size=0.22, random_state=rs, stratify=y2) # 22% of the remaining 90% gives ~ 20%


tuning.to_csv(tuning_output, sep='\t', index=False)
train.to_csv(training_output, sep='\t', index=False)
test.to_csv(test_output, sep='\t', index=False)




