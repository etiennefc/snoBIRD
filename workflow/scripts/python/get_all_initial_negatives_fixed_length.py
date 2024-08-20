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
pseudosno = pd.read_csv(snakemake.input.human_snoRNA_pseudogenes, sep='\t')
pseudosno_mouse = pd.read_csv(snakemake.input.mouse_snoRNA_pseudogenes, sep='\t')
pseudosno_droso = pd.read_csv(snakemake.input.droso_snoRNA_pseudogenes, sep='\t')
rfam_pseudo = pd.read_csv(snakemake.input.rfam_pseudo, sep='\t')
rfam_clans = pd.read_csv(snakemake.input.rfam_clans, sep='\t', names=['rfam_clan_id', 'rfam_family_id'])

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
pseudosno['species'] = pseudosno['species_name']
pseudosno_mouse['species'] = pseudosno_mouse['species_name']
pseudosno_droso['species'] = pseudosno_droso['species_name']

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
for neg_df in [ncRNA, haca, shuffle_sno, intronic, intergenic, exonic, pseudosno, pseudosno_mouse, pseudosno_droso]:
    df_ = neg_df[['gene_id', 'gene_biotype', 'species', f'extended_{fixed_length}nt_sequence']]
    all_dfs.append(df_)

all_negatives = pd.concat(all_dfs)

# Drop duplicate negatives (only in intergenic regions that have multiple NNNNNN)
all_negatives = all_negatives.drop_duplicates(subset=[f'extended_{fixed_length}nt_sequence'])

# Keep all shuffled snoRNA sequence (1:1 ratio compared to positives)
final_negatives = [all_negatives[all_negatives['gene_biotype'].isin(
                        ['shuffled_expressed_CD_snoRNA'])]]

# Keep a ratio of 1:1 for H/ACA, tRNA, snRNA and premiRNA compared to positives (total 4:1)
for mncRNA in ['HACA_snoRNA', 'pre_miRNA', 'tRNA', 'snRNA']:
    mncRNA_df = all_negatives[all_negatives.gene_biotype == mncRNA]
    mncRNA_df = shuffle(shuffle(mncRNA_df, random_state=rs), random_state=rs)  # double shuffle
    selected_mncRNA_df = mncRNA_df.sample(n=len(positives), random_state=rs)
    final_negatives.append(selected_mncRNA_df)

# We tried 3:1 exon, 6:1 for intron and intergenic before, now try to keep
# a ratio of 50:1 for exonic, intronic and intergenic regions (total 150: 1, 10X more)
for region in ['random_intronic_region', 'random_intergenic_region', 'random_exonic_region']:
    region_df = all_negatives[all_negatives.gene_biotype == region]
    n_total = 50 * len(positives)
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


# Split snoRNA pseudogenes by keeping those in the same Rfam clan (first, and then Rfam family) in the same set
print(f'Initial nb of C/D snoRNA pseudogenes: {len(rfam_pseudo)}')
# Keep only 1 copy per snoRNA duplicates with exactly the same extended sequences
rfam_pseudo = rfam_pseudo.drop_duplicates(subset=[f'extended_{fixed_length}nt_sequence'])
print(f'After dropping duplicate extended sequences: {len(rfam_pseudo)}')

# Exclude all snoRNAs over 95th percentile because they 
# are too long and would affect the model's accuracy
rfam_pseudo['length'] = rfam_pseudo['end'] - rfam_pseudo['start'] + 1
rfam_pseudo = rfam_pseudo[rfam_pseudo.length <= int(fixed_length) - 30]  # -15 nt flanking on both sides of the snoRNA
rfam_pseudo = rfam_pseudo.drop(columns=['length'])
print(f'After <= {fixed_length}nt threshold: {len(rfam_pseudo)}')
print(f'Of {len(rfam_pseudo)} C/D, {len(rfam_pseudo[rfam_pseudo["rfam_family_id"].isna()])} without Rfam id, {len(rfam_pseudo[~rfam_pseudo["rfam_family_id"].isna()])} with Rfam id')

# For the snoRNAs with a Rfam id, separate between those with 1 vs 2 to 10 members vs >10 members per family
sno_w_rfam = shuffle(rfam_pseudo[~rfam_pseudo['rfam_family_id'].isna()], 
                    random_state=rs).reset_index(drop=True)
sno_families_1_member, sno_big_families, sno_families = {}, {}, {}
for i, group in enumerate(sno_w_rfam.groupby('rfam_family_id')):
    if len(group[1]) == 1:  # 1 member per family
        sno_families_1_member[group[0]] = group[1]
    elif len(group[1]) > 10:  # big families
        sno_big_families[group[0]] = group[1]
    else: # 2 to 10 members per family
        sno_families[group[0]] = group[1]

print(f'\nOf {len(rfam_pseudo[~rfam_pseudo["rfam_family_id"].isna()])} with Rfam id:')
print(sum([len(i) for id_, i in sno_families_1_member.items()]), '1_member')
print(sum([len(i) for id_, i in sno_families.items()]), '2_10_members')
print(sum([len(i) for id_, i in sno_big_families.items()]), '>10 members')

# For snoRNAs with > 10 members per family, reduce the number to 10 snoRNAs to avoid overrepresentation of some families
# (maximize the number of different species snoRNAs then randomly choose for the remaining snoRNAs)
def largest_factor(num):
    """ Returns the maximal number if times that num is contained within 10."""
    for i in range(10, 0, -1):  # go backward from num-1 to 1
        if i * num <= 10:
            return i * num
            
reduced_big_families = []
for id, df in sno_big_families.items():
    temp_df = df.copy()
    species_dict = dict(coll.Counter(df.species_name)) # nb sno per species (one family per dict)
    species_nb_per_family = len(species_dict.keys())
    largest = largest_factor(species_nb_per_family) # largest number of sno that can be equally chosen from all species
    remaining = 10 - largest  # remaining number of sno that will be randomly chosen within the rest
    i = 0  # i is the number of chosen sno (max(i) = largest)
    j = 0  # j is the row index in a species sno_df 
    sno_ids, sno_rows = [], []
    for n in range(1, largest+1):  # iterate until largest (to ensure that the maximum nb of sno is chosen)
        for name in species_dict.keys():  # iterate across the species present in that rfam sno family
            species_df = df[df['species_name'] == name]
            if (i < largest) & (j <= len(species_df) - 1):  # i < largest and row index must be present within given df 
                random_rows = species_df.sample(n=len(species_df),  # sample randomly all rows in species_df
                                random_state=rs).reset_index(drop=True)
                sno_id = random_rows.loc[j, 'gene_id']  # get sno_id at index j
                if sno_id not in sno_ids:  # if sno not already chosen
                    sno_row = list(random_rows.iloc[j, :].reset_index(drop=True)) # select that snoRNA
                    sno_rows.append(sno_row)
                    sno_ids.append(sno_id)
                    i += 1
        j += 1
    
    # Create list of 1 df from appended sno rows
    sno_rows_df_list = [pd.DataFrame((sno_rows), 
                        columns=['gene_id', 'gene_name', 'gene_biotype', 'species_name', 'chr', 'strand', 'start', 'end', 
                                'sequence', 'extended_sequence', 
                                f'extended_{fixed_length}nt_sequence', 'rfam_family_id'])]
    # Complete with randomly chosen snoRNAs if 10 snoRNAs are not already chosen per family
    if remaining > 0:
        remaining_rows = temp_df[~temp_df['gene_id'].isin(sno_ids)].sample(n=remaining,
                                random_state=rs).reset_index(drop=True)
        sno_rows_df_list.append(remaining_rows)
    reduced_big_families.append(sno_rows_df_list)

# Concat the rows per big family into one df per family, then add them to remaining sno_families dict
big_fam_final_df = [pd.concat(sub) for sub in reduced_big_families]
remaining_sno_dict = sno_families.copy()
for big_fam_df in big_fam_final_df:
    rfam_id = pd.unique(big_fam_df['rfam_family_id'])[0]
    remaining_sno_dict[rfam_id] = big_fam_df


print(f'After reducing the >10 members to 10: {sum([len(df) for df in big_fam_final_df])} instead of {sum([len(i) for id_, i in sno_big_families.items()])}')
print(f'So total of 2_10 members = {sum([len(df) for df in big_fam_final_df]) + sum([len(i) for id_, i in sno_families.items()])}', 
        f' ({sum([len(i) for id_, i in sno_families.items()])}+{sum([len(df) for df in big_fam_final_df])})')
print(f'So total of 1_10 members = {sum([len(df) for df in big_fam_final_df]) + sum([len(i) for id_, i in sno_families.items()]) + sum([len(i) for id_, i in sno_families_1_member.items()])}', 
        f'({sum([len(df) for df in big_fam_final_df]) + sum([len(i) for id_, i in sno_families.items()])} + {sum([len(i) for id_, i in sno_families_1_member.items()])})\n')
s1 = sum([len(df) for df in big_fam_final_df]) + sum([len(i) for id_, i in sno_families.items()]) + sum([len(i) for id_, i in sno_families_1_member.items()])
print(f'TOTAL OF C/D PSEUDOGENES TO DISPATCH: {s1 + len(rfam_pseudo[rfam_pseudo["rfam_family_id"].isna()])} ({s1} + {len(rfam_pseudo[rfam_pseudo["rfam_family_id"].isna()])})\n')

# Merge the rfam clan id to the 1 and 2-10 snoRNA families df
# Concat the dfs 
sno_fam = pd.concat([df for id, df in remaining_sno_dict.items()])
sno_fam = sno_fam.merge(rfam_clans, how='left', on='rfam_family_id')

one_fam = pd.concat([df for id, df in sno_families_1_member.items()])
one_fam = one_fam.merge(rfam_clans, how='left', on='rfam_family_id')
sno_fam = pd.concat([sno_fam, one_fam])

print(f'Nb of C/D pseudogenes with Rfam Clan id in the 1_member and 2_10_members: {len(sno_fam[~sno_fam["rfam_clan_id"].isna()])}')
print(f'Choose pseudo-randomly ~{round(len(sno_fam[~sno_fam["rfam_clan_id"].isna()])*0.1)} C/D pseudogenes (10%) in the tuning ', 
        f'and ~{round(len(sno_fam[~sno_fam["rfam_clan_id"].isna()])* 0.2)} (20%) in the test set (by keeping the same clan in the same set)')

# First dispatch snoRNAs with the same Rfam clan id in the same set
sno_w_clan = shuffle(sno_fam[~sno_fam['rfam_clan_id'].isna()], random_state=rs)  
big_clans_tuning, big_clans_test = [], []
used_clan_id = []
print(f'Based on :\n{coll.Counter(sno_w_clan["rfam_clan_id"])}')
tuning_occurences_per_clan, test_occurences_per_clan = [1, 1, 6], [1, 2, 4, 10]
print(f'Clans with the following number of C/D pseudogenes in it were chosen randomly for the\n-Tuning: {tuning_occurences_per_clan}--> sum={sum(tuning_occurences_per_clan)}')
print(f'-Test: {test_occurences_per_clan}--> sum={sum(test_occurences_per_clan)}')
print(f'-Train: the rest --> sum={len(sno_fam[~sno_fam["rfam_clan_id"].isna()]) - sum(tuning_occurences_per_clan) - sum(test_occurences_per_clan)}\n')

# The sum of occurences  was chosen 
# pseudo-randomly to make sure that we have ~10%  and 20%  of sno in clans 
# that are distributed proportionally to the % we want in the tuning, train, test sets (10, 70 and 20%)
tun_nb, test_nb = 0, 0
for occ in tuning_occurences_per_clan:
    for j, group in enumerate(sno_w_clan.groupby('rfam_clan_id', sort=False)):
        id_clan = str(group[0])
        if (id_clan not in used_clan_id) & (len(group[1]) == occ):
            used_clan_id.append(id_clan)
            big_clans_tuning.append(group[1])
            tun_nb += 1
            break

for occ in test_occurences_per_clan:
    for j, group in enumerate(sno_w_clan.groupby('rfam_clan_id', sort=False)):
        id_clan = str(group[0])
        if (id_clan not in used_clan_id) & (len(group[1]) == occ):
            used_clan_id.append(id_clan)
            big_clans_test.append(group[1])
            test_nb += 1
            break

tun_clan_df = pd.concat(big_clans_tuning)
test_clan_df = pd.concat(big_clans_test)
train_clan_df = sno_w_clan[~sno_w_clan.rfam_clan_id.isin(used_clan_id)]
rfam_ids_chosen_clans = list(pd.concat([tun_clan_df, test_clan_df, train_clan_df]).rfam_family_id)


# For the snoRNAs without a Rfam family or with a rfam id but only 1 member per family and no Rfam clan), 
# shuffle and split them randomly within the 3 sets (10 %, 70% and 20 % for tuning, 
# training and test set respectively)
sno_no_or_1_rfam = rfam_pseudo[(rfam_pseudo['rfam_family_id'].isna()) | 
                            (rfam_pseudo['rfam_family_id'].isin(sno_families_1_member.keys()))]
sno_no_or_1_rfam = sno_no_or_1_rfam[~sno_no_or_1_rfam.rfam_family_id.isin(rfam_ids_chosen_clans)]

sno_no_or_1_rfam = shuffle(sno_no_or_1_rfam, random_state=rs).reset_index(drop=True) 

sno_no_or_1_rfam_rest, sno_no_or_1_rfam_tuning = train_test_split(sno_no_or_1_rfam, train_size=0.9,  
                                                test_size=0.1, random_state=rs)
sno_no_or_1_rfam_train, sno_no_or_1_rfam_test = train_test_split(sno_no_or_1_rfam_rest, train_size=0.78, # 78% of the remaining 90% gives ~ 70%
                                                test_size=0.22, random_state=rs) # 22% of the remaining 90% gives ~ 20%

print('Dispatch randomly across the 3 sets (10%, 70%, 20%) a total of')
print(f'{len(sno_no_or_1_rfam)} C/D pseudogenes without Rfam id ({len(rfam_pseudo[rfam_pseudo["rfam_family_id"].isna()])}) or '+
        f'with 1_member in Rfam family but no Clan id ({len(sno_no_or_1_rfam[~sno_no_or_1_rfam["rfam_family_id"].isna()])})\n')

# For snoRNAs with Rfam id and between 2 to 10 members per family (which are not in a rfam clan previously dispatched)
remaining_sno_nb = sum([len(v) for k, v in remaining_sno_dict.items() if k not in rfam_ids_chosen_clans])

print('\n""At this point, all C/D pseudogenes wo rfam id, with clan id or with 1 member per family have been dispatched""\n\n')
print(f'Remaining C/D pseudogenes to dispatch, i.e. 2_10_members (without Clan id):\n',
        f'{remaining_sno_nb} (Total 1_10_members - CD_w_clan - 1_member_w_family_wo_clan)')

# Count the remaining number of sno needed in each set
tuning_nb_remaining = floor(0.1 * remaining_sno_nb)  
train_nb_remaining = round(0.7 * remaining_sno_nb)  
test_nb_remaining = round(0.2 * remaining_sno_nb)  

print(f'10 % in the tuning: {tuning_nb_remaining}')
print(f'70 % in the training: {train_nb_remaining}')
print(f'20 % in the test: {test_nb_remaining}')

## Distribute the families pseudo-randomly with respect to the proportion of snoRNAs per set
fam_len_dict = {k:len(v) for k,v in remaining_sno_dict.items() if k not in rfam_ids_chosen_clans}
print('Based on:')
print(fam_len_dict)
tuning, train, test = [], [], []


# To get to the remaining C/D pseudogenes to add in the tuning set,
# randomly choose a combination of different family sizes
# This combination was manually picked to ensure the right total
def sample_cd(dictio, nb_per_family, given_list, big_df, n_samples, rs):
    """ Select from dictio all families of n (nb_per_family) snoRNAs, retrieve all 
        snoRNAs  of that family from big_df and pick randomly n_samples (i.e n different 
        families of size nb_per_family). Ensure that the family was not already selected 
        in given_list from another dataset"""
    selected_fams = [id for id, val in dictio.items() if (len(val) == nb_per_family) & (id not in given_list)]
    df_ = big_df[big_df['rfam_family_id'].isin(selected_fams)].drop_duplicates(subset=['rfam_family_id'])
    sno_selected = df_.sample(n=n_samples, random_state=rs)
    return sno_selected.rfam_family_id.tolist()

tun_fam_nb = [2, 3]
tuning_occurences = [1, 1]
tuning_ids = []
d = {k:v for k,v in remaining_sno_dict.items() if k not in rfam_ids_chosen_clans}
for i, number in enumerate(tun_fam_nb):
    ids = sample_cd(d, number, tuning_ids, rfam_pseudo, tuning_occurences[i], rs)
    tuning += ids
    for id in ids:
        tuning_ids.append(id)

filtered_sno = [df.gene_id.tolist() for id, df in d.items()]
filtered_sno = [item for sublist in filtered_sno for item in sublist]
filtered_df = rfam_pseudo[rfam_pseudo['gene_id'].isin(filtered_sno)]
tuning_df = filtered_df[filtered_df['rfam_family_id'].isin(tuning)]  

# For test set, randomly choose a combination of different family sizes
# This combination of was manually picked to ensure the right total
test_fam_nb = [2, 6]
test_occurences = [3, 1]
for i, number in enumerate(test_fam_nb):
    ids = sample_cd(d, number, tuning_ids, rfam_pseudo, test_occurences[i], rs)
    test += ids
    for id in ids:
        tuning_ids.append(id)

test_df = filtered_df[filtered_df['rfam_family_id'].isin(test)]  # 62 C/D pseudogenes

# For training set, select the remaining snoRNAs not in the test nor tuning sets
train_df = filtered_df[~filtered_df['rfam_family_id'].isin(
                        test+tuning+rfam_ids_chosen_clans+list(sno_no_or_1_rfam_tuning.rfam_family_id))]  # 201 C/D pseudogenes


print('Families with the following number of members were chosen randomly for the')
print(f'-Tuning: {tun_fam_nb} (1 fam of 2 and 1 fam of 3) --> sum={len(tuning_df)}')
print(f'-Test: {test_fam_nb} (3 fam of 2, 1 fam of 6) --> sum={len(test_df)}')
print(f'-Training: the rest --> sum={len(train_df)}')
# Concat the sets composed of families of 2-10 members to their respective set 
# composed of families with 0 or 1 rfam id
final_tuning = shuffle(pd.concat([tun_clan_df, sno_no_or_1_rfam_tuning, tuning_df]), 
                        random_state=rs).reset_index(drop=True)
final_tuning = final_tuning[['gene_id', 'gene_biotype', 'species_name', f'extended_{fixed_length}nt_sequence']]
final_train = shuffle(pd.concat([train_clan_df, sno_no_or_1_rfam_train, train_df]), 
                        random_state=rs).reset_index(drop=True)
final_train = final_train[['gene_id', 'gene_biotype', 'species_name', f'extended_{fixed_length}nt_sequence']]
final_test = shuffle(pd.concat([test_clan_df, sno_no_or_1_rfam_test, test_df]), 
                        random_state=rs).reset_index(drop=True)
final_test = final_test[['gene_id', 'gene_biotype', 'species_name', f'extended_{fixed_length}nt_sequence']]




# Concat all in a df
final_all = pd.concat([final_tuning, final_train, final_test])

print(f'\n\nLength of C/D pseudogenes in final tuning ({len(final_tuning)})')
print(f'Length of C/D pseudogenes in final training ({len(final_train)})')
print(f'Length of C/D pseudogenes in final test ({len(final_test)})')
print(f'FINAL TOTAL NUMBER OF C/D pseudogenes DISPATCHED: {len(final_all)}\n\n')

# Concat with other negatives
final_tuning = pd.concat([final_tuning, tuning_neg])
final_train = pd.concat([final_train, train_neg])
final_test = pd.concat([final_test, test_neg])
final_tuning['dataset'] = 'Tuning'
final_train['dataset'] = 'Training'
final_test['dataset'] = 'Test'


final_tuning.to_csv(snakemake.output.tuning, sep='\t', index=False)
final_train.to_csv(snakemake.output.training, sep='\t', index=False)
final_test.to_csv(snakemake.output.test, sep='\t', index=False)
final_all = pd.concat([final_tuning, final_train, final_test])
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

