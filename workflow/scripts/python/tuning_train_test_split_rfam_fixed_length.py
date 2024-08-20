#!/usr/bin/python3
import pandas as pd
import collections as coll 
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split

seed = snakemake.params.random_seed
fixed_length = snakemake.wildcards.fixed_length  # fixed window length with 15 nt up/down stream of snoRNA

# Load dfs
extended_fixed_length_seq = pd.read_csv(snakemake.input.extended_sno_seq, sep='\t')
ext_dict = dict(zip(extended_fixed_length_seq.gene_id, extended_fixed_length_seq[f'extended_{fixed_length}nt_sequence']))
rfam_clans = pd.read_csv(snakemake.input.rfam_clans, sep='\t', 
                names=['rfam_clan_id', 'rfam_family_id'])
sno_literature = pd.read_csv(snakemake.input.sno_literature, sep='\t')
sno_literature = sno_literature[['gene_id', 'species_name']]

# Remove Drosophila C/D from Huang et al. 2005 (since they overlap with those from Sklias et al 2024)
sno_literature = sno_literature[sno_literature['species_name'] != 'D_melanogaster']

# Load Drosophila expressed C/D (in TGIRT-Seq Sklias et al. 2024)
droso = pd.read_csv(snakemake.input.sno_tgirt_droso[0], sep='\t')
droso = droso[['gene_id', 'species_name']]

tgirt_dfs = [droso]
for path in snakemake.input.sno_tgirt:
    sno_tgirt = pd.read_csv(path, sep='\t')
    sno_tgirt = sno_tgirt[['gene_id', 'species_name']]
    tgirt_dfs.append(sno_tgirt)
sno_species_df = pd.concat(tgirt_dfs + [sno_literature])

# Simplify species name and merge it to sno_rfam df
sno_species_df['species_name'] = sno_species_df['species_name'].str.replace('mus_', 'M_')
sno_species_df['species_name'] = sno_species_df['species_name'].str.replace('homo_', 'H_')
sno_species_df['species_name'] = sno_species_df['species_name'].str.replace('saccharomyces_', 'S_')
sno_rfam = pd.read_csv(snakemake.input.sno_rfam, sep='\t')
sno_rfam = sno_rfam.drop(columns=['species_name'])
sno_rfam = sno_rfam.merge(sno_species_df, how='left', on='gene_id')

# Create extended_{fixed_length}nt_sequence column
sno_rfam[f'extended_{fixed_length}nt_sequence'] = sno_rfam['gene_id'].map(ext_dict)

print(f'Initial nb of snoRNAs: {len(sno_rfam)}')
# Keep only 1 copy per snoRNA duplicates with exactly the same extended sequences
sno_rfam = sno_rfam.drop_duplicates(subset=[f'extended_{fixed_length}nt_sequence'])
print(f'After dropping duplicate extended sequences: {len(sno_rfam)}')





# Exclude all snoRNAs over 95th percentile because they 
# are too long and would affect the model's accuracy
sno_rfam['length'] = sno_rfam['end'] - sno_rfam['start'] + 1
sno_rfam = sno_rfam[sno_rfam.length <= int(fixed_length) - 30]  # -15 nt flanking on both sides of the snoRNA
sno_rfam = sno_rfam.drop(columns=['length'])
print(f'After <= {fixed_length}nt threshold: {len(sno_rfam)}')
print(f'Of {len(sno_rfam)} C/D, {len(sno_rfam[sno_rfam["rfam_family_id"].isna()])} without Rfam id, {len(sno_rfam[~sno_rfam["rfam_family_id"].isna()])} with Rfam id')

# For the snoRNAs with a Rfam id, separate between those with 1 vs 2 to 10 members vs >10 members per family
sno_w_rfam = shuffle(sno_rfam[~sno_rfam['rfam_family_id'].isna()], 
                    random_state=seed).reset_index(drop=True)
sno_families_1_member, sno_big_families, sno_families = {}, {}, {}
for i, group in enumerate(sno_w_rfam.groupby('rfam_family_id')):
    if len(group[1]) == 1:  # 1 member per family
        sno_families_1_member[group[0]] = group[1]
    elif len(group[1]) > 10:  # big families
        sno_big_families[group[0]] = group[1]
    else: # 2 to 10 members per family
        sno_families[group[0]] = group[1]

print(f'\nOf {len(sno_rfam[~sno_rfam["rfam_family_id"].isna()])} with Rfam id:')
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
                                random_state=seed).reset_index(drop=True)
                sno_id = random_rows.loc[j, 'gene_id']  # get sno_id at index j
                if sno_id not in sno_ids:  # if sno not already chosen
                    sno_row = list(random_rows.iloc[j, :].reset_index(drop=True)) # select that snoRNA
                    sno_rows.append(sno_row)
                    sno_ids.append(sno_id)
                    i += 1
        j += 1
    
    # Create list of 1 df from appended sno rows
    sno_rows_df_list = [pd.DataFrame((sno_rows), 
                        columns=['gene_id', 'chr', 'strand', 'start', 'end', 
                                'sequence', 'extended_sequence', 'rfam_family_id', 
                                'species_name', f'extended_{fixed_length}nt_sequence'])]
    # Complete with randomly chosen snoRNAs if 10 snoRNAs are not already chosen per family
    if remaining > 0:
        remaining_rows = temp_df[~temp_df['gene_id'].isin(sno_ids)].sample(n=remaining,
                                random_state=seed).reset_index(drop=True)
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
print(f'TOTAL OF C/D TO DISPATCH: {s1 + len(sno_rfam[sno_rfam["rfam_family_id"].isna()])} ({s1} + {len(sno_rfam[sno_rfam["rfam_family_id"].isna()])})\n')

# Merge the rfam clan id to the 1 and 2-10 snoRNA families df
# Concat the dfs 
sno_fam = pd.concat([df for id, df in remaining_sno_dict.items()])
sno_fam = sno_fam.merge(rfam_clans, how='left', on='rfam_family_id')

one_fam = pd.concat([df for id, df in sno_families_1_member.items()])
one_fam = one_fam.merge(rfam_clans, how='left', on='rfam_family_id')
sno_fam = pd.concat([sno_fam, one_fam])

print(f'Nb of C/D with Rfam Clan id in the 1_member and 2_10_members: {len(sno_fam[~sno_fam["rfam_clan_id"].isna()])}')
print(f'Choose pseudo-randomly ~{round(len(sno_fam[~sno_fam["rfam_clan_id"].isna()])*0.1)} C/D (10%) in the tuning ', 
        f'and ~{round(len(sno_fam[~sno_fam["rfam_clan_id"].isna()])* 0.2)} (20%) in the test set (by keeping the same clan in the same set)')

# First dispatch snoRNAs with the same Rfam clan id in the same set
sno_w_clan = shuffle(sno_fam[~sno_fam['rfam_clan_id'].isna()], random_state=seed)  # length = 326
big_clans_tuning, big_clans_test = [], []
used_clan_id = []
print(f'Based on :\n{coll.Counter(sno_w_clan["rfam_clan_id"])}')
tuning_occurences_per_clan, test_occurences_per_clan = [1, 2, 3, 4, 5, 7, 11], [2, 2, 3, 3, 4, 8, 9, 10, 12, 13]
print(f'Clans with the following number of C/D in it were chosen randomly for the\n-Tuning: {tuning_occurences_per_clan}--> sum={sum(tuning_occurences_per_clan)}')
print(f'-Test: {test_occurences_per_clan}--> sum={sum(test_occurences_per_clan)}')
print(f'-Train: the rest --> sum={len(sno_fam[~sno_fam["rfam_clan_id"].isna()]) - sum(tuning_occurences_per_clan) - sum(test_occurences_per_clan)}\n')

# The sum of occurences (33 and 65 for tuning and test) was chosen 
# pseudo-randomly to make sure that we have ~ 10% (33/326) and 20% (65/326) of sno in clans 
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
sno_no_or_1_rfam = sno_rfam[(sno_rfam['rfam_family_id'].isna()) | 
                            (sno_rfam['rfam_family_id'].isin(sno_families_1_member.keys()))]
sno_no_or_1_rfam = sno_no_or_1_rfam[~sno_no_or_1_rfam.rfam_family_id.isin(rfam_ids_chosen_clans)]

sno_no_or_1_rfam = shuffle(sno_no_or_1_rfam, random_state=seed).reset_index(drop=True) 

sno_no_or_1_rfam_rest, sno_no_or_1_rfam_tuning = train_test_split(sno_no_or_1_rfam, train_size=0.9,  
                                                test_size=0.1, random_state=seed)
sno_no_or_1_rfam_train, sno_no_or_1_rfam_test = train_test_split(sno_no_or_1_rfam_rest, train_size=0.78, # 78% of the remaining 90% gives ~ 70%
                                                test_size=0.22, random_state=seed) # 22% of the remaining 90% gives ~ 20%

print('Dispatch randomly across the 3 sets (10%, 70%, 20%) a total of')
print(f'{len(sno_no_or_1_rfam)} C/D without Rfam id ({len(sno_rfam[sno_rfam["rfam_family_id"].isna()])}) or '+
        f'with 1_member in Rfam family but no Clan id ({len(sno_no_or_1_rfam[~sno_no_or_1_rfam["rfam_family_id"].isna()])})\n')

# For snoRNAs with Rfam id and between 2 to 10 members per family (which are not in a rfam clan previously dispatched)
remaining_sno_nb = sum([len(v) for k, v in remaining_sno_dict.items() if k not in rfam_ids_chosen_clans])

print('\n""At this point, all C/D wo rfam id, with clan id or with 1 member per family have been dispatched""\n\n')
print(f'Remaining C/D to dispatch, i.e. 2_10_members (without Clan id):\n',
        f'{remaining_sno_nb} (Total 1_10_members - CD_w_clan - 1_member_w_family_wo_clan)')

# Count the remaining number of sno needed in each set
tuning_nb_remaining = round(0.1 * remaining_sno_nb)  # 31 C/D
train_nb_remaining = round(0.7 * remaining_sno_nb)  # 218 C/D
test_nb_remaining = round(0.2 * remaining_sno_nb)  # 62 C/D

print(f'10 % in the tuning: {tuning_nb_remaining}')
print(f'70 % in the training: {train_nb_remaining}')
print(f'20 % in the test: {test_nb_remaining}')

## Distribute the families pseudo-randomly with respect to the proportion of snoRNAs per set
fam_len_dict = {k:len(v) for k,v in remaining_sno_dict.items() if k not in rfam_ids_chosen_clans}
print('Based on:')
print(fam_len_dict)
tuning, train, test = [], [], []


# To get to the 31 remaining C/D to add in the tuning set,
# randomly choose a family of 10, 6, 5, 4, 3, 3 snoRNAs (31 snoRNAs)
# This combination of 10, 6, 5, 4, 3, 3 was manually picked to ensure a total of 31
def sample_cd(dictio, nb_per_family, given_list, big_df, n_samples, rs):
    """ Select from dictio all families of n (nb_per_family) snoRNAs, retrieve all 
        snoRNAs  of that family from big_df and pick randomly n_samples (i.e n different 
        families of size nb_per_family). Ensure that the family was not already selected 
        in given_list from another dataset"""
    selected_fams = [id for id, val in dictio.items() if (len(val) == nb_per_family) & (id not in given_list)]
    df_ = big_df[big_df['rfam_family_id'].isin(selected_fams)].drop_duplicates(subset=['rfam_family_id'])
    sno_selected = df_.sample(n=n_samples, random_state=rs)
    return sno_selected.rfam_family_id.tolist()

tun_fam_nb = [10, 6, 5, 4, 3]
tuning_occurences = [1, 1, 1, 1, 2]
tuning_ids = []
d = {k:v for k,v in remaining_sno_dict.items() if k not in rfam_ids_chosen_clans}
for i, number in enumerate(tun_fam_nb):
    ids = sample_cd(d, number, tuning_ids, sno_rfam, tuning_occurences[i], seed)
    tuning += ids
    for id in ids:
        tuning_ids.append(id)

filtered_sno = [df.gene_id.tolist() for id, df in d.items()]
filtered_sno = [item for sublist in filtered_sno for item in sublist]
filtered_df = sno_rfam[sno_rfam['gene_id'].isin(filtered_sno)]
tuning_df = filtered_df[filtered_df['rfam_family_id'].isin(tuning)]  # 31 C/D

# For test set, randomly choose a family of 10, 9, 8, 7, 6, 5, 4, 3, 3, 3, 2, 2 snoRNAs (62 snoRNAs)
# This combination of 10, 9, 8, 7, ... 2 was manually picked to ensure a total of 62
test_fam_nb = [10, 9, 8, 7, 6, 5, 4, 3, 2]
test_occurences = [1, 1, 1, 1, 1, 1, 1, 3, 2]
for i, number in enumerate(test_fam_nb):
    ids = sample_cd(d, number, tuning_ids, sno_rfam, test_occurences[i], seed)
    test += ids
    for id in ids:
        tuning_ids.append(id)

test_df = filtered_df[filtered_df['rfam_family_id'].isin(test)]  # 62 C/D

# For training set, select the remaining snoRNAs not in the test nor tuning sets
train_df = filtered_df[~filtered_df['rfam_family_id'].isin(
                        test+tuning+rfam_ids_chosen_clans+list(sno_no_or_1_rfam_tuning.rfam_family_id))]  # 201 C/D


print('Families with the following number of members were chosen randomly for the')
print(f'-Tuning: {tun_fam_nb} (2 fam of 3) --> sum={len(tuning_df)}')
print(f'-Tuning: {test_fam_nb} (3 fam of 3, 2 fam of 2) --> sum={len(test_df)}')
print(f'-Training: the rest --> sum={len(train_df)}')
# Concat the sets composed of families of 2-10 members to their respective set 
# composed of families with 0 or 1 rfam id
final_tuning = shuffle(pd.concat([tun_clan_df, sno_no_or_1_rfam_tuning, tuning_df]), 
                        random_state=seed).reset_index(drop=True)
final_train = shuffle(pd.concat([train_clan_df, sno_no_or_1_rfam_train, train_df]), 
                        random_state=seed).reset_index(drop=True)
final_test = shuffle(pd.concat([test_clan_df, sno_no_or_1_rfam_test, test_df]), 
                        random_state=seed).reset_index(drop=True)

final_tuning.to_csv(snakemake.output.tuning, sep='\t', index=False)
final_train.to_csv(snakemake.output.training, sep='\t', index=False)
final_test.to_csv(snakemake.output.test, sep='\t', index=False)


# Concat all in a df
final_tuning['dataset'] = 'Tuning'
final_train['dataset'] = 'Training'
final_test['dataset'] = 'Test'
final_all = pd.concat([final_tuning, final_train, final_test])
final_all.to_csv(snakemake.output.all_positives, sep='\t', index=False)

print(f'\n\nLength of positives in final tuning ({len(final_tuning)})')
print(f'Length of positives in final training ({len(final_train)})')
print(f'Length of positives in final test ({len(final_test)})')
print(f'FINAL TOTAL NUMBER OF C/D DISPATCHED: {len(final_all)}\n\n')
print('Proportion of species in tuning, training and test sets:\n')
print(coll.Counter(final_tuning.species_name), '\n')
print(coll.Counter(final_train.species_name), '\n')
print(coll.Counter(final_test.species_name))




