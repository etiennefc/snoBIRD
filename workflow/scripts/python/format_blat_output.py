#!/usr/bin/python3
import pandas as pd
import numpy as np
import itertools as it
from pybedtools import BedTool
import subprocess as sp

extension = snakemake.params.extension
cols = ['match', 'mismatch', 'rep_match', 'Ns', 'query_gap_count', 'query_gap_bases',
        'target_gap_count', 'target_gap_bases', 'strand', 'gene_name', 'query_size',
        'query_start', 'query_end', 'chr', 'target_size', 'start', 'end', 'block_count', 
        'block_sizes', 'query_starts', 'target_starts']
# where query is a snoRNA and target is a location in the species genome
df = pd.read_csv(snakemake.input.blat, sep='\t', skiprows=5, names=cols)

# Custom filter to remove snoRNAs not on real assembled chromosomes **TO UPDATE FOR OTHER SPECIES
df['chr'] = df['chr'].astype(str)
df = df[~df['chr'].str.contains('QNV|RZJ')]

# Remove alignments that have more than 5 nt gaps in the target sequence
df = df[df['target_gap_bases'] <= 5]

# Keep only the match with the highest matched nt per gene_name 
df = df.loc[df.groupby('gene_name').match.idxmax()]

# Remove duplicate snoRNAs with the same start/end (keep only the one with highest match)
df = df.loc[df.groupby(['chr', 'strand','start', 'end']).match.idxmax()]

df = df[['gene_name', 'chr', 'strand', 'start', 'end']]


# Merge snoRNAs that overlap and that are almost exactly at the same position (+- 3nt start and end)
sno_to_remove, sno_merged_to_add = [], []
for i, group in enumerate(df.groupby(['chr', 'strand'])):
    total_df = group[1].set_index('gene_name')
    d = total_df.to_dict()
    comb_start = [a for a in it.combinations(d['start'].values(), r=2)] # combinations of all pairs of starts
    comb_end = [a for a in it.combinations(d['end'].values(), r=2)] # combinations of all pairs of ends
    for i, start_pair in enumerate(comb_start):
        # find pairs of overlapping snoRNAs with similar starts/ends (+- 3nt)
        if (abs(start_pair[1] - start_pair[0]) <= 3) & (abs(comb_end[i][1] - comb_end[i][0]) <= 3): 
            if start_pair[1] == start_pair[0]:  # if the 2 snos have the same start
                temp_d = d.copy()
                sno1 = [i for i in d['start'] if d['start'][i]==start_pair[0]][0]  # sno1 name
                del temp_d['start'][[i for i in d['start'] if d['start'][i]==start_pair[0]][0]]
                sno2 = [i for i in temp_d['start'] if temp_d['start'][i]==start_pair[1]][0]  # sno2 name
                start1, end1 = temp_d['start'][sno2], temp_d['end'][sno1] 
                start2, end2 = temp_d['start'][sno2], temp_d['end'][sno2]
                start, end = start1, max([end1, end2]) # select the most downstream end
                sno_name_concat = f'{sno1}_{sno2}'
                sno_to_remove.append(sno1)
                sno_to_remove.append(sno2)
                merge_sno_row = [sno_name_concat, temp_d['chr'][sno1], temp_d['strand'][sno1], start, end]
                sno_merged_to_add.append(merge_sno_row)
            elif comb_end[i][1] == comb_end[i][0]:  # if the 2 snos have the same end
                temp_d = d.copy()
                sno1 = [i for i in d['end'] if d['end'][i]==comb_end[i][0]][0]  # sno1 name
                del temp_d['end'][[i for i in d['end'] if d['end'][i]==comb_end[i][0]][0]]
                sno2 = [i for i in temp_d['end'] if temp_d['end'][i]==comb_end[i][1]][0]  # sno2 name
                start1, end1 = temp_d['start'][sno1], temp_d['end'][sno2] 
                start2, end2 = temp_d['start'][sno2], temp_d['end'][sno2]
                start, end = min([start1, start2]), end2  # select the most upstream start
                sno_name_concat = f'{sno1}_{sno2}'
                sno_to_remove.append(sno1)
                sno_to_remove.append(sno2)
                merge_sno_row = [sno_name_concat, temp_d['chr'][sno1], temp_d['strand'][sno1], start, end]
                sno_merged_to_add.append(merge_sno_row)
            else: # if snos have similar starts and ends but not exactly   
                sno1 = [i for i in d['start'] if d['start'][i]==start_pair[0]][0]
                sno2 = [i for i in d['start'] if d['start'][i]==start_pair[1]][0]
                start1, end1 = d['start'][sno1], d['end'][sno1] 
                start2, end2 = d['start'][sno2], d['end'][sno2]
                start, end = min([start1, start2]), max([end1, end2]) # select the most upstream start and most downstream end
                sno_name_concat = f'{sno1}_{sno2}'
                sno_to_remove.append(sno1)
                sno_to_remove.append(sno2)
                merge_sno_row = [sno_name_concat, d['chr'][sno1], d['strand'][sno1], start, end]
                sno_merged_to_add.append(merge_sno_row)

# Remove rows that will be merged
df = df[~df['gene_name'].isin(sno_to_remove)]

# Concat merged_sno_df to df
merged_sno_df = pd.DataFrame(sno_merged_to_add, columns=['gene_name', 'chr', 'strand', 'start', 'end'])
final_df = pd.concat([df, merged_sno_df])

# Create fake bed
final_df['dot'], final_df['dot2'], final_df['source'], final_df['feature'], final_df['attr'] = '.', '.', 'article', 'gene', 'gene_id "test"'
sno_bed = final_df[['chr', 'start', 'end', 'gene_name', 'dot', 'strand', 'source', 'feature', 'dot2', 'attr']]
sno_bed.to_csv('temp_sno_bed_seq_'+snakemake.wildcards.sno_fasta, index=False, header=None, sep='\t')
bed = BedTool('temp_sno_bed_seq_'+snakemake.wildcards.sno_fasta)

# Update snoRNA sequence from the most recent genome 
fasta = bed.sequence(fi=snakemake.input.genome, nameOnly=True, s=True)
seq_dict = {}
with open(fasta.seqfn, 'r') as fasta_file:
    for line in fasta_file:
        if '>' in line:
            sno_name = line.strip('\n').replace('>', '').replace('(+)', '').replace('(-)', '')
        else:
            seq = line.strip('\n')
            seq_dict[sno_name] = seq


# Add gene_id column (gene_id will be for ex: SNORD45_M_mul_Zhang_2010)
species, article = str(snakemake.wildcards.sno_fasta).split('_cd_')
genus, species = species.split('_')
species_short = genus[0].capitalize() + '_' + species[0:3] # e.g. M_mul instead of macaca_mulatta
final_df['gene_id'] = final_df['gene_name'] + f'_{species_short}_{article}'

# Add genome_version, species name and classification, validation method and article source columns
final_df = final_df[['gene_id', 'gene_name', 'chr', 'strand', 'start', 'end']]
attributes_dict = snakemake.params.dataset_attribute
for k,v in attributes_dict.items():
    final_df[k] = v

# Add sequence col to final_df
final_df['sequence'] = final_df['gene_name'].map(seq_dict)

# Remove offset of 1 nt created by bedtools when converting to bed
final_df['start'] = final_df['start'] + 1

# Create an extended_sequence by adding 15 nt up- and downstream of each snoRNA
tempo_df = final_df.copy()
tempo_df[['dot1', 'dot2', 'source', 'feature']] = '.', '.', 'ensembl', 'gene'
tempo_df[['chr', 'start', 'end', 'gene_id', 'dot1', 'strand', 'source', 
        'feature', 'dot2', 'gene_name']].to_csv('temp_cd_'+snakemake.wildcards.sno_fasta+'.tsv', 
                                                sep='\t', index=False, header=False)
cd_bed = BedTool('temp_cd_'+snakemake.wildcards.sno_fasta+'.tsv')
extended_cd = cd_bed.slop(g=snakemake.input.chr_size, l=extension+1, r=extension)
sequences = extended_cd.sequence(fi=snakemake.input.genome, s=True, nameOnly=True)
extended_seq_dict = {}
with open(sequences.seqfn, 'r') as f:
    for line in f:
        if '>' in line:
            cd_id = line.replace('\n', '').replace('>', '').replace('(+)', '').replace('(-)', '')
        else:
            seq_ = line.replace('\n', '')
            extended_seq_dict[cd_id] = seq_
final_df['extended_sequence'] = final_df['gene_id'].map(extended_seq_dict)


sp.call('rm temp_sno_bed_seq_'+snakemake.wildcards.sno_fasta, shell=True)
sp.call('rm temp_cd_'+snakemake.wildcards.sno_fasta+'.tsv', shell=True)


final_df.to_csv(snakemake.output.df, sep='\t', index=False)
