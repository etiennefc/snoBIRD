#!/usr/bin/python3
import pandas as pd
import numpy as np
import itertools as it
from pybedtools import BedTool
import subprocess as sp
import re

cols = ['match', 'mismatch', 'rep_match', 'Ns', 'query_gap_count', 'query_gap_bases',
        'target_gap_count', 'target_gap_bases', 'strand', 'gene_name', 'query_size',
        'query_start', 'query_end', 'chr', 'target_size', 'start', 'end', 'block_count', 
        'block_sizes', 'query_starts', 'target_starts']

fixed_length = int(snakemake.wildcards.fixed_length)
chr_size = snakemake.input.chr_size 
genome = snakemake.input.genome

# where query is a snoRNA and target is a location in the species genome
df = pd.read_csv(snakemake.input.blat, sep='\t', skiprows=5, names=cols)

# Custom filter to remove snoRNAs not on real assembled chromosomes 
df['chr'] = df['chr'].astype(str)
df = df[~df['chr'].str.contains('QNV|RZJ|GL|JH5|ML1|KI|MU1')]

# Remove alignments that have more than 5 nt gaps in the target sequence
df = df[df['target_gap_bases'] <= 5]

# Keep only the match with the highest matched nt per gene_name 
df = df.loc[df.groupby('gene_name').match.idxmax()]

# Remove duplicate snoRNAs with the same start/end (keep only the one with highest match)
df = df.loc[df.groupby(['chr', 'strand','start', 'end']).match.idxmax()]
df['gene_id'] = df['gene_name']
df = df[['gene_id', 'gene_name', 'chr', 'strand', 'start', 'end']]

# Merge snoRNAs that overlap 
final_df = df.copy()
final_df['dot'], final_df['dot2'], final_df['source'], final_df['feature'], final_df['attr'] = '.', '.', 'article', 'gene', 'gene_id "test"'
sno_bed = final_df[['chr', 'start', 'end', 'gene_name', 'dot', 'strand', 'source', 'feature', 'dot2', 'attr']]
sno_bed.to_csv('temp_sno_bed_seq_'+snakemake.wildcards.species, index=False, header=None, sep='\t')
sp.call('sort -k1,1 -k2,2n temp_sno_bed_seq_'+snakemake.wildcards.species+' > temp_sno_bed_seq_'+snakemake.wildcards.species+'.sorted.bed', shell=True)
bed = BedTool('temp_sno_bed_seq_'+snakemake.wildcards.species+'.sorted.bed')
merged_bed = bed.merge(s=True, c=[4, 5, 6, 8, 9], o='distinct').saveas('temp_sno_bed_seq_'+snakemake.wildcards.species+'.sorted_merged.bed')

# Update snoRNA sequence from the most recent genome 
fasta = merged_bed.sequence(fi=snakemake.input.genome, nameOnly=True, s=True)
seq_dict = {}
with open(fasta.seqfn, 'r') as fasta_file:
    for line in fasta_file:
        if '>' in line:
            sno_name = line.strip('\n').replace('>', '').replace('(+)', '').replace('(-)', '')
        else:
            seq = line.strip('\n')
            seq_dict[sno_name] = seq


# Retrieve the location and sequence of merged entries
merged_df = pd.read_csv('temp_sno_bed_seq_'+snakemake.wildcards.species+'.sorted_merged.bed', sep='\t', 
                        names=['chr', 'start', 'end', 'gene_name', 'dot', 'strand', 'feature', 'dot2'])

# Add sequence col to final_df
merged_df['sequence'] = merged_df['gene_name'].map(seq_dict)


# Extend each H/ACA snoRNA sequence to obtain 211 nt (i.e. longest expressed CD + 15 nt up/downstream)
ext_seq_dict = {}
for i_row in merged_bed:
    length_haca = int(i_row[2]) - int(i_row[1])
    difference = fixed_length - length_haca
    if difference >= 0:
        if difference % 2 == 0: # even number: split the remaining nt equally each side of the ncRNA
            l_extension = int(difference / 2)
            r_extension = int(l_extension)
        else: # odd number: split the remaining nt almost equally each side of the ncRNA (1 nt more on the 3'end)
            l_extension = int((difference - 1) / 2)
            r_extension = int(l_extension + 1)
        extended_haca = BedTool(str(i_row), from_string=True).slop(g=chr_size, r=r_extension, l=l_extension)
        extended_sequences_fa = extended_haca.sequence(fi=genome, nameOnly=True, s=True)
    
        with open(extended_sequences_fa.seqfn, 'r') as fa:
            for line in fa:
                if '>' in line:
                    haca_name = line.strip('\n').replace('>', '').replace('(+)', '').replace('(-)', '')
                    haca_name = re.sub("\(.*\)", "", haca_name)
                else:
                    seq = line.strip('\n')
                    ext_seq_dict[haca_name] = seq


merged_df[f'extended_{fixed_length}nt_sequence'] = merged_df['gene_name'].map(ext_seq_dict)

# Remove offset of 1 nt created by bedtools when converting to bed
merged_df['start'] = merged_df['start'] + 1

# Add genome_version, species name and classification
attributes_dict = snakemake.params.dataset_attribute
attributes_dict = {k.split('_cd_')[0]:v for k,v in attributes_dict.items()} # simplify species name
attributes_dict = attributes_dict[snakemake.wildcards.species]  # select only dict of species of interest
interesting_att = ['genome_version', 'species_name', 'species_classification']
attributes_dict = {k:v for k,v in attributes_dict.items() if k in interesting_att}
for k,v in attributes_dict.items():
    merged_df[k] = v

merged_df = merged_df[['gene_name', 'chr', 'strand', 'start', 'end', 
                        'genome_version', 'species_name', 'species_classification', 
                        'sequence', f'extended_{fixed_length}nt_sequence']]

merged_df['gene_name'] = merged_df['gene_name'].str.replace(',', '_')

# Remove long H/ACA (> 211 nt) so that their extended sequence is NaN
merged_df = merged_df[~merged_df[f'extended_{fixed_length}nt_sequence'].isna()]

sp.call('rm temp_sno_bed_seq_'+snakemake.wildcards.species+'*', shell=True)
merged_df.to_csv(snakemake.output.df, sep='\t', index=False)
