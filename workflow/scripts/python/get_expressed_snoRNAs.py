#!/usr/bin/python3
import pandas as pd
import collections as coll
import numpy as np
from pybedtools import BedTool
from gtfparse import read_gtf
import subprocess as sp
from statistics import median

tpm_df = pd.read_csv(snakemake.input.tpm_df, sep='\t')
sno_df = tpm_df[tpm_df['gene_biotype'] == 'snoRNA']
sno_type_df = pd.read_csv(snakemake.input.sno_type_df, sep='\t')
gtf = read_gtf(snakemake.input.gtf)
genome = snakemake.input.genome
species = str(snakemake.wildcards.species)
extension = snakemake.params.extension

# Create bed of all snoRNAs from the gtf
sno_gtf = gtf[(gtf['feature'] == 'gene') & (gtf['gene_biotype'] == 'snoRNA')]
sno_gtf['dot'], sno_gtf['dot2'] = '.', '.'
bed = sno_gtf[['seqname', 'start', 'end', 'gene_id', 'dot', 
                    'strand', 'source', 'feature', 'dot2', 'gene_biotype']]
bed['start'] = bed['start'] - 1  # tweak because bed is 0 based whereas gtf is 1-based
bed.to_csv(f'sno_bed_temp_{species}.tsv', sep='\t', index=False, header=False)
sno_bed = BedTool(f'sno_bed_temp_{species}.tsv')

# Get sequence of snoRNAs
fasta = sno_bed.sequence(fi=genome, nameOnly=True, s=True)
d = {}
with open(fasta.seqfn, 'r') as f:
    for line in f:
        if '>' in line:
            sno_id = line.replace('>', '').replace('\n', '')
            sno_id = sno_id.replace('(+)', '').replace('(-)', '')
        else:
            seq = line.replace('\n', '')
            d[sno_id] = seq


# Select only C/D box snoRNAs
sno_type_df = sno_type_df[['gene_id', 'sno_type']]
sno_df = sno_df.merge(sno_type_df, how='left', on='gene_id')  # scaRNA present in sno_type_df are automatically dropped (Thus 1045 and not 1056 snoRNAs)
sno_df = sno_df[sno_df['sno_type'] == 'C/D']
print(sno_df)
print(len(sno_df))
print(coll.Counter(sno_type_df.sno_type))
print(coll.Counter(sno_df.sno_type))

# For this to work, columns of same condition must be following each other in df
# Get unique conditions acroos sno_df columns
cols = list(sno_df.columns)
conditions = []
for col in cols:
    if ('gene_' not in col) & ('sno' not in col):
        condition, number = col.rsplit('_', maxsplit=1)
        if condition not in conditions:
            conditions.append(condition)
print(conditions)

# Get the snoRNAs for which at least one average condition is > 1 TPM
sno_df = sno_df.reset_index(drop=True)
expressed_sno = []
for i in range(0, len(sno_df)):
    for cond in conditions:
        temp_df = sno_df.copy()
        temp_df[f'{cond}_avg'] = temp_df.filter(regex=f'^{cond}').mean(axis=1) # compute avg per condition
        if temp_df.loc[i, f'{cond}_avg'] > 1:
            if cond == 'Plasma':  # plasma samples have a lot of variability between the 17, so a second filter needs to be applied
                plasma_df = sno_df.filter(regex=f'^{cond}').loc[i, :]
                plasma_tpm = plasma_df.values.tolist()  # get all plasma TPM of a given i snoRNA
                # if at least half (9/17) of plasma samples are above 1 TPM, use Plasma cols to define if sno is expressed
                if median(plasma_tpm) > 1:  #same as if len([it for it in plasma_tpm if it > 1]) >= 9:
                    # This adds ~7 plasma-enriched snoRNAs (mostly snoU13, U3 and U8)
                    expressed_sno.append(str(temp_df.loc[i, 'gene_id']))  
                    break
            else:
                expressed_sno.append(str(temp_df.loc[i, 'gene_id']))
                break

# Get species, sno location and sequence
sno_df['species_name'] = str(snakemake.wildcards.species)
sno_df = sno_df.merge(sno_gtf[['seqname', 'strand', 'start', 'end', 'gene_id']], how='left', on='gene_id')
sno_df = sno_df.rename(columns={'seqname': 'chr'})
sno_df['sequence'] = sno_df['gene_id'].map(d)

# Get sno extended sequence (15 nt up- and downstream of the sno)
extended_sno = sno_bed.slop(g=snakemake.input.chr_size, r=extension, l=extension)
sequences_sno = extended_sno.sequence(fi=snakemake.input.genome, s=True, nameOnly=True)
extended_seq_dict = {}
with open(sequences_sno.seqfn, 'r') as f:
    for line in f:
        if '>' in line:
            cd_id = line.replace('\n', '').replace('>', '').replace('(+)', '').replace('(-)', '')
        else:
            seq_ = line.replace('\n', '')
            extended_seq_dict[cd_id] = seq_
sno_df['extended_sequence'] = sno_df['gene_id'].map(extended_seq_dict)

# Select only expressed snoRNAs
expressed_sno_df = sno_df[sno_df['gene_id'].isin(expressed_sno)]
print(len(expressed_sno_df))
print(expressed_sno_df)
expressed_sno_df.to_csv(snakemake.output.expressed_sno_df, sep='\t', index=False)

# Get snoRNA pseudogenes in human only 
if snakemake.wildcards.species == 'homo_sapiens':
    sno_pseudogenes = sno_df[~sno_df['gene_id'].isin(expressed_sno)]
    sno_pseudogenes.to_csv(snakemake.params.human_pseudosno, sep='\t', index=False)

# Get snoRNA pseudogenes in mouse only 
if snakemake.wildcards.species == 'mus_musculus':
    sno_pseudogenes = sno_df[~sno_df['gene_id'].isin(expressed_sno)]
    sno_pseudogenes.to_csv(snakemake.params.mouse_pseudosno, sep='\t', index=False)

# Remove temp file
sp.call(f'rm sno_bed_temp_{species}.tsv', shell=True)
