#!/usr/bin/python3
import pandas as pd
import collections as coll
from pybedtools import BedTool
from gtfparse import read_gtf
import subprocess as sp

tpm_df = pd.read_csv(snakemake.input.tpm_df, sep='\t')
#gtf = read_gtf(snakemake.input.gtf)
gtf = pd.read_csv(snakemake.input.gtf, sep='\t', skiprows=4, names=['chr', 'source', 'feature', 'start', 
                                                        'end', 'dot', 'strand', 'dot2', 'attributes'])
gtf = gtf[gtf['feature'] == 'gene']
gtf[['a', 'gene_id', 'c', 'd', 'e', 'f', 'g', 'h', 'i']] = gtf['attributes'].str.split('"', expand=True)
genome = snakemake.input.genome
extension = snakemake.params.extension
species = snakemake.wildcards.species

# Clean up df columns
tpm_df = tpm_df.rename(columns={'FBgn ID': 'gene_id', 'Expression TGIRT': 'expression_TGIRT', 'New gene symbol': 'gene_name'})
tpm_df = tpm_df[['gene_id', 'gene_name', 'expression_TGIRT']]
tpm_df['gene_name'] = tpm_df['gene_name'].fillna(tpm_df['gene_id'])


# Get expressed snoRNAs (expressed in Head, Ovaries and/or S2R cell line (H, O, S))
tpm_df.loc[tpm_df['expression_TGIRT'].isin(['HOS', 'HO']), 'gene_biotype'] = 'expressed_CD_snoRNA'
tpm_df['gene_biotype'] = tpm_df['gene_biotype'].fillna('snoRNA_pseudogene')
tpm_df = tpm_df.drop(columns='expression_TGIRT')

# Create bed of these snoRNAs from genomic location info found in the gtf

sno_df = tpm_df.merge(gtf[['gene_id', 'chr', 'start', 'end', 'strand', 'dot']], on='gene_id', how='left')
sno_df = sno_df[['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'gene_biotype', 'gene_name']]
sno_df['start'] = sno_df['start'] - 1  # tweak because bed is 0 based whereas gtf is 1-based
sno_df[['start', 'end']] = sno_df[['start', 'end']].astype(int)

sno_df.to_csv(f'sno_bed_temp_{species}.tsv', sep='\t', index=False, header=False)
sno_bed = BedTool(f'sno_bed_temp_{species}.tsv')


# Get sequence of C/D snoRNAs
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

sno_df['sequence'] = sno_df['gene_id'].map(d)



# Get species name as a column
sno_df['species_name'] = str(snakemake.wildcards.species)

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
sno_df = sno_df[['gene_id', 'gene_name', 'gene_biotype', 'species_name', 'chr', 
                'strand', 'start', 'end', 'sequence', 'extended_sequence']]

# Select only expressed snoRNAs
expressed_sno_df = sno_df[sno_df['gene_biotype'] == 'expressed_CD_snoRNA']
expressed_sno_df.to_csv(snakemake.output.expressed_sno_df, sep='\t', index=False)

# Get snoRNA pseudogenes in drosophila only 
if snakemake.wildcards.species == 'drosophila_melanogaster':
    sno_pseudogenes = sno_df[sno_df['gene_biotype'] == 'snoRNA_pseudogene']
    sno_pseudogenes.to_csv(snakemake.output.droso_pseudosno, sep='\t', index=False)



# Remove temp file
sp.call(f'rm sno_bed_temp_{species}.tsv', shell=True)
