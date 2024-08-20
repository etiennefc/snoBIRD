#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp 
import re

output = snakemake.output.pseudogenes
chr_size = snakemake.input.mouse_chr_size[0] 
genome = snakemake.input.mouse_genome
fixed_length = int(snakemake.wildcards.fixed_length)
pseudogene_df = pd.read_csv(snakemake.params.mouse_snoRNA_pseudogenes, sep='\t')
pseudogene_df_copy = pseudogene_df.copy()
pseudogene_df_copy['dot'] = '.'
pseudogene_df_copy = pseudogene_df_copy[['chr', 'start', 'end', 'gene_id', 
                                'dot', 'strand', 'gene_biotype']]
pseudogene_df_copy.to_csv('mouse_temp_pseudogene.bed', sep='\t', index=False, header=False)
sp.call('sort -k1,1 -k2,2n mouse_temp_pseudogene.bed > mouse_temp_pseudogene_sorted.bed', shell=True)

# Load pseudogene bed and get extended sequence 
pseudogene_bed = BedTool('mouse_temp_pseudogene_sorted.bed')

ext_seq_dict = {}
for i_row in pseudogene_bed:
    length_sno_pseudogene = int(i_row[2]) - int(i_row[1])
    difference = fixed_length - length_sno_pseudogene
    if difference >= 0:
        if difference % 2 == 0: # even number: split the remaining nt equally each side of the ncRNA
            l_extension = int(difference / 2)
            r_extension = int(l_extension)
        else: # odd number: split the remaining nt almost equally each side of the ncRNA (1 nt more on the 3'end)
            l_extension = int((difference - 1) / 2)
            r_extension = int(l_extension + 1)
        extended_sno_pseudogene = BedTool(str(i_row), from_string=True).slop(g=chr_size, r=r_extension, l=l_extension)
        extended_sequences_fa = extended_sno_pseudogene.sequence(fi=genome, nameOnly=True, s=True)
    
        with open(extended_sequences_fa.seqfn, 'r') as fa:
            for line in fa:
                if '>' in line:
                    sno_pseudogene_name = line.strip('\n').replace('>', '').replace('(+)', '').replace('(-)', '')
                    sno_pseudogene_name = re.sub("\(.*\)", "", sno_pseudogene_name)
                else:
                    seq = line.strip('\n')
                    ext_seq_dict[sno_pseudogene_name] = seq

# Add extended sequence to pseudogen_df
pseudogene_df[f'extended_{fixed_length}nt_sequence'] = pseudogene_df['gene_id'].map(ext_seq_dict)

# Remove offset of 1 nt created by bedtools when converting to bed
pseudogene_df['start'] = pseudogene_df['start'] + 1

# Remove sno pseudogenes that are longer than fixed_length (i.e. those with no extended sequence obtained)
pseudogene_df = pseudogene_df[~pseudogene_df[f'extended_{fixed_length}nt_sequence'].isna()]

# Select only relevant columns 
pseudogene_df['gene_biotype'] = pseudogene_df['gene_biotype'].replace('snoRNA', 'snoRNA_pseudogene')
pseudogene_df['species_name'] = 'mus_musculus'
pseudogene_df = pseudogene_df[['gene_id', 'gene_name', 'gene_biotype', 'species_name', 'chr', 'strand', 
                                'start', 'end', 'sequence', 'extended_sequence', 
                                f'extended_{fixed_length}nt_sequence']]

pseudogene_df.to_csv(output, sep='\t', index=False)

sp.call('rm mouse_temp_pseudogene*', shell=True)
