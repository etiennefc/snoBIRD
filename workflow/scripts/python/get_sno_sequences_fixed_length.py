#!/usr/bin/python3
import pandas as pd
import subprocess as sp
import collections as coll
from pybedtools import BedTool
import subprocess as sp 

fixed_length = int(snakemake.wildcards.fixed_length)
species_short_name_dict = snakemake.params.species_short_name
genomes = snakemake.input.genomes
chr_size = snakemake.input.chr_size
sno_literature = pd.read_csv(snakemake.input.sno_literature, sep='\t')
sno_literature = sno_literature[['gene_id', 'chr', 'strand', 'start', 
                    'end', 'sequence', 'extended_sequence', 'species_name']]

# Remove Drosophila C/D from Huang et al. 2005 (since they overlap with those from Sklias et al 2024)
sno_literature = sno_literature[sno_literature['species_name'] != 'D_melanogaster']
sno_literature['species_name'] = sno_literature['species_name'].map(species_short_name_dict)

# Load Drosophila expressed C/D (in TGIRT-Seq)
droso = pd.read_csv(snakemake.input.sno_tgirt_droso[0], sep='\t')
droso = droso[['gene_id', 'chr', 'strand', 'start', 'end', 'sequence', 
            'extended_sequence', 'species_name']]

# Load human, mouse and S. cerevisiae dfs (C/D expressed in TGIRT-Seq)
sno_dfs = [sno_literature, droso]
for path in snakemake.input.sno_tgirt:
    df = pd.read_csv(path, sep='\t')
    sp_name = path.split('/')[-1].split('_expressed_')[0]
    df['species_name'] = sp_name
    df = df[['gene_id', 'chr', 'strand', 'start', 'end', 'sequence', 
            'extended_sequence', 'species_name']]
    sno_dfs.append(df)

# Concat all cd dfs together
all_cd_df = pd.concat(sno_dfs)

# Extend to fixed_length nt total around each snoRNA 
# (fixed_length nt (i.e. 198 nt) = length of the longest 
# expressed CD below the 95th percentile (168 nt) + 15 nt up/down stream)
ext_seq_dict = {}
for species in pd.unique(all_cd_df.species_name):
    genome_path = [path for path in genomes if species in path][0]
    chr_size_path = [path for path in chr_size if species in path][0]
    species_df = all_cd_df[all_cd_df.species_name == species]
    species_df['dot'], species_df['dot2'], species_df['feature'] = '.', '.', 'gene'
    species_df['gene_biotype'] = 'expressed_cd'

    species_df['length'] = species_df.end - species_df.start # sno length
    species_df['diff'] = fixed_length - species_df.length  # remaining nt to extend and split betwen the sno 5'/3' ends
    for row in species_df.iterrows():
        row = pd.DataFrame(row[1]).T.reset_index(drop=True) # convert row to dataframe
        if row.loc[0, 'diff'] % 2 == 0: # even number: split the remaining nt equally each side of the sno
            l_extension = int(row['diff'] / 2)
            r_extension = int(l_extension)
        else: # odd number: split the remaining nt almost equally each side of the sno (1 nt more on the 3'end)
            l_extension = int((row['diff'] - 1) / 2)
            r_extension = int(l_extension + 1)
        bed = row[['chr', 'start', 'end', 'gene_id', 'dot', 'strand']]
        bed.to_csv(f'cd_fixed_bed_temp_{species}.tsv', sep='\t', index=False, header=False)
        sno_bed = BedTool(f'cd_fixed_bed_temp_{species}.tsv')
        extended_sno_bed = sno_bed.slop(g=chr_size_path, l=l_extension, r=r_extension)  # extend sequence to obtain 211 nt
        extended_sequence_fa = extended_sno_bed.sequence(fi=genome_path, s=True, nameOnly=True)

        with open(extended_sequence_fa.seqfn, 'r') as fa:
            for line in fa:
                if '>' in line:
                    sno_name = line.strip('\n').replace('>', '').replace('(+)', '').replace('(-)', '')
                else:
                    seq = line.strip('\n')
                    ext_seq_dict[sno_name] = seq

        sp.call(f'rm cd_fixed_bed_temp_{species}.tsv', shell=True)



# Create extended_211nt_sequence columns
all_cd_df[f'extended_{fixed_length}nt_sequence'] = all_cd_df['gene_id'].map(ext_seq_dict)
all_cd_df = all_cd_df.drop(columns='extended_sequence')
all_cd_df.to_csv(snakemake.output.df, sep='\t', index=False)

# Create fasta of all extended expressed C/D sno sequences with fixed length of 211 nt
d = dict(zip(all_cd_df.gene_id, all_cd_df[f'extended_{fixed_length}nt_sequence']))

with open(snakemake.output.fa, 'w') as f:
    for id, sequence in d.items():
        f.write(f'>{id}\n')
        f.write(f'{sequence}\n')
