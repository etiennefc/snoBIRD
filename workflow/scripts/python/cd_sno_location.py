#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp

cols = ['gene_id', 'chr', 'strand', 
        'start', 'end', 'sequence', 
        'extended_sequence', 'species_name']
bed = pd.read_csv(snakemake.input.sno_df, sep='\t')
box_df = pd.read_csv(snakemake.input.box_location, sep='\t')
cd_fa = snakemake.input.sno_sequences
chr_size_file = snakemake.input.chr_size
genome = snakemake.input.genome
species_name_dict = snakemake.params.species_name
output_fa_c = snakemake.output.fifteen_upstream_c_fa
output_fa_d = snakemake.output.fifteen_downstream_d_fa

# Get sequences in dict
d = {}
with open(cd_fa, 'r') as file:
    for line in file:
        if line.startswith('>'):
            sno_id = line.replace('>', '').replace('\n', '')
        else:
            seq = line.replace('\n', '')
            d[sno_id] = seq
print(box_df)
print(bed)

# Merge sno genomic location and sequence to box_df
box_df = box_df.merge(bed[['chr', 'start', 'end', 'strand', 'gene_id', 'species_name']], 
                    on='gene_id', how='left')
box_df[['start', 'end']] = box_df[['start', 'end']].astype(int)
box_df['sno_sequence'] = box_df['gene_id'].map(d)
print(box_df)

# Get start of C box for sno on both strands
temp_bed = box_df.copy()
temp_bed = temp_bed[temp_bed['C_sequence'] != 'NNNNNNN']  # drop sno for which no C box was found
temp_bed['dot'] = '.'
pos, neg = temp_bed[temp_bed['strand'] == '+'], temp_bed[temp_bed['strand'] == '-']
pos['start_mod'] = pos.start + pos.C_start - 1
pos = pos.rename(columns={'start': 'sno_start', 'start_mod': 'start'})
neg['end_mod'] = neg.end - neg.C_start + 1
neg = neg.rename(columns={'end': 'sno_end', 'end_mod': 'end'})
c = ['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'species_name', 'C_sequence']
final_temp_bed = pd.concat([pos[c], neg[c]])

# Replace short species name for full name
final_temp_bed['species_name'] = final_temp_bed['species_name'].replace(species_name_dict)
final_temp_bed.to_csv('temp_total.bed', index=False, header=False, sep='\t')

# Get start of D box for sno on both strands
temp_bed2 = box_df.copy()
temp_bed2 = temp_bed2[temp_bed2['D_sequence'] != 'NNNN']  # drop sno for which no D box was found
temp_bed2['dot'] = '.'
pos2, neg2 = temp_bed2[temp_bed2['strand'] == '+'], temp_bed2[temp_bed2['strand'] == '-']
pos2['start_mod'] = pos2.start + pos2.C_start - 1
pos2 = pos2.rename(columns={'start': 'sno_start', 'start_mod': 'start'})
neg2['end_mod'] = neg2.end - neg2.C_start + 1
neg2 = neg2.rename(columns={'end': 'sno_end', 'end_mod': 'end'})
d = ['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'species_name', 'D_sequence']
final_temp_bed_d = pd.concat([pos2[d], neg2[d]])

# Replace short species name for full name
final_temp_bed_d['species_name'] = final_temp_bed_d['species_name'].replace(species_name_dict)
final_temp_bed_d.to_csv('temp_total2.bed', index=False, header=False, sep='\t')

# Get the 15 nt up_downstream of each C_D box per species
cd_bed = BedTool('temp_total.bed')
cd_bed2 = BedTool('temp_total2.bed')
for species in species_name_dict.values():
    cd_bed_species = BedTool(line for line in cd_bed if line[6] == species)
    cd_bed_species2 = BedTool(line for line in cd_bed2 if line[6] == species)
    chr_size_sp = [path for path in chr_size_file if species in path][0]
    genome_sp = [path for path in genome if species in path][0]
    flanking = cd_bed_species.flank(g=chr_size_sp, l=15, r=0, s=True) 
    fifteen_nt_seq = flanking.sequence(fi=genome_sp, nameOnly=True, s=True)
    flanking2 = cd_bed_species2.flank(g=chr_size_sp, l=0, r=15, s=True) 
    fifteen_nt_seq2 = flanking2.sequence(fi=genome_sp, nameOnly=True, s=True)

    with open(fifteen_nt_seq.seqfn, 'r') as f, open(output_fa_c, 'a') as out_file_c:
        for line in f:
            if line.startswith('>'):
                id_ = line.replace('>', '').split('(')[0]
                out_file_c.write(f'>{id_}\n')
            else:
                flanking_seq = line.replace('\n', '')
                out_file_c.write(line)
    
    with open(fifteen_nt_seq2.seqfn, 'r') as f2, open(output_fa_d, 'a') as out_file_d:
        for line in f2:
            if line.startswith('>'):
                id_ = line.replace('>', '').split('(')[0]
                out_file_d.write(f'>{id_}\n')
            else:
                flanking_seq = line.replace('\n', '')
                out_file_d.write(line)




sp.call('rm temp_total*.bed', shell=True)