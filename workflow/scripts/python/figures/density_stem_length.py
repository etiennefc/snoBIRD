#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
import subprocess as sp 
import re
import regex

fixed_length = snakemake.wildcards.fixed_length
expressed_cd = 'temp_expressed_cd.fa'
expressed_cd_df = pd.read_csv(snakemake.input.positives_df, sep='\t')
path = snakemake.output.density

def cut_sequence(seq):
    # Get the 20 first and 20 last nt of a given sequence
    first, last = seq[:20], seq[-20:]
    length = len(seq)
    return first, last, length

def find_d_box(seq):
    """ Find exact D box (CUGA), if not present, find D box with 1 or max 2
        substitutions. Return also the start and end position of that box as
        1-based values. If no D box is found, return a 'NNNN' empty D box and 0
        as start and end of D box."""
    first_20, last_20, length_seq = cut_sequence(seq)
    len_d_box = 4
    # First, find exact D box (CUGA) within 20 nt of the snoRNA 3' end
    if re.search('CUGA', last_20) is not None:  # find exact D box
        *_, last_possible_d = re.finditer('CUGA', last_20)
        d_motif = last_possible_d.group(0)  # if multiple exact D boxes found, keep the D box closest to 3' end
        d_start = (length_seq - 20) + last_possible_d.start() + 1
        d_end = (length_seq - 20) + last_possible_d.end()
        return d_motif, d_start, d_end
    else:  # find not exact D box (up to max 50% of substitution allowed (i.e. 2 nt))
        for sub in range(1, int(len_d_box/2 + 1)):  # iterate over 1 to 2 substitutions allowed
            d_motif = regex.findall("(CUGA){s<="+str(sub)+"}", last_20, overlapped=True)
            if len(d_motif) >= 1:  # if we have a match, break and keep that match (1 sub privileged over 2 subs)
                d_motif = d_motif[-1]  # if multiple D boxes found, keep the the D box closest to 3' end
                d_start = (length_seq - 20) + last_20.rindex(d_motif) + 1
                d_end = d_start + len(d_motif) - 1
                return d_motif, d_start, d_end  # this exits the global else statement
        # If no D box is found, return NNNN and 0, 0 as D box sequence, start and end
        d_motif, d_start, d_end = 'NNNN', 0, 0
        return d_motif, d_start, d_end


def find_c_box(seq):
    """ Find exact C box (RUGAUGA, where R is A or G), if not present, find C
        box with 1,2 or max 3 substitutions. Return also the start and end
        position of that box as 1-based values. If no C box is found, return a
        'NNNNNNN' empty C box and 0 as start and end of C box."""
    first_20, last_20, length_seq = cut_sequence(seq)
    len_c_box = 7
    # First, find exact C box (RUGAUGA) within 20 nt of the snoRNA 5' end
    if re.search('(A|G)UGAUGA', first_20) is not None:  # find exact C box
        i = 1
        for possible_c in re.finditer('(A|G)UGAUGA', first_20):
            if i <= 1:  # select first matched group only (closest RUGAUGA to 5' end of snoRNA)
                c_motif = possible_c.group(0)
                c_start = possible_c.start() + 1
                c_end = possible_c.end()
                i += 1
                return c_motif, c_start, c_end  # this exits the global if statement
    else:  # find not exact C box (up to max 3 substitution allowed)
        for sub in range(1, int((len_c_box-1)/2 + 1)):  # iterate over 1 to 3 substitutions allowed
            c_motif = regex.findall("((A|G)UGAUGA){s<="+str(sub)+"}", first_20, overlapped=True)
            if len(c_motif) >= 1:  # if we have a match, break and keep that match (1 sub privileged over 2 subs)
                c_motif = c_motif[0][0]  # if multiple C boxes found, keep the the C box closest to 5' end
                c_start = first_20.find(c_motif) + 1
                c_end = c_start + len(c_motif) - 1
                return c_motif, c_start, c_end  # this exits the global else statement
        # If no C box is found, return NNNNNNN and 0, 0 as C box sequence, start and end
        c_motif, c_start, c_end = 'NNNNNNN', 0, 0
        return c_motif, c_start, c_end

def generate_df(fasta, func, motif_name):
    """ From a fasta of snoRNA sequences, find a given motif (C or D) using
        predefined function func and output the motif sequence, start and end
        as a df."""
    # Get motif, start and end position inside dict
    box_dict = {}
    with open(fasta, 'r') as f:
        sno_id = ''
        for line in f:
            if line.startswith('>'):
                id = line.lstrip('>').rstrip('\n')
                sno_id = id
            else:
                seq = line.rstrip('\n')
                motif, start, end = func(seq)
                box_dict[sno_id] = [motif, start, end]

    # Create dataframe from box_dict
    box = pd.DataFrame.from_dict(box_dict, orient='index',
                                columns=[f'{motif_name}_sequence', f'{motif_name}_start',
                                        f'{motif_name}_end'])
    box = box.reset_index()
    box = box.rename(columns={"index": "gene_id"})
    return box


# Create temp fasta
len_dict = {}
with open(expressed_cd, 'w') as f:
    for id, seq in dict(zip(expressed_cd_df.gene_id, expressed_cd_df.sequence)).items():
        f.write(f'>{id}\n{seq}\n')
        len_dict[id] = len(seq)

# Find C and D box positions in the expressed snoRNA
df_c = generate_df(expressed_cd, find_c_box, 'C')
df_c = df_c[df_c['C_sequence'] != 'NNNNNNN']
df_d = generate_df(expressed_cd, find_d_box, 'D')
df_d = df_d[df_d['D_sequence'] != 'NNNN']

# Merge sno length to df_d
df_d['sno_length'] = df_d['gene_id'].map(len_dict)

# Get the nb of nt between the start of the sno and the start of the C box
df_c['up_stem_len'] = df_c.C_start - 1

# Get the nb of nt between the end of the D box and the end of the sno
df_d['down_stem_len'] = df_d.sno_length - df_d.D_end

# Create density plot
title = 'Distribution of the distance between the start/end\nof the snoRNA and the start/end of the C/D boxes'
ft.density_x([df_c.up_stem_len, df_d.down_stem_len], 'Distance between box and\nsnoRNA boundaries (nt)', 
            'Density', 'linear', title, ['blue', 'red'], ['C box distance', 'D box distance'], path, 
            xvline=5, yminvline=0, ymaxvline=0.22)

print(f'C distance median: {df_c.up_stem_len.median()} nt')
print(f'C distance mean: {df_c.up_stem_len.mean()} nt')
print(f'D distance median: {df_d.down_stem_len.median()} nt')
print(f'D distance mean: {df_d.down_stem_len.mean()} nt')


sp.call(f'rm {expressed_cd}', shell=True)