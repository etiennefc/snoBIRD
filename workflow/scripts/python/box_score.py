#!/usr/bin/python3
import pandas as pd
import re
import regex
import matplotlib.pyplot as plt
import seaborn as sns 
from functools import reduce
from math import ceil, floor
from figures import functions as ft
import subprocess as sp

""" Find the most probable C, D, C' and D' boxes in a sequence (if they exist). 
    Return box score (0:perfect box score, higher:more degenerate motifs)."""
expressed_cd_df = pd.read_csv(snakemake.input.positives_fa, sep='\t')
negatives_df = pd.concat(([pd.read_csv(path, sep='\t') for path in snakemake.input.negatives_fa]))
output_pos, output_neg = snakemake.output.positives, snakemake.output.negatives
fixed_length = int(snakemake.wildcards.fixed_length)
first_, last_ = floor(fixed_length/2), ceil(fixed_length/2)  # to get index of 1st and 2nd half of sequence

# Create temporary fastas
expressed_cd = "expressed_cd.fa"
rest = "negatives.fa"
with open(expressed_cd, 'w') as f:
    for id, seq in dict(zip(expressed_cd_df.gene_id, expressed_cd_df[f'extended_{fixed_length}nt_sequence'])).items():
        f.write(f'>{id}\n{seq}\n')

with open(rest, 'w') as f:
    for id, seq in dict(zip(negatives_df.gene_id, negatives_df[f'extended_{fixed_length}nt_sequence'])).items():
        f.write(f'>{id}\n{seq}\n')



def cut_sequence(seq):
    # Get the 105 first (without the first 15 nt which are external terminal stems) 
    # and 106 last nt (without the last 15 nt which are external terminal stems) of a given sequence
    seq = seq.replace('U', 'T')  # convert RNA to DNA
    first, last = seq[15:first_], seq[-last_:-15]
    length = len(seq)
    return first, last, length

def find_c_box(seq):
    """ Find exact C box (RUGAUGA, where R is A or G), if not present, find C
        box with 1,2 or max 3 substitutions. Return also the start and end
        position of that box as 1-based values. If no C box is found, return a
        'NNNNNNN' empty C box and 0 as start and end of C box."""
    first_range, last_range, length_seq = cut_sequence(seq)
    len_c_box = 7
    # First, find exact C box (RUGAUGA) within range to the snoRNA 5' end
    if re.search('(A|G)TGATGA', first_range) is not None:  # find exact C box
        i = 1
        for possible_c in re.finditer('(A|G)TGATGA', first_range):
            if i <= 1:  # select first matched group only (closest RUGAUGA to 5' end of snoRNA)
                c_motif = possible_c.group(0)
                c_start = possible_c.start() + 1 + 15
                c_end = c_start + len_c_box - 1
                i += 1
                return c_motif, c_start, c_end  # this exits the global if statement
    else:  # find not exact C box (up to max 3 substitution allowed)
        for sub in range(1, int((len_c_box-1)/2 + 1)):  # iterate over 1 to 3 substitutions allowed
            c_motif = regex.findall("((A|G)TGATGA){s<="+str(sub)+"}", first_range, overlapped=True)
            if len(c_motif) >= 1:  # if we have a match, break and keep that match (1 sub privileged over 2 subs)
                c_motif = c_motif[0][0]  # if multiple C boxes found, keep the the C box closest to 5' end
                c_start = first_range.find(c_motif) + 1 + 15
                c_end = c_start + len_c_box - 1
                return c_motif, c_start, c_end  # this exits the global else statement
        # If no C box is found, return NNNNNNN and 0, 0 as C box sequence, start and end
        c_motif, c_start, c_end = 'NNNNNNN', 0, 0
        return c_motif, c_start, c_end


def find_d_box(seq, c_end):
    """ Given a C box position, find exact D box (CUGA), if not present, find D box with 1 or max 2
        substitutions. Return also the start and end position of that box as
        1-based values. If no D box is found, return a 'NNNN' empty D box and 0
        as start and end of D box."""
    first_range, last_range, length_seq = cut_sequence(seq)
    len_d_box = 4
    d_motif = ''
    # First, find exact D box (CUGA) within range to the snoRNA 3' end
    if re.search('CTGA', last_range) is not None:  # find exact D box
        all_perfect_d = [d for d in re.finditer('CTGA', last_range) if d.start() + len(first_range) + 15 > c_end + 30]  # >30nt from the C box (smallest distance between C and D boxes in known snoRNAs)
        if len(all_perfect_d) > 0:
            last_d = all_perfect_d[-1]
            d_motif = last_d.group(0)  # if multiple exact D boxes found, keep the D box closest to 3' end
            d_start = (length_seq - len(last_range) -15) + last_d.start() + 1
            d_end = d_start + len_d_box - 1
            return d_motif, d_start, d_end
    if d_motif == '':  # find not exact D box (up to max 50% of substitution allowed (i.e. 2 nt))
        for sub in range(1, int(len_d_box/2 + 1)):  # iterate over 1 to 2 substitutions allowed
            d_motifs = regex.findall("(CTGA){s<="+str(sub)+"}", last_range, overlapped=True)
            all_possible_d = [d for d in d_motifs if (length_seq - len(last_range) + last_range.rindex(d) + 1) > c_end + 30]  # >30nt from the C box
            if len(all_possible_d) > 0:  # if we have a match, break and keep that match (1 sub privileged over 2 subs)
                d_motif = all_possible_d[-1]  # if multiple D boxes found, keep the the D box closest to 3' end
                d_start = (length_seq - len(last_range) - 15) + last_range.rindex(d_motif) + 1
                d_end = d_start + len_d_box - 1
                return d_motif, d_start, d_end  # this exits the global else statement
        # If no D box is found, return NNNN and 0, 0 as D box sequence, start and end
        d_motif, d_start, d_end = 'NNNN', 0, 0
        return d_motif, d_start, d_end


def find_c_prime_d_prime_hamming(seq, c_end, d_start):
    """ Find best C'/D' pair that minimizes Hamming distance compared to consensus C'/D' motif. """
    # Get the nucleotides between the end of the C box (+2 nt space) and the start of the D box (-2 nt space)
    middle_seq = seq[(c_end + 2):(d_start-2)]
    len_c_prime_box, len_d_prime_box = 7, 4
    

    # Find all possible C' boxes and their start/end and compute their Hamming distance compared to consensus motif
    hamming_c_prime, all_c_primes, all_c_primes_start, all_c_primes_end, temp_motif = [], [], [], [], ''
    for sub in range(0, int((len_c_prime_box-1)/2 + 1)):
        c_prime_motif = regex.findall("((A|G)TGATGA){s<="+str(sub)+"}", middle_seq, overlapped=True)
        if len(c_prime_motif) >= 1:
            for motif in c_prime_motif:
                if motif[0] != temp_motif:  # to avoid repeated motif between 0 and 1 substitution allowed
                    all_c_primes.append(motif[0])
                    hamming_c_prime.append(sub)
                    c_prime_start = middle_seq.index(motif[0]) + c_end + 1 + 2  # +2 because of the 2 nt space created in middle_seq
                    c_prime_end = c_prime_start + len_c_prime_box - 1
                    all_c_primes_start.append(c_prime_start)
                    all_c_primes_end.append(c_prime_end)                    
                    temp_motif = motif[0]              
    if len(all_c_primes) == 0:  # if no C' box was found, hamming distance is 7, C' motif is NNNNNNN and start and end are 0
        hamming_c_prime, all_c_primes, all_c_primes_start, all_c_primes_end = [7], ['NNNNNNN'], [0], [0]
    
    
    # Find all possible D' boxes and their start/end and compute their Hamming distance compared to consensus motif
    hamming_d_prime, all_d_primes, all_d_primes_start, all_d_primes_end, temp_motif = [], [], [], [], ''
    for sub in range(0, int((len_d_prime_box)/2 + 1)):
        d_prime_motif = regex.findall("(CTGA){s<="+str(sub)+"}", middle_seq, overlapped=True)
        if len(d_prime_motif) >= 1:
            for motif in d_prime_motif:
                if motif != temp_motif:  # to avoid repeated motif between 0 and 1 substitution allowed
                    all_d_primes.append(motif)
                    hamming_d_prime.append(sub)
                    d_prime_start = middle_seq.index(motif) + c_end + 1 + 2
                    d_prime_end = d_prime_start + len_d_prime_box - 1
                    all_d_primes_start.append(d_prime_start)
                    all_d_primes_end.append(d_prime_end)
                    temp_motif = motif                
    if len(all_d_primes) == 0:  # if no D' box was found, hamming distance is 4, D' motif is NNNN and start and end are 0
        hamming_d_prime, all_d_primes, all_d_primes_start, all_d_primes_end = [4], ['NNNN'], [0], [0]
        
    # Find all possible D'-C' pairs where C' is downstream of D' by at least 2 nt (i.e. at least 1 nt between D' and C' boxes) 
    # and return the best pair according to the lowest total Hamming distance (if two pairs have the same Hamming distance, the closest one to 5' of snoRNA is chosen)
    total_hamming, d_prime_index, c_prime_index = 10000000000, 10000000000, 10000000000  # these are dummy high numbers 
    for i, d_prime in enumerate(all_d_primes):
        for j, c_prime in enumerate(all_c_primes):
            # C' downstream of D' by at least 2 nt or if D' is found but not C' (the case where no D' is found but C' is found is also included: d_prime_end = 0 is smaller than any existing c_prime_start)
            if (all_d_primes_end[i] +2 <= all_c_primes_start[j]) | ((all_d_primes[i] != 'NNNN') & (all_c_primes[j] == 'NNNNNNN')):   
                temp_total_hamming = hamming_d_prime[i] + hamming_c_prime[j]
                if temp_total_hamming < total_hamming:
                    total_hamming = temp_total_hamming
                    d_prime_index = i
                    c_prime_index = j

            # if no D' nor C' box are found
            elif (all_d_primes[i] != 'NNNN') & (all_c_primes[j] == 'NNNNNNN'):
                temp_total_hamming = hamming_d_prime[i] + hamming_c_prime[j]
                d_prime_index, c_prime_index = 0, 0


    # If only one C' box is found and multiple D' are found  but are overlapping, keep the D' box with the 
    # lowest Hamming distance (closest to 5' if multiple D' have the same Hamming) and return NNNNNNN as the C'
    if (d_prime_index == 10000000000) | (c_prime_index == 10000000000):
        all_c_primes, all_c_primes_start, all_c_primes_end, c_prime_index = ['NNNNNNN'], [0], [0], 0
        d_prime_index = hamming_d_prime.index(min(hamming_d_prime))
    c_prime_motif, c_prime_start, c_prime_end = all_c_primes[c_prime_index], all_c_primes_start[c_prime_index], all_c_primes_end[c_prime_index]
    d_prime_motif, d_prime_start, d_prime_end = all_d_primes[d_prime_index], all_d_primes_start[d_prime_index], all_d_primes_end[d_prime_index]
    return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end



def generate_df_c_d(fasta):
    """ From a fasta of snoRNA sequences, find a given motif (C and D) using
        predefined function functions and output the motif sequence, start and end
        as separate dicts."""
    # Get motif, start and end position inside dict
    box_dict_c, box_dict_d = {}, {}
    with open(fasta, 'r') as f:
        sno_id = ''
        for line in f:
            if line.startswith('>'):
                id = line.lstrip('>').rstrip('\n')
                sno_id = id
            else:
                seq = line.rstrip('\n')
                c_motif, c_start, c_end = find_c_box(seq)
                box_dict_c[sno_id] = [c_motif, c_start, c_end]

                d_motif, d_start, d_end = find_d_box(seq, c_end)
                box_dict_d[sno_id] = [d_motif, d_start, d_end]


    return box_dict_c, box_dict_d


def generate_df_prime(fasta, c_box_dict, d_box_dict):
    """ From a fasta of snoRNA sequences, find a given motif (C' or D') using
        predefined function find_c_prime_d_prime_hamming and output the motif sequence,
        start and end as a df. USe c_box_dict and d_box_dict to define where the 
        C and D box start (and thereby the search space for the C'/D' boxes)."""
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
                if d_box_dict[sno_id][0] == "NNNN":
                    d_start = len(seq) - 4  # at least 4 nt to make a whole D' box from the end of the seq
                else:
                    d_start = d_box_dict[sno_id][1] 
                c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end = find_c_prime_d_prime_hamming(seq, c_box_dict[sno_id][2], d_start)
                box_dict[sno_id] = [c_prime_motif, c_prime_start,
                                    c_prime_end, d_prime_motif, d_prime_start,
                                    d_prime_end]

    # Create dataframe from box_dict
    box = pd.DataFrame.from_dict(box_dict, orient='index',
                                columns=['C_prime_sequence', 'C_prime_start',
                                        'C_prime_end', 'D_prime_sequence',
                                        'D_prime_start', 'D_prime_end'])
    box = box.reset_index()
    box = box.rename(columns={"index": "gene_id"})
    return box


def find_all_boxes(fasta):
    """ Find C, D, C' and D' boxes in given fasta using generate_df and concat
        resulting dfs horizontally."""
    c_box_dict, d_box_dict = generate_df_c_d(fasta)
    df_c = pd.DataFrame.from_dict(c_box_dict, orient='index',
                                columns=['C_sequence', 'C_start',
                                        'C_end'])
    df_c = df_c.reset_index()
    df_c = df_c.rename(columns={"index": "gene_id"})
    df_d = pd.DataFrame.from_dict(d_box_dict, orient='index',
                                columns=['D_sequence', 'D_start',
                                        'D_end'])
    df_d = df_d.reset_index()
    df_d = df_d.rename(columns={"index": "gene_id"})
    df_c_prime_d_prime = generate_df_prime(fasta, c_box_dict, d_box_dict)

    df_final = reduce(lambda  left,right: pd.merge(left,right,on=['gene_id'],
                                            how='outer'),
                                            [df_c, df_d, df_c_prime_d_prime])
    return df_final


def hamming(found_motif, consensus_motif):
    """ Find the hamming distance of a found motif compared to the consensus
        motif. """
    hamming = 0
    if consensus_motif.startswith('R'):  # To deal with C/C' motif
        if found_motif[0] in ['A', 'G']:
            consensus_motif = found_motif[0] + consensus_motif[1:]
    for i, char in enumerate(found_motif):
        if char != consensus_motif[i]:
            hamming += 1
    return hamming



# Find box score for expressed C/D
df = find_all_boxes(expressed_cd)
dictio = df.set_index('gene_id').filter(regex='_sequence$').to_dict('index')
consensus_dict = {'D_sequence': 'CTGA', 'C_sequence': 'RTGATGA', 
                'D_prime_sequence': 'CTGA', 'C_prime_sequence': 'RTGATGA'}

score_dict = {}
for gene_id, seq_dict in dictio.items():
    box_score = 0
    for box_name, seq in seq_dict.items():
        box_score += hamming(seq, consensus_dict[box_name])
    score_dict[gene_id] = box_score

df['box_score'] = df['gene_id'].map(score_dict)


# Find box score for the rest of examples (negatives and pseudosno)
ddf = find_all_boxes(rest)
ddictio = ddf.set_index('gene_id').filter(regex='_sequence$').to_dict('index')

score_dict = {}
for gene_id, seq_dict in ddictio.items():
    box_score = 0
    for box_name, seq in seq_dict.items():
        box_score += hamming(seq, consensus_dict[box_name])
    score_dict[gene_id] = box_score

ddf['box_score'] = ddf['gene_id'].map(score_dict)


# Save both dfs
df.to_csv(output_pos, sep='\t', index=False)
ddf.to_csv(output_neg, sep='\t', index=False)


sp.call('rm expressed_cd.fa negatives.fa', shell=True)

