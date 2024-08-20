#!/usr/bin/python3
import pandas as pd
import re
import regex
from functools import reduce

""" Find C, D, C' and D' boxes of each snoRNA (if they exist) and their position"""
cd_fa = snakemake.input.cd_fa
output = snakemake.output.c_d_box_location

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


def find_c_prime_d_prime(seq):
    """ Find exact C' box (RUGAUGA, where R is A or G), if not present, find C'
        box with 1,2 or max 3 substitutions. Return also the start and end
        position of that box as 1-based values. If no C' box is found, return a
        'NNNNNNN' empty C' box and 0 as start and end of C' box. Of all potential
        C' boxes, always pick the C' closest to 3' end. Find after the closest
        D' box that must be upstream of the found C' box (exact D' (CUGA), then
        1 and 2 substitutions allowed)."""
    # Get the nucleotides between the 20th and last 20th nucleotides
    middle_seq = seq[20:-20]
    len_c_prime_box, len_d_prime_box = 7, 4

    # First, find exact C' box (RUGAUGA) closest to 3' end
    if re.search('(A|G)UGAUGA', middle_seq) is not None:  # find exact C' box
        *_, last_possible_c_prime = re.finditer('(A|G)UGAUGA', middle_seq)
        c_prime_motif = last_possible_c_prime.group(0)
        c_prime_start = last_possible_c_prime.start() + 21
        c_prime_end = last_possible_c_prime.end() + 20
        # Find closest exact D' box upstream of found C' box (upstream by at least 1 nt between D' and C')
        if re.search('CUGA', seq[20:c_prime_start]) is not None:  # find exact D' box
            *_, last_possible_d_prime = re.finditer('CUGA', seq[20:c_prime_start])
            d_prime_motif = last_possible_d_prime.group(0)  # if multiple exact D' boxes found, keep the D' box closest to 3' (i.e. closest to C' box)
            d_prime_start = last_possible_d_prime.start() + 21
            d_prime_end = last_possible_d_prime.end() + 20
            return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end  # this exits the global if statement
        else:  # find not exact D' box (up to max 50% of substitution allowed (i.e. 2 nt))
            for sub in range(1, int(len_d_prime_box/2 + 1)):  # iterate over 1 to 2 substitutions allowed
                d_prime_motif = regex.findall("(CUGA){s<="+str(sub)+"}", seq[20:c_prime_start], overlapped=True)
                if len(d_prime_motif) >= 1:  # if we have a match, break and keep that match (1 sub privileged over 2 subs)
                    d_prime_motif = d_prime_motif[-1]  # if multiple D boxes found, keep the the D box closest to 3' (i.e. closest to C' box)
                    d_prime_start = seq[20:c_prime_start].rindex(d_prime_motif) + 21
                    d_prime_end = d_prime_start + len(d_prime_motif) - 1
                    return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end  # this exits the global if statement
            # If no D' box is found, return NNNN and 0, 0 as D' box sequence, start and end
            d_prime_motif, d_prime_start, d_prime_end = 'NNNN', 0, 0
            return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end

    else: # find not exact C' box (up to max 3 substitution allowed)
        for sub in range(1, int((len_c_prime_box-1)/2 + 1)):  # iterate over 1 to 3 substitutions allowed
            c_prime_motif = regex.findall("((A|G)UGAUGA){s<="+str(sub)+"}", middle_seq, overlapped=True)
            if len(c_prime_motif) >= 1:  # if we have a match, break and keep that match (1 sub privileged over 2 subs)
                c_prime_motif = c_prime_motif[0][0]  # if multiple C' boxes found, keep the the C box closest to 3' end
                c_prime_start = middle_seq.rfind(c_prime_motif) + 21
                c_prime_end = c_prime_start + len(c_prime_motif) - 1
                # Find closest exact D' box upstream of found C' box (upstream by at least 1 nt between D' and C')
                if re.search('CUGA', seq[20:c_prime_start]) is not None:  # find exact D' box
                    *_, last_possible_d_prime = re.finditer('CUGA', seq[20:c_prime_start])
                    d_prime_motif = last_possible_d_prime.group(0)  # if multiple exact D' boxes found, keep the D' box closest to 3' (i.e. closest to C' box)
                    d_prime_start = last_possible_d_prime.start() + 21
                    d_prime_end = last_possible_d_prime.end() + 20
                    return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end  # this exits the global else statement
                else:  # find closest not exact D' box (up to max 50% of substitution allowed (i.e. 2 nt))
                    for sub in range(1, int(len_d_prime_box/2 + 1)):  # iterate over 1 to 2 substitutions allowed
                        d_prime_motif = regex.findall("(CUGA){s<="+str(sub)+"}", seq[20:c_prime_start], overlapped=True)
                        if len(d_prime_motif) >= 1:  # if we have a match, break and keep that match (1 sub privileged over 2 subs)
                            d_prime_motif = d_prime_motif[-1]  # if multiple D boxes found, keep the the D box closest to 3' (i.e. closest to C' box)
                            d_prime_start = seq[20:c_prime_start].rindex(d_prime_motif) + 21
                            d_prime_end = d_prime_start + len(d_prime_motif) - 1
                            return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end  # this exits the global else statement
                    # If no D' box is found, return NNNN and 0, 0 as D' box sequence, start and end
                    d_prime_motif, d_prime_start, d_prime_end = 'NNNN', 0, 0
                    return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end  # this exits the global else statement
        # If no C' and no D' box are found, return NNNNNNN/NNNN and 0, 0 as C'/D' box sequence, start and end
        d_prime_motif, d_prime_start, d_prime_end = 'NNNN', 0, 0
        c_prime_motif, c_prime_start, c_prime_end = 'NNNNNNN', 0, 0
        return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end


def find_d_prime_c_prime(seq):
    """ Find exact D' box (CUGA) first; if not present, find D'
        box with 1 or 2 max substitutions. Return also the start and end
        position of that box as 1-based values. If no D' box is found, return a
        'NNNN' empty D' box and 0 as start and end of C box. Of all potential
        D' boxes, always pick the D' closest to 5' end. Find after the closest
        C' box that must be downstream of the found D' box (exact C' (RUGAUGA), then
        1, 2 and 3 substitutions allowed)."""
    # Get the nucleotides between the 20th and last 20th nucleotides
    middle_seq = seq[20:-20]
    len_c_prime_box, len_d_prime_box = 7, 4

    # First, find exact D' box (CUGA) closest to 5' end
    if re.search('CUGA', middle_seq) is not None:  # find exact D' box
        i = 1
        for possible_d_prime in re.finditer('CUGA', middle_seq):
            if i <= 1:  # select first matched group only (closest CUGA to 5' end of snoRNA)
                d_prime_motif = possible_d_prime.group(0)
                d_prime_start = possible_d_prime.start() + 21
                d_prime_end = possible_d_prime.end() + 20
                i += 1

        # Find closest exact C' box downstream of found D' box (downstream by at least 2 nt (1 nt between D' and C'))
        if re.search('(A|G)UGAUGA', seq[d_prime_end+2:-20]) is not None:  # find exact C' box
            j = 1
            for possible_c_prime in re.finditer('(A|G)UGAUGA', seq[d_prime_end+2:-20]):
                if j <= 1:  # select first matched group only (closest RUGAUGA to found D' box)
                    c_prime_motif = possible_c_prime.group(0)
                    c_prime_start = possible_c_prime.start() + d_prime_end
                    c_prime_end = possible_c_prime.end() + d_prime_end
                    j += 1
                    return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end  # this exits the global if statement

        else:  # find not exact C' box (up to max 3 substitutions allowed)
            for sub in range(1, int((len_c_prime_box-1)/2 + 1)):  # iterate over 1 to 3 substitutions allowed
                c_prime_motif = regex.findall("((A|G)UGAUGA){s<="+str(sub)+"}", seq[d_prime_end+2:-20], overlapped=True)
                if len(c_prime_motif) >= 1:  # if we have a match, break and keep that match (1 sub privileged over 2 subs)
                    c_prime_motif = c_prime_motif[0][0]  # if multiple C' boxes found, keep the the C' box closest to 5' (i.e. closest to D' box)
                    c_prime_start = seq[d_prime_end+2:-20].index(c_prime_motif) + d_prime_end
                    c_prime_end = c_prime_start + len(c_prime_motif) - 1
                    return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end  # this exits the global if statement
            # If no C' box is found, return NNNNNNN and 0, 0 as C' box sequence, start and end
            c_prime_motif, c_prime_start, c_prime_end = 'NNNNNNN', 0, 0
            return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end

    else: # find not exact D' box (up to max 2 substitutions allowed)
        for sub in range(1, int((len_d_prime_box)/2 + 1)):  # iterate over 1 to 2 substitutions allowed
            d_prime_motif = regex.findall("(CUGA){s<="+str(sub)+"}", middle_seq, overlapped=True)
            if len(d_prime_motif) >= 1:  # if we have a match, break and keep that match (1 sub privileged over 2 subs)
                d_prime_motif = d_prime_motif[0]  # if multiple D' boxes found, keep the the D' box closest to 5' end
                d_prime_start = middle_seq.find(d_prime_motif) + 21
                d_prime_end = d_prime_start + len(d_prime_motif) - 1
                # Find closest exact C' box downstream of found D' box (downstream by at least 2 nt (1 nt between D' and C'))
                if re.search('(A|G)UGAUGA', seq[d_prime_end+2:-20]) is not None:  # find exact C' box
                    k = 1
                    for possible_c_prime in re.finditer('(A|G)UGAUGA', seq[d_prime_end+2:-20]):
                        if k <= 1:  # select first matched group only (closest RUGAUGA to found D' box)
                            c_prime_motif = possible_c_prime.group(0)
                            c_prime_start = possible_c_prime.start() + d_prime_end
                            c_prime_end = possible_c_prime.end() + d_prime_end
                            k += 1
                            return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end  # this exits the global else statement
                else:  # find closest not exact C' box (up to max 3 substitutions allowed)
                    for sub in range(1, int((len_c_prime_box-1)/2 + 1)):  # iterate over 1 to 3 substitutions allowed
                        c_prime_motif = regex.findall("((A|G)UGAUGA){s<="+str(sub)+"}", seq[d_prime_end+2:-20], overlapped=True)
                        if len(c_prime_motif) >= 1:  # if we have a match, break and keep that match (1 sub privileged over 2 subs)
                            c_prime_motif = c_prime_motif[0][0]  # if multiple C' boxes found, keep the the C' box closest to 5' (i.e. closest to D' box)
                            c_prime_start = seq[d_prime_end+2:-20].index(c_prime_motif) + d_prime_end
                            c_prime_end = c_prime_start + len(c_prime_motif) - 1
                            return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end  # this exits the global else statement
                    # If no C' box is found, return NNNN and 0, 0 as C' box sequence, start and end
                    c_prime_motif, c_prime_start, c_prime_end = 'NNNNNNN', 0, 0
                    return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end  # this exits the global else statement
        # If no D' and no C' box are found, return NNNN/NNNNNNN and 0, 0 as D'/C' box sequence, start and end
        d_prime_motif, d_prime_start, d_prime_end = 'NNNN', 0, 0
        c_prime_motif, c_prime_start, c_prime_end = 'NNNNNNN', 0, 0
        return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end


def find_c_prime_d_prime_hamming(seq):
    """ Find best C'/D' pair that minimizes Hamming distance compared to consensus C'/D' motif. """
    # Get the nucleotides between the 20th and last 20th nucleotides
    middle_seq = seq[20:-20]
    len_c_prime_box, len_d_prime_box = 7, 4

    # Find all possible C' boxes and their start/end and compute their Hamming distance compared to consensus motif
    hamming_c_prime, all_c_primes, all_c_primes_start, all_c_primes_end, temp_motif = [], [], [], [], ''
    for sub in range(0, int((len_c_prime_box-1)/2 + 1)):
        c_prime_motif = regex.findall("((A|G)UGAUGA){s<="+str(sub)+"}", middle_seq, overlapped=True)
        if len(c_prime_motif) >= 1:
            for motif in c_prime_motif:
                if motif[0] != temp_motif:  # to avoid repeated motif between 0 and 1 substitution allowed
                    all_c_primes.append(motif[0])
                    hamming_c_prime.append(sub)
                    c_prime_start = middle_seq.index(motif[0]) + 21
                    c_prime_end = c_prime_start + len(motif[0]) - 1
                    all_c_primes_start.append(c_prime_start)
                    all_c_primes_end.append(c_prime_end)                    
                    temp_motif = motif[0]              
    if len(all_c_primes) == 0:  # if no C' box was found, hamming distance is 7, C' motif is NNNNNNN and start and end are 0
        hamming_c_prime, all_c_primes, all_c_primes_start, all_c_primes_end = [7], ['NNNNNNN'], [0], [0]
    
    # Find all possible D' boxes and their start/end and compute their Hamming distance compared to consensus motif
    hamming_d_prime, all_d_primes, all_d_primes_start, all_d_primes_end, temp_motif = [], [], [], [], ''
    for sub in range(0, int((len_d_prime_box)/2 + 1)):
        d_prime_motif = regex.findall("(CUGA){s<="+str(sub)+"}", middle_seq, overlapped=True)
        if len(d_prime_motif) >= 1:
            for motif in d_prime_motif:
                if motif != temp_motif:  # to avoid repeated motif between 0 and 1 substitution allowed
                    all_d_primes.append(motif)
                    hamming_d_prime.append(sub)
                    d_prime_start = middle_seq.index(motif) + 21
                    d_prime_end = d_prime_start + len(motif) - 1
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
            if (all_d_primes_end[i] <= all_c_primes_start[j] + 2) | ((all_d_primes[i] != 'NNNN') & (all_c_primes[j] == 'NNNNNNN')):   
                temp_total_hamming = hamming_d_prime[i] + hamming_c_prime[j]
                if temp_total_hamming < total_hamming:
                    total_hamming = temp_total_hamming
                    d_prime_index = i
                    c_prime_index = j

            # if no D' nor C' box are found
            elif (all_d_primes[i] != 'NNNN') & (all_c_primes[j] == 'NNNNNNN'):
                temp_total_hamming = hamming_d_prime[i] + hamming_c_prime[j]
                d_prime_index, c_prime_index = 0, 0


    # If only one C' box is found and multiple D' are found (i.e. for 3 snoRNAs) but are overlapping, keep the D' box with the 
    # lowest Hamming distance (closest to 5' if multiple D' have the same Hamming) and return NNNNNNN as the C'
    if (d_prime_index == 10000000000) | (c_prime_index == 10000000000):
        all_c_primes, all_c_primes_start, all_c_primes_end, c_prime_index = ['NNNNNNN'], [0], [0], 0
        d_prime_index = hamming_d_prime.index(min(hamming_d_prime))
    c_prime_motif, c_prime_start, c_prime_end = all_c_primes[c_prime_index], all_c_primes_start[c_prime_index], all_c_primes_end[c_prime_index]
    d_prime_motif, d_prime_start, d_prime_end = all_d_primes[d_prime_index], all_d_primes_start[d_prime_index], all_d_primes_end[d_prime_index]
    return c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end



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


def generate_df_prime(fasta):
    """ From a fasta of snoRNA sequences, find a given motif (C' or D') using
        predefined function find_c_prime_d_prime_hamming and output the motif sequence,
        start and end as a df."""
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
                c_prime_motif, c_prime_start, c_prime_end, d_prime_motif, d_prime_start, d_prime_end = find_c_prime_d_prime_hamming(seq)
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


def find_all_boxes(fasta, path):
    """ Find C, D, C' and D' boxes in given fasta using generate_df and concat
        resulting dfs horizontally."""
    df_c = generate_df(fasta, find_c_box, 'C')
    df_d = generate_df(fasta, find_d_box, 'D')
    df_c_prime_d_prime = generate_df_prime(fasta)

    df_final = reduce(lambda  left,right: pd.merge(left,right,on=['gene_id'],
                                            how='outer'),
                                            [df_c, df_d, df_c_prime_d_prime])
    df_final.to_csv(path, sep='\t', index=False)


find_all_boxes(cd_fa, output)
