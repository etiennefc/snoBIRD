#!/usr/bin/python3
import shap
import numpy as np
import regex
import re
import pandas as pd
import subprocess as sp
from scipy.signal import find_peaks, savgol_filter


""" Various functions used by SnoBIRD to make its 
    predictions and compute SHAP values"""

def seq2kmer(seq, k):
    """
    Convert original str sequence to kmers (str of all kmers)
    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers

def kmer_score_per_nt(seq, scores, k):
    """
    Get the avg SHAP value per nt across 
    the kmers in which that nt is found.
    """
    seq_scores = []
    L = len(seq)
    firsts = [i for i in range(k-1)]
    lasts = [i for i in range(L-k, L)]
    for i, nt in enumerate(seq):
        if i in firsts:
            sc = scores[0:i+1]
        elif i in lasts:
            sc = scores[-(L-i):]
        else:
            sc = scores[i+1-k:i+1]
        avg_score = np.mean(sc)
        seq_scores.append(avg_score)
    return seq_scores

def batch_generator(sequences, batch_size):
    """
    Create a generator of sequences using a 
    given batch size.
    """
    for i in range(0, len(sequences), batch_size):
        yield sequences[i:i+batch_size]

def shap_batch(text_batch, gene_ids, pipeline):
    """
    Compute SHAP values for a given batch of predicted 
    sequences by the first SnoBIRD model.
    """
    pred_label_dict = {'LABEL_0': 'Other', 'LABEL_1': 'CD_snoRNA'}
    
    # Perform prediction on the batch
    prediction_batch = pipeline([seq2kmer(text, 6) for text in text_batch])
    predicted_labels = [pred['label'] for pred in prediction_batch]
    probabilities = [pred['score'] for pred in prediction_batch]
    
    # Compute SHAP values for the batch
    explainer = shap.Explainer(pipeline, seed=42, output_names=['Other', 
                            'CD_snoRNA'], algorithm='partition', max_evals=500)
    shap_values_batch = explainer([seq2kmer(t, 6) for t in text_batch])
    shap_values_batch = shap_values_batch.values
    
    # Process SHAP values for each sequence in the batch
    avg_shap_per_nt_batch = []
    for i, (shap_values, gene_id, predicted_label, prob) in enumerate(
        zip(shap_values_batch, gene_ids, predicted_labels, probabilities)):
        pos_shap_vals = shap_values[:, 1]
        CLS = pos_shap_vals[0]
        SEP = pos_shap_vals[-1]
        other_tokens = pos_shap_vals[1:-2]
        avg_shap_per_nt = kmer_score_per_nt(text_batch[i], other_tokens, 6)
        avg_shap_per_nt = [gene_id, pred_label_dict[predicted_label], prob, 
                            CLS] + avg_shap_per_nt + [SEP]
        avg_shap_per_nt_batch.append(avg_shap_per_nt)
    
    return avg_shap_per_nt_batch







""" Functions to find the boxes in C/D snoRNAs"""
len_c_box, len_d_box = 7, 4

def find_c_box(seq, middle_index, window_range):
    """ 
    Find exact C box (RTGATGA, where R is A or G), if not present, find C
    box with 1,2 or max 3 substitutions. Return also the start and end
    position of that box as 1-based values with regards to the whole window 
    sequence, not just the subsequence in which we look for the C box. 
    If no C box is found, return a 'NNNNNNN' empty C box and None as start 
    and end of C box.
    """
    shift = middle_index - window_range
    # First, find exact C box (RTGATGA) within snoRNA sequence
    if re.search('(A|G)TGATGA', seq) is not None:  # find exact C box
        i = 1
        for possible_c in re.finditer('(A|G)TGATGA', seq): 
            # select first matched group only (closest RTGATGA to 5' 
            # end of snoRNA)
            if i <= 1: 
                c_motif = possible_c.group(0)
                c_start = possible_c.start() + 1 + shift
                c_end = possible_c.end() + shift
                i += 1
                return c_motif, c_start, c_end  
    else:  # find not exact C box (up to max 3 substitution allowed)
        # iterate over 1 to 3 substitutions allowed
        for sub in range(1, int((len_c_box-1)/2 + 1)):  
            c_motif = regex.findall("((A|G)TGATGA){s<="+str(sub)+"}", seq, 
                        overlapped=True)
            # if we have a match, break and keep that match 
            # (1 sub privileged over 2 subs)
            if len(c_motif) >= 1:  
                # if multiple C boxes found, keep the C box closest to 5' end
                c_motif = c_motif[0][0]  
                c_start = seq.find(c_motif) + 1 + shift
                c_end = c_start + len(c_motif) - 1 
                return c_motif, c_start, c_end  
        # If no C box is found, return NNNNNNN and None, None as C 
        # box sequence, start and end
        c_motif, c_start, c_end = 'NNNNNNN', None, None
        return c_motif, c_start, c_end

def find_d_box(seq, middle_index, window_range):
    """ 
    Find exact D box (CTGA), if not present, find D box with 1 or max 2
    substitutions. Return also the start and end position of that box as 
    1-based values with regards to the whole window sequence, not just 
    the subsequence in which we look for the D box. If no D box is found, 
    return a 'NNNN' empty D box and None as start and end of D box.
    """
    shift = middle_index - window_range
    
    # First, find exact D box (CTGA) within snoRNA sequence
    if re.search('CTGA', seq) is not None:  # find exact D box
        *_, last_possible_d = re.finditer('CTGA', seq)
        # if multiple exact D boxes found, keep the D box closest to 3' end
        d_motif = last_possible_d.group(0)  
        d_start = last_possible_d.start() + 1 + shift
        d_end = last_possible_d.end() + shift
        return d_motif, d_start, d_end
    else:  # find not exact D box (max 50% of substitution allowed (i.e. 2 nt))
        # iterate over 1 to 2 substitutions allowed
        for sub in range(1, int(len_d_box/2 + 1)):  
            d_motif = regex.findall("(CTGA){s<="+str(sub)+"}", seq, 
                        overlapped=True)
            # if we have a match, break and keep that match 
            # (1 sub privileged over 2 subs)
            if len(d_motif) >= 1:  
                # if multiple D boxes found, keep the D box closest to 3' end
                d_motif = d_motif[-1]  
                d_start = seq.rindex(d_motif) + 1 + shift
                d_end = d_start + len(d_motif) - 1 
                return d_motif, d_start, d_end
        # If no D box is found, return NNNN and None, None as D 
        # box sequence, start and end
        d_motif, d_start, d_end = 'NNNN', None, None
        return d_motif, d_start, d_end


def find_perfect_d_box(seq, middle_index, window_range):
    """ 
    Find only exact D box (CTGA). Return also the start and end position of 
    that box as 1-based values with regards to the whole window sequence, not 
    just the subsequence in which we look for the D box.
    """
    shift = middle_index - window_range
    
    # First, find exact D box (CTGA) within snoRNA sequence
    if re.search('CTGA', seq) is not None:  # find exact D box
        *_, last_possible_d = re.finditer('CTGA', seq)
        # if multiple exact D boxes found, keep the D box closest to 3' end
        d_motif = last_possible_d.group(0)  
        d_start = last_possible_d.start() + 1 + shift
        d_end = last_possible_d.end() + shift
    else:
        d_motif, d_start, d_end = 'NNNN', 0, 0
    return d_motif, d_start, d_end


def find_c_prime_d_prime_hamming(seq, c_start, d_end):
    """ 
    Find best C'/D' pair that minimizes Hamming 
    distance compared to consensus C'/D' motif. 
    """
    
    len_c_prime_box, len_d_prime_box = 7, 4
    len_c_box, len_d_box = 7, 4
    middle = len(seq) // 2
    # Get the nucleotides between the middle of the sequence and
    # the start of the D box (-2 nt space)
    c_prime_range_seq = seq[middle:int(d_end - len_d_box -10)]
    # Get the nucleotides between the end of the C box (+2 nt space) and
    # the middle of the sequence
    d_prime_range_seq = seq[int(c_start + len_c_box + 10):middle]

    # Find all possible C' boxes and their start/end and compute their 
    # Hamming distance compared to consensus motif
    hamming_c_prime, all_c_primes, all_c_primes_start = [], [], [] 
    all_c_primes_end, temp_motif =  [], ''
    for sub in range(0, int((len_c_prime_box-1)/2 + 1)):
        c_prime_motif = regex.findall("((A|G)TGATGA){s<="+str(sub)+"}", 
                    c_prime_range_seq, overlapped=True)
        if len(c_prime_motif) >= 1:
            for motif in c_prime_motif:
                # To avoid repeated motif between 0 and 1 substitution allowed
                if motif[0] != temp_motif:  
                    all_c_primes.append(motif[0])
                    hamming_c_prime.append(sub)
                    c_prime_start = c_prime_range_seq.index(
                        motif[0]) + middle + 1  # +1 to become 1-based 
                    c_prime_end = c_prime_start + len_c_prime_box - 1
                    all_c_primes_start.append(c_prime_start)
                    all_c_primes_end.append(c_prime_end)                    
                    temp_motif = motif[0]              
    if len(all_c_primes) == 0:  # if no C' box was found, hamming dist is 7, 
                                # C' motif is NNNNNNN and start and end are 0
        hamming_c_prime, all_c_primes, all_c_primes_start, all_c_primes_end = (
                        [7], ['NNNNNNN'], [0], [0]
                        )
    
    # Find all possible D' boxes and their start/end and compute their Hamming 
    # distance compared to consensus motif
    hamming_d_prime, all_d_primes, all_d_primes_start = [], [], [] 
    all_d_primes_end, temp_motif = [], ''
    for sub in range(0, int((len_d_prime_box)/2 + 1)):
        d_prime_motif = regex.findall("(CTGA){s<="+str(sub)+"}", 
                        d_prime_range_seq, overlapped=True)
        if len(d_prime_motif) >= 1:
            for motif in d_prime_motif:
                # To avoid repeated motif between 0 and 1 substitution allowed
                if motif != temp_motif:  
                    all_d_primes.append(motif)
                    hamming_d_prime.append(sub)
                    d_prime_start = d_prime_range_seq.index(
                                    motif) + c_start + len_c_box + 10 + 1
                    # +2 for the 2nt buffer after C box; +1 to become 1-based
                    d_prime_end = d_prime_start + len_d_prime_box - 1
                    all_d_primes_start.append(d_prime_start)
                    all_d_primes_end.append(d_prime_end)
                    temp_motif = motif                
    if len(all_d_primes) == 0:  # if no D' box was found, hamming dist is 4, 
                                # D' motif is NNNN and start and end are 0
        hamming_d_prime, all_d_primes, all_d_primes_start, all_d_primes_end = (
            [4], ['NNNN'], [0], [0]
            )
        
    # Find all possible D'-C' pairs where C' is downstream of D' by at 
    # least 2 nt (i.e. at least 1 nt between D' and C' boxes) 
    # and return the best pair according to the lowest total Hamming distance 
    # (if two pairs have the same Hamming distance, the closest one to 5' of 
    # snoRNA is chosen)
    total_hamming, d_prime_index, c_prime_index = (
        10000000000, 10000000000, 10000000000  # these are dummy high numbers 
    )
    for i, d_prime in enumerate(all_d_primes):
        for j, c_prime in enumerate(all_c_primes):
            # C' downstream of D' by at least 2 nt or if D' is found but not C' 
            # (the case where no D' is found but C' is found is also included: 
            # d_prime_end = 0 is smaller than any existing c_prime_start)
            if (
                all_d_primes_end[i] +2 <= all_c_primes_start[j]) | \
                ((all_d_primes[i] != 'NNNN') & (all_c_primes[j] == 'NNNNNNN')):   
                temp_total_hamming = hamming_d_prime[i] + hamming_c_prime[j]
                if temp_total_hamming < total_hamming:
                    total_hamming = temp_total_hamming
                    d_prime_index = i
                    c_prime_index = j

            # if no D' nor C' box are found
            elif (all_d_primes[i] != 'NNNN') & (all_c_primes[j] == 'NNNNNNN'):
                temp_total_hamming = hamming_d_prime[i] + hamming_c_prime[j]
                d_prime_index, c_prime_index = 0, 0

    # If only one C' box is found and multiple D' are found  but are 
    # overlapping, keep the D' box with the lowest Hamming distance (closest to
    # 5' if multiple D' have the same Hamming) and return NNNNNNN as the C'
    if (d_prime_index == 10000000000) | (c_prime_index == 10000000000):
        all_c_primes, all_c_primes_start, all_c_primes_end, c_prime_index = (
            ['NNNNNNN'], [0], [0], 0
            )
        d_prime_index = hamming_d_prime.index(min(hamming_d_prime))
    c_prime_motif, c_prime_start, c_prime_end = (all_c_primes[c_prime_index], 
        all_c_primes_start[c_prime_index], all_c_primes_end[c_prime_index]
        )
    d_prime_motif, d_prime_start, d_prime_end = (all_d_primes[d_prime_index], 
        all_d_primes_start[d_prime_index], all_d_primes_end[d_prime_index]
        )
    return (c_prime_motif, c_prime_start, c_prime_end,  
            d_prime_motif, d_prime_start, d_prime_end
            )


def hamming(motif, consensus='RTGATGA'):
    """ 
    Compute Hamming distance with regards to consensus motif.
    """
    score = 0
    for i, s in enumerate(motif):
        if i == 0:
            if consensus == 'RTGATGA':
                if s not in ['A', 'G']:
                    score +=1
            else:  # for D box
                if s != consensus[i]:
                    score +=1
        else:
            if s != consensus[i]:
                score +=1
    return score


def symm_region(position, sequence, box_len, max_length=194):
    """ 
    From a subsequence position in a given predicted window, 
    find its symmetrical region with regards to the middle 
    of the window. If the subsequence contains a C box, 
    then its symmetrical region should contain a D box and 
    vice-versa if the subsequence contains a D box.
    """
    middle = max_length / 2
    symmetrical_pos = max_length - position - 1
    symmetrical_region = sequence[
        symmetrical_pos-box_len:symmetrical_pos+box_len+1
    ]
    return symmetrical_region, symmetrical_pos


def find_boxes_SHAP(row, fixed_length, shap_cols, C_range, D_range, 
                    len_c_box=len_c_box, len_d_box=len_d_box):
    """
    Based on highest SHAP values out of the first model of SnoBIRD, define 
    best C and D boxes per snoRNA.
    """
    ext_seq = row[f'extended_{fixed_length}nt_sequence']
    gene_id = row['gene_id']
    all_shap = list(row[shap_cols])
    cls_ = [row['CLS']]
    sep_ = [row['SEP']]
    prob = row['probability_CD']  # 1st model prediction prob

    # Use Savitzky-Golay filter to smoothen the SHAP value 
    # distribution into clear peaks
    smoothed = savgol_filter(np.array(cls_ + all_shap + sep_), 11, 2)
    peak_indexes, prop_dict = find_peaks(smoothed, threshold=0.00001, 
                                        height=0.005)
    peak_heights = prop_dict['peak_heights']

    # Find potential C peaks and potential D peaks
    # Sort these peaks by descending value of peak value (most important first)
    potential_peaks = {
        f'peak_{peak}': peak_heights[i]
        for i, peak in enumerate(peak_indexes)
        if ((peak - len_c_box in C_range) & (peak + len_c_box in C_range)) |
           ((peak - len_d_box in D_range) & (peak + len_d_box in D_range))
    }
    potential_peaks = dict(sorted(potential_peaks.items(), 
                        key=lambda item: item[1], reverse=True))
    potential_C = {k: v for k, v in potential_peaks.items() if int(
                                        k.split('_')[1]) in C_range}
    potential_D = {k: v for k, v in potential_peaks.items() if int(
                                        k.split('_')[1]) in D_range}
    
    best_c_motif, best_c_start, best_c_end, best_score_c = (
                    'NNNNNNN', None, None, 7
                    )
    best_d_motif, best_d_start, best_d_end, best_score_d = (
                    'NNNN', None, None, 4
                    )
    best_box_score = best_score_c + best_score_d

    # Go over all peaks, find the best C with its corresponding D or 
    # best D with its corresponding C box,
    # corresponding meaning being in a symmetrical range
    for pot_peak, pot_val in potential_peaks.items():
        peak_index = int(pot_peak.split('_')[1])
        
        if pot_peak in potential_C.keys():
            c_region = ext_seq[peak_index - len_c_box:peak_index + len_c_box]
            c_motif, c_start, c_end = find_c_box(
                                        c_region, peak_index, len_c_box)
            score_c = hamming(c_motif)
            # We use len_c_box in symm_region even to find the D box in order 
            # to have a range length comparable to the C box range length 
            # (otherwise not wide enough to find D boxes)
            d_region, d_region_pos = symm_region(
                                        peak_index, ext_seq, len_c_box)
            d_motif, d_start, d_end = find_d_box(
                                        d_region, d_region_pos, len_d_box + 3)
                    # + 3 to remove the difference between C & D motif lengths
            score_d = hamming(d_motif, consensus='CTGA')
            box_score = score_c + score_d
            if box_score < best_box_score:
                # Update best box C and D
                best_c_motif, best_c_start, best_c_end, best_score_c = (
                            c_motif, c_start, c_end, score_c
                            )
                best_d_motif, best_d_start, best_d_end, best_score_d = (
                            d_motif, d_start, d_end, score_d
                            )       
                best_box_score = box_score
        elif pot_peak in potential_D.keys():
            d_region = ext_seq[peak_index - len_d_box:peak_index + len_d_box]
            d_motif, d_start, d_end = find_d_box(
                                        d_region, peak_index, len_d_box)
            score_d = hamming(d_motif, consensus='CTGA')
            c_region, c_region_pos = symm_region(
                                        peak_index, ext_seq, len_c_box)
            c_motif, c_start, c_end = find_c_box(
                                        c_region, c_region_pos, len_c_box)
            score_c = hamming(c_motif)
            box_score = score_c + score_d
            if box_score < best_box_score:
                # Update best box C and D
                best_c_motif, best_c_start, best_c_end, best_score_c = (
                            c_motif, c_start, c_end, score_c
                            )
                best_d_motif, best_d_start, best_d_end, best_score_d = (
                            d_motif, d_start, d_end, score_d
                            )
                best_box_score = box_score

    # a C and D box were found
    if (best_c_motif != 'NNNNNNN') & (best_d_motif != 'NNNN'): 
        if best_d_motif != 'CTGA':
            # Try to find an exact D box in the 25 nt following the previously
            # found D box, otherwise a lot of real D boxes are missed
            max_range = min(fixed_length - 20, best_d_end + 25)
            d_region = ext_seq[best_d_start:max_range]
            window_ = max_range - best_d_start
            d_motif, d_start, d_end = find_perfect_d_box(
                                        d_region, max_range, window_)
            if d_motif == 'CTGA':
                score_d = hamming(d_motif, consensus='CTGA')
                best_d_motif, best_d_start, best_d_end, best_score_d = (
                            d_motif, d_start, d_end, score_d)
                best_box_score = best_score_c + best_score_d

        return (best_c_motif, best_c_start, best_c_end, best_d_motif, 
                best_d_start, best_d_end, gene_id, best_score_c, best_score_d, 
                best_box_score, prob)

    # If either or both C and D boxes were not found, use more lenient range 
    # values to find the boxes (don't count previously found peaks with 
    # stringent range values). Usually the missed peaks are C peaks closer to 
    # the middle of the snoRNA sequence. More lenient because we don't add or 
    # substract box length to the peak to see if it fits in the given range
    else:
        lenient_peaks = {
            f'peak_{peak}': peak_heights[i]
            for i, peak in enumerate(peak_indexes)
            if ((peak in C_range) | (peak in D_range)) &
               ((f'peak_{peak}' not in potential_C.keys()) & 
               (f'peak_{peak}' not in potential_D.keys()))
        }
        lenient_peaks = dict(sorted(lenient_peaks.items(), 
                            key=lambda item: item[1], reverse=True))
        lenient_C = {k: v for k, v in lenient_peaks.items() if int(
                                                k.split('_')[1]) in C_range}
        lenient_D = {k: v for k, v in lenient_peaks.items() if int(
                                                k.split('_')[1]) in D_range}

        if len(lenient_peaks.keys()) == 0:  # no new peaks were found
            # Keep the previous best C and D boxes if they were found, but 
            # first try to find an exact D box in the 25 nt following the
            # previously found D box, otherwise lots of real D boxes are missed
            if (best_d_motif != 'CTGA') & (best_d_motif != 'NNNN'):
                max_range = min(fixed_length - 20, best_d_end + 25)
                d_region = ext_seq[best_d_start:max_range]
                window_ = max_range - best_d_start
                d_motif, d_start, d_end = find_perfect_d_box(
                                        d_region, max_range, window_)
                if d_motif == 'CTGA':
                    score_d = hamming(d_motif, consensus='CTGA')
                    best_d_motif, best_d_start, best_d_end, best_score_d = (
                                d_motif, d_start, d_end, score_d)
                    best_box_score = best_score_c + best_score_d
            return (best_c_motif, best_c_start, best_c_end, best_d_motif, 
                    best_d_start, best_d_end, gene_id, best_score_c, 
                    best_score_d, best_box_score, prob)
        else:   # iterate over new lenient peaks to find the best 
                # potential C and D boxes
            for pot_peak, pot_val in lenient_peaks.items():
                peak_index = int(pot_peak.split('_')[1])

                if pot_peak in lenient_C.keys():
                    c_region = ext_seq[
                                peak_index - len_c_box:peak_index + len_c_box]
                    c_motif, c_start, c_end = find_c_box(
                                                c_region, peak_index, len_c_box
                                                )
                    score_c = hamming(c_motif)
                    # We use len_c_box in symm_region even to find the D box in 
                    # order to have a range length comparable to the C box 
                    # range length (otherwise not wide enough to find D boxes)
                    d_region, d_region_pos = symm_region(
                                                peak_index, ext_seq, len_c_box)
                    d_motif, d_start, d_end = find_d_box(
                                        d_region, d_region_pos, len_d_box + 3)
                    # + 3 to remove the difference between C & D motif lengths
                    score_d = hamming(d_motif, consensus='CTGA')
                    box_score = score_c + score_d
                    if box_score < best_box_score:
                        # Update best box C and D
                        best_c_motif, best_c_start = c_motif, c_start
                        best_c_end, best_score_c = c_end, score_c
                        best_d_motif, best_d_start = d_motif, d_start
                        best_d_end, best_score_d = d_end, score_d
                        best_box_score = box_score
                elif pot_peak in lenient_D.keys():
                    d_region = ext_seq[
                                peak_index - len_d_box:peak_index + len_d_box]
                    d_motif, d_start, d_end = find_d_box(
                                                d_region, peak_index, len_d_box
                                                )
                    score_d = hamming(d_motif, consensus='CTGA')
                    c_region, c_region_pos = symm_region(
                                                peak_index, ext_seq, len_c_box)
                    c_motif, c_start, c_end = find_c_box(
                                            c_region, c_region_pos, len_c_box)
                    score_c = hamming(c_motif)
                    box_score = score_c + score_d
                    if box_score < best_box_score:
                        # Update best box C and D
                        best_c_motif, best_c_start = c_motif, c_start
                        best_c_end, best_score_c = c_end, score_c
                        best_d_motif, best_d_start = d_motif, d_start
                        best_d_end, best_score_d = d_end, score_d
                        best_box_score = box_score

            # Try as a last resort to find an exact D box in the 25 nt 
            # following the previously found D box, otherwise lots of real 
            # D boxes are missed
            if (best_d_motif != 'CTGA') & (best_d_motif != 'NNNN'):
                max_range = min(fixed_length - 20, best_d_end + 25)
                d_region = ext_seq[best_d_start:max_range]
                window_ = max_range - best_d_start
                d_motif, d_start, d_end = find_perfect_d_box(
                                        d_region, max_range, window_)
                if d_motif == 'CTGA':
                    score_d = hamming(d_motif, consensus='CTGA')
                    best_d_motif, best_d_start, best_d_end, best_score_d = (
                                d_motif, d_start, d_end, score_d)
                    best_box_score = best_score_c + best_score_d
            return (best_c_motif, best_c_start, best_c_end, best_d_motif, 
                    best_d_start, best_d_end, gene_id, best_score_c, 
                    best_score_d, best_box_score, prob)


def find_all_boxes(df, fixed_length, shap_cols, C_range, D_range, flanking_nt):
    """
    Find C, D, C' and D' boxes using find_boxes_SHAP and 
    find_c_prime_d_prime_hamming.
    """
    # Find C and D boxes based on SHAP
    results = df.apply(lambda row: find_boxes_SHAP(row, fixed_length, 
                                        shap_cols, C_range, D_range), axis=1)
    results_df = pd.DataFrame(results.tolist(), columns=['C_MOTIF', 'C_START', 
                'C_END', 'D_MOTIF', 'D_START', 'D_END', 'gene_id', 'score_c', 
                'score_d', 'box_score', 'probability_CD'])

    # Fix C_start/ D_end for examples where one of C or D boxes couldn't be 
    # found by getting symmetrical region to the only found box in order to
    # delimit sno start/end
    results_df.loc[(results_df['C_MOTIF'] == 'NNNNNNN') & 
                    (results_df['D_MOTIF'] != 'NNNN') & 
                    (results_df['C_START'].isnull()), 
                    'C_START'] = fixed_length - results_df['D_END'] - 1
    results_df.loc[(results_df['C_MOTIF'] != 'NNNNNNN') & 
                    (results_df['D_MOTIF'] == 'NNNN') & 
                    (results_df['D_START'].isnull()), 
                    'D_END'] = fixed_length - results_df['C_START'] - 1

    # Fix C_start/D_end for examples where no C and D box could be found
    # C box should start at least after flanking nt
    results_df['C_START'] = results_df['C_START'].fillna(flanking_nt)  
    # D box should end at most before flanking nt
    results_df['D_END'] = results_df['D_END'].fillna(
                                                fixed_length - flanking_nt)

    # Add sequence to results_df
    results_df = results_df.merge(df[['gene_id', 
                    f'extended_{fixed_length}nt_sequence']], how='left', 
                    on='gene_id')

    # Find C' and D' boxes
    prime_cols = ['C_PRIME_MOTIF', 'C_PRIME_START', 'C_PRIME_END', 
                'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END']
    results_df[prime_cols] = results_df.apply(
            lambda row: find_c_prime_d_prime_hamming(
                        row[f'extended_{fixed_length}nt_sequence'], 
                        row['C_START'], row['D_END']), axis=1, 
                        result_type='expand')
    results_df['score_c_prime'] = results_df['C_PRIME_MOTIF'].apply(
                                                        lambda x: hamming(x))
    results_df['score_d_prime'] = results_df['D_PRIME_MOTIF'].apply(
                                        lambda x: hamming(x, consensus='CTGA'))
    results_df['box_score'] = (results_df.score_c + results_df.score_d + 
                        results_df.score_c_prime + results_df.score_d_prime)
    return results_df


def correct_box_pos(row, motifs):
    """
    Convert box position that are initially based on the extended window length
    to the actual snoRNA length within that window.
    """
    strand = row['strand']
    window_start = int(row['start_window'])
    window_end = int(row['end_window'])  # end of window as genomic location
    start = int(row['start'])  # start of snoRNA as genomic location
    end = int(row['end'])
    correction = 1
    for motif in motifs:
        if (row[f'{motif}_MOTIF'] == 'NNNNNNN') | (
            row[f'{motif}_MOTIF'] == 'NNNN'):
            row[f'{motif}_START'] = 0
            row[f'{motif}_END'] = 0
        elif strand == '+':
            row[f'{motif}_START'] = window_start + row[
                                    f'{motif}_START'] - start + correction
            row[f'{motif}_END'] = window_start + row[
                                    f'{motif}_END'] - start + correction
        else: # - strand
            row[f'{motif}_START'] = row[f'{motif}_START'] - (window_end - end) 
            row[f'{motif}_END'] = row[f'{motif}_END'] - (window_end - end) 

    return row


def get_seq(row, fixed_length=194, extended=False, extension=15):
    """
    Get snoRNA predicted sequence based on C and D boxes positions.
    """
    seq = row[f'extended_{fixed_length}nt_sequence']
    C_start = int(row['C_START'])
    D_end = int(row['D_END'])
    s = seq[C_start - 5 -1: D_end+5]

    if extended == True:
        s = seq[max(0, C_start - 5 -1 - extension): min(
                                                D_end+5+extension, len(seq))]

    return s


def get_sno_location(row, fixed_length=194):
    """
    Get snoRNA genomic location based on C and D boxes positions.
    """
    C_start = int(row['C_START'])
    D_end = int(row['D_END'])
    chrom, strand = row['chr_window'], row['strand_window']
    start, end = int(row['start_window']), int(row['end_window'])

    if strand == '+':
        # + 5 before the box 
        new_start = start + C_start - 5 
        # + 5 after the box 
        new_end = start + D_end + 5 
        new_cols = {'chr': chrom, 'start': new_start, 
                    'end': new_end, 'strand': strand}
    
    # Negative strand needs to be corrected because C_start and D_end are 
    # always from 5' to 3' in the sense of the strand, whereas the start and 
    # end of the extended window are based on the + strand even if the gene 
    # is on the - strand
    else:  
        # +1 to account for the 1-based D_end and C_start 
        new_start = start + (fixed_length - D_end) - 5 + 1
        new_end = end - C_start + 5 + 1
        new_cols = {'chr': chrom, 'start': new_start, 
                    'end': new_end, 'strand': strand}

    return pd.Series(new_cols)





""" Functions to deal with window filtering of predictions."""
def merge_rows(row1, row2):
    """
    Function to merge two rows if they overlap reciprocally by at least 50%.
    """
    overlap_start = max(row1['start'], row2['start'])
    overlap_end = min(row1['end'], row2['end'])
    overlap_len = max(0, overlap_end - overlap_start)
    max_len = max(row1['len'], row2['len'])
    
    # test only for the longest block (the shortest will be also true 
    # if the longest overlaps at least at 50%)
    reciprocal_overlap = (overlap_len / max_len) >= 0.5  
    
    if (row1['chrom'] == row2['chrom']) & (
        row1['strand'] == row2['strand']) & (reciprocal_overlap == True):
        merged_row = {
            'chrom': row1['chrom'],
            'start': row1['start'],
            'end': row2['end'],
            'name': str(row1['name']) + '_' + str(row2['name']),
            'score': (row1.score + row2.score) / 2,
            'strand': row1['strand'], 
            'len': row2['end'] - row1['start'] + 1
        }
        return merged_row
    else:
        return None


def centered_window(row, fixed_length):
    """
    From a large window, find the centered window of smaller fixed_length.
    """
    interval_length = row['end'] - row['start'] + 1
    midpoint = row['start'] + (interval_length // 2)
    
    # For odd_interval length, it shifts 1 nt more to the left 
    # than to the right
    start = midpoint - (fixed_length // 2) - 1
    end = midpoint + (fixed_length // 2) - 1

    return pd.Series([start, end])





""" Functions to get the feature values of snoRNAs (snoRNA stability, 
    terminal stability, etc.) and to filter based on these."""

def sno_mfe(df, seq_col, len_col):
    """
    Compute the Minimal Free Energy (MFE) of the snoRNA predicted secondary 
    structure using RNAFold, and normalize that mfe by the snoRNA length.
    """
    seq_dict = dict(zip(df.gene_id, df[seq_col]))
    
    # Create fasta required for RNAfold 
    with open('temp_sno_mfe.fa', 'w') as f:
        for id_, predicted_seq in seq_dict.items():
            f.write(f'>{id_}\n{predicted_seq}\n')

    # Run RNAfold
    sp.call("sed 's/>/>SNOBIRD_/g' temp_sno_mfe.fa > SNOBIRD_rna_fold_mfe.fa", 
            shell=True)
    sp.call(f'RNAfold --infile=SNOBIRD_rna_fold_mfe.fa --outfile=SNOBIRD.mfe', 
            shell=True)
    sp.call("sed -i 's/SNOBIRD_//g' SNOBIRD.mfe", shell=True)
    
    # Get MFE in dict
    mfe_dict = {}
    with open('SNOBIRD.mfe', 'r') as f:
        for line in f:
            if line.startswith('>'):
                gene_id = line.strip('>\n')
            elif line.startswith('.') | line.startswith('('):
                mfe = float(line.split(' ', maxsplit=1)[1].strip(' \n()'))
                mfe_dict[gene_id] = mfe

    # Create cols for structure stability and normalized MFE by the 
    #  predicted length
    df['sno_stability'] = df['gene_id'].map(mfe_dict)
    df['normalized_sno_stability'] = df['sno_stability'] / df[len_col]
    df = df.drop(columns=['sno_stability'])

    sp.call('rm temp_sno_mfe.fa SNOBIRD_rna_fold_mfe.fa *.ps SNOBIRD.mfe', 
            shell=True)
    return df


def flanking_nt(box_location, sequence, motif='C', flank_len=20):
    """
    Get flanking nt surrounding before/after the C/D boxes to compute 
    afterwards the terminal stem stability.
    """
    # (15 external stem nt + 5 nt internal to the snoRNA)
    box_location = int(box_location)
    if motif == 'C':
        seq_ = sequence[box_location - flank_len - 1:box_location-1]
    if motif == 'D':
        seq_ = sequence[box_location: box_location + flank_len]
    return seq_


def terminal_stem(df, seq_col, extended_col=False, extended_col_name=None):
    """
    Compute the Minimal Free Energy (MFE) of the terminal stem of snoRNAs 
    using RNACofold, also get the terminal stem score and create the 
    terminal_combined metrics.
    """
    # Create flanking seq columns
    df['left_20nt_flanking'] = df.apply(lambda row: flanking_nt(
                                    row['C_START'], row[seq_col]), axis=1)
    df['right_20nt_flanking'] = df.apply(lambda row: flanking_nt(row['D_END'], 
                                    row[seq_col], motif='D'), axis=1)
    if extended_col == True:
        df['left_20nt_flanking'] = df[extended_col_name].apply(
                                lambda x: x[0:20])
        df['right_20nt_flanking'] = df[extended_col_name].apply(
                                lambda x: x[-20:])

    flanking_dict = dict(zip(df.gene_id, zip(
                        df.left_20nt_flanking, df.right_20nt_flanking)))
    
    # Create fasta required for RNAcofold 
    # (reversed left seq + '&' + reversed right seq in fasta format)
    with open('SNOBIRD_terminal_stem.fa', 'w') as f:
        for id_, flanking_nts in flanking_dict.items():
            f.write(f'>{id_}\n')
            # Reverse both flanking seq so that it is 
            # correctly represented in the RNAcofold plot
            reverse_left = flanking_nts[0][::-1]
            reverse_right = flanking_nts[1][::-1]
            f.write(f'{reverse_left}&{reverse_right}\n')

    # Run RNAcofold
    sp.call(f'RNAcofold < SNOBIRD_terminal_stem.fa > SNOBIRD_terminal.mfe', 
            shell=True)
    


    # Get terminal_stem stability and length score
    # The terminal stem length score is defined as 
    # score = paired_nt - intramol_paired_nt - gap_number
    terminal_mfe_dict = {}
    terminal_stem_length_dict = {}
    with open('SNOBIRD_terminal.mfe', 'r') as f:
        for line in f:
            if line.startswith('>'):
                gene_id = line.strip('>\n')
            elif line.startswith('.') | line.startswith('('):
                terminal_mfe = float(line.split(' ', 
                                maxsplit=1)[1].strip(' \n()'))
                terminal_mfe_dict[gene_id] = terminal_mfe

                # From the dot bracket, extract the number of paired nt
                dot_bracket = line.split(' ', maxsplit=1)[0]
                # We select only the left flanking region (20 nt)
                dot_bracket = dot_bracket[0:20].replace('\n', '')  
                paired_base = dot_bracket.count('(')
                # These following ')' are intramolecular paired nt
                # (i.e. nt pair within left sequence only)
                intra_paired = dot_bracket.count(')')  
                # This finds all overlapping (and non-overlapping) gaps of 
                # 1 to 19 nt inside the left flanking region
                gaps = re.findall(r'(?=(\(\.{1,19}\())', dot_bracket)
                number_gaps = ''.join(gaps)  # join all gaps in one string
                # Count the number of nt gaps in sequence
                number_gaps = len(re.findall('\.', number_gaps))  

                stem_length_score = paired_base - intra_paired - number_gaps
                if stem_length_score < 0:
                    stem_length_score = 0
                terminal_stem_length_dict[gene_id] = stem_length_score

    df['terminal_stem_stability'] = df['gene_id'].map(terminal_mfe_dict)
    df['terminal_stem_length_score'] = df['gene_id'].map(
                                        terminal_stem_length_dict)
    df['terminal_combined'] = (
            df.terminal_stem_stability * df.terminal_stem_length_score
            )
    df = df.drop(columns=['terminal_stem_stability', 
                'terminal_stem_length_score', 'left_20nt_flanking', 
                'right_20nt_flanking'])
    sp.call('rm SNOBIRD_terminal* *.ps', shell=True)
    return df


def feature_filters(df, prob_column, old_pred_column, new_pred_column, 
                    terminal_combined = -25, sno_mfe = -0.2, box_score = 5, 
                    score_c = 2, score_d = 1, prob_thresh = 0.999):
    """
    Based on the snoRNA features and probability of prediction of the second 
    model, filter the predictions so that we maximize the precision on the 
    snoRNA pseudogene class.
    """
    # Change prediction for the actual label instead of number
    dictio = {0: 'CD_snoRNA_pseudogene', 1: 'expressed_CD_snoRNA'}
    df[new_pred_column] = df[old_pred_column].replace(dictio)

    # Apply 1st filter post-snoBIRD's 2nd model based on prediction probability
    # (minor change to prediction if pred probability > 0.999 (high prob); 
    # major changes if pred probability is <= 0.999 (low prob))
    minor_change = df[df[prob_column] > prob_thresh]
    major_change = df[df[prob_column] <= prob_thresh]
    changed_id = []

    # FIRST FILTER: on prediction with low probability (<=0.999)
    # If C box has >=2 mutations and D box has >=1 mutation --> snoRNA_pseudo
    major_change.loc[(~major_change['gene_id'].isin(changed_id)) & (
                    major_change['score_c'] >= score_c) & (
                    major_change['score_d'] >= score_d), 
                    new_pred_column] = 'CD_snoRNA_pseudogene'
    changed_id.extend(list(major_change[(~major_change['gene_id'].isin(
                            changed_id)) & (
                            major_change['score_c'] >= score_c) & (
                            major_change['score_d'] >= score_d)]['gene_id']))
    
    # SECOND FILTER: on prediction with low probability (<=0.999)
    # If terminal combined < -25 --> expressed_CD_snoRNA
    major_change.loc[(~major_change['gene_id'].isin(changed_id)) & (
                    major_change['terminal_combined'] <= terminal_combined), 
                    new_pred_column] = 'expressed_CD_snoRNA'
    changed_id.extend(list(major_change[(~major_change['gene_id'].isin(
                    changed_id)) & (
                    major_change['terminal_combined'] <= terminal_combined)][
                                                                'gene_id']))
    
    # THIRD FILTER: on prediction with low probability (<=0.999)
    # If <= 2 features are "favorable" for expression --> CD_snoRNA_pseudogene
    major_change.loc[(~major_change['gene_id'].isin(changed_id)) & (((
            major_change['terminal_combined'] <= terminal_combined).astype(
                int) + 
            (major_change['normalized_sno_stability'] <= sno_mfe).astype(int) + 
            (major_change['box_score'] <= box_score).astype(int) + (
            major_change['score_c'] <= score_c).astype(int) + (
            major_change['score_d'] < score_d).astype(int)) <= 2), 
            new_pred_column] = 'CD_snoRNA_pseudogene'
    changed_id.extend(list(major_change.loc[(~major_change['gene_id'].isin(
                    changed_id)) & (((major_change[
                    'terminal_combined'] <= terminal_combined).astype(int) + 
                (major_change['normalized_sno_stability'] <= sno_mfe).astype(
                    int) + 
                (major_change['box_score'] <= box_score).astype(int) + (
                major_change['score_c'] <= score_c).astype(int) + (
                major_change['score_d'] < score_d).astype(int)) <= 2)][
                                                                'gene_id']))
    
    # FOURTH FILTER: on prediction with high probability (>0.999) 
    # (only on CD_snoRNA_pseudogene predictions)
    # If >= 4/5 features are "favorable" for expression --> expressed_CD_snoRNA
    minor_change['last_predicted_label'] = minor_change[new_pred_column]
    minor_change.loc[(
                minor_change.last_predicted_label == 'CD_snoRNA_pseudogene') & 
                (((
                minor_change['terminal_combined'] <= terminal_combined).astype(
                    int) + 
                (minor_change['normalized_sno_stability'] <= sno_mfe).astype(
                    int) + 
                (minor_change['box_score'] <= box_score).astype(int) + (
                minor_change['score_c'] <= score_c).astype(int) + (
                minor_change['score_d'] < score_d).astype(int)) >= 4), 
                new_pred_column] = 'expressed_CD_snoRNA'
    

    concat_df = pd.concat([minor_change, major_change])

    # LAST FILTER: if NNNNNN as C or NNNN as D box, 
    # automatically assign as CD_snoRNA_pseudogene
    concat_df.loc[(concat_df['C_MOTIF'] == 'NNNNNNN') | (
        concat_df['D_MOTIF'] == 'NNNN'), 
        new_pred_column] = 'CD_snoRNA_pseudogene'
    
    concat_df = concat_df.sort_values(['chr', 'start', 'strand'])


    return concat_df





