#!/usr/bin/python3
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
from math import log2
import numpy as np

""" Create logo of C and D boxes from fasta of either expressed or
    not expressed C/D box snoRNAs."""

fasta_c = snakemake.input.c_box_fasta
fasta_d = snakemake.input.d_box_fasta
logo_outputs = [snakemake.output.c_logo, snakemake.output.d_logo]


# Get all box sequences (not sno_id) in a list
# Get the name of group of C/D to redirect figure to correct output

for j, fasta in enumerate([fasta_c, fasta_d]):
    d = {}
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                id = line.replace('>','').replace('\n', '')
            else:
                d[id] = line.replace('\n', '')
    seqs = list(d.values())

    #Get a count and probability matrix to create the logo
    counts_matrix = logomaker.alignment_to_matrix(seqs)
    prob_matrix = logomaker.transform_matrix(counts_matrix, from_type='counts',
                                            to_type='probability')

    # Create logo wo blanks
    rc = {'ytick.labelsize': 32}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    logo = logomaker.Logo(prob_matrix, color_scheme='classic')
    logo.ax.set_ylabel("Frequency", fontsize=35)
    if j == 0:  # for C box
        logo.ax.set_xlabel("Relative position to C box", fontsize=35)
        xticklabels = ["-"+str(i) for i in range(len(prob_matrix), 0, -1)]
    else:  # for D box
        logo.ax.set_xlabel("Relative position to D box", fontsize=35)
        xticklabels = [str(i) for i in range(1, len(prob_matrix) + 1)]
    logo.ax.set_xticks(np.arange(len(prob_matrix)))
    logo.ax.set_xticklabels(xticklabels, fontsize=20)
    plt.savefig(logo_outputs[j], bbox_inches='tight', dpi=600)



