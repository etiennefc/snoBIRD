#!/usr/bin/python3
import pandas as pd
import subprocess as sp 
import re

fixed_length = snakemake.wildcards.fixed_length
positives_len = snakemake.input.positives_len_fa
negatives_len = snakemake.input.negatives_len_fa
output_positives_fa = snakemake.output.positives_fa
output_negatives_fa = snakemake.output.negatives_fa
output_positives_tsv = snakemake.output.positives_tsv
output_negatives_tsv = snakemake.output.negatives_tsv


# Load dfs into dicts
positives_df = pd.read_csv(snakemake.input.positives_windows, sep='\t')
negatives_df = pd.concat([pd.read_csv(path, sep='\t') for path in snakemake.input.negatives_windows])
pos_dict = dict(zip(positives_df.gene_id, positives_df[f'extended_{fixed_length}nt_sequence']))
neg_dict = dict(zip(negatives_df.gene_id, negatives_df[f'extended_{fixed_length}nt_sequence']))


# Get the terminal stem sequence of positives
with open(positives_len, 'r') as f, open('terminal_temp_pos.fa', 'w') as f2:
    for line in f:
        if line.startswith('>'):
            id_ = line.replace('>', '').replace('\n', '')
        elif re.match(r"^A|T|C|G|U|N", line):
            line = line.replace('\n', '')
            whole_window = pos_dict[id_].replace('T', 'U')  # convert to RNA
            sno_start = re.search(line, whole_window).start()
            if sno_start < 15:                                # if the sno start is within the first 10-14 nt of the window
                left_external_nt = whole_window[0:sno_start]  # i.e. if the C box is within the first 15-19 nt of the window
                left_internal_nt = whole_window[sno_start:sno_start+5+(15-sno_start)]
                line_left = left_external_nt + left_internal_nt
            else:
                left_internal_5nt = line[0:5]
                left_external_15nt = whole_window[sno_start-15:sno_start]
                line_left = left_external_15nt + left_internal_5nt

            sno_end = re.search(line, whole_window).end()
            if sno_end >= (int(fixed_length)-15):                # if the sno end is within the last 14-10 nt of the window
                right_external_nt = whole_window[sno_end:]  # i.e. if the D box is within the last 19-15 nt of the window
                right_internal_nt = whole_window[sno_end-5-(15-len(right_external_nt)):sno_end]
                line_right = right_internal_nt + right_external_nt
            else:
                right_internal_5nt = whole_window[sno_end-5:sno_end]
                right_external_15nt = whole_window[sno_end:sno_end+15]
                line_right = right_internal_5nt + right_external_15nt

            if len(line_right) != 20:
                print(id_, len(line_left), len(line_right))
            # Reverse the order of both flanking sequences with [::-1]
            co_seq = line_left[::-1] + "&" + line_right[::-1]  # this order is how RNAcofold will accurately try to see the
                                                # best base pairing between the left and right extended flanking regions of snoRNAs
            f2.write(f'>{id_}\n{co_seq}\n')


# Get the terminal stem sequence of negatives
with open(negatives_len, 'r') as f, open('terminal_temp_neg.fa', 'w') as f2:
    for line in f:
        if line.startswith('>'):
            id_ = line.replace('>', '').replace('\n', '')
        elif re.match(r"^A|T|C|G|U|N", line):
            line = line.replace('\n', '')
            whole_window = neg_dict[id_].replace('T', 'U')  # convert to RNA
            sno_start = re.search(line, whole_window).start()
            if sno_start < 15:                                # if the sno start is within the first 10-14 nt of the window
                left_external_nt = whole_window[0:sno_start]  # i.e. if the C box is within the first 15-19 nt of the window
                left_internal_nt = whole_window[sno_start:sno_start+5+(15-sno_start)]
                line_left = left_external_nt + left_internal_nt
            else:
                left_internal_5nt = line[0:5]
                left_external_15nt = whole_window[sno_start-15:sno_start]
                line_left = left_external_15nt + left_internal_5nt

            sno_end = re.search(line, whole_window).end()
            if sno_end >= (int(fixed_length)-15):                # if the sno end is within the last 14-10 nt of the window
                right_external_nt = whole_window[sno_end:]  # i.e. if the D box is within the last 19-15 nt of the window
                right_internal_nt = whole_window[sno_end-5-(15-len(right_external_nt)):sno_end]
                line_right = right_internal_nt + right_external_nt
            else:
                right_internal_5nt = whole_window[sno_end-5:sno_end]
                right_external_15nt = whole_window[sno_end:sno_end+15]
                line_right = right_internal_5nt + right_external_15nt

            # Reverse the order of both flanking sequences with [::-1]
            co_seq = line_left[::-1] + "&" + line_right[::-1]  # this order is how RNAcofold will accurately try to see the
                                                # best base pairing between the left and right extended flanking regions of snoRNAs
            f2.write(f'>{id_}\n{co_seq}\n')


# Use RNAcofold to get terminal stem stability of positives and negatives
sp.call(f'RNAcofold < terminal_temp_pos.fa > {output_positives_fa}', shell=True)
sp.call(f'RNAcofold < terminal_temp_neg.fa > {output_negatives_fa}', shell=True)

# Extract the MFE of the predicted structure
sp.call(f"""grep -E ">" {output_positives_fa} | sed 's/>//g' > cofold_temp_id_pos && """+
        f"""grep -oE "\-*[0-9]+\.[0-9]*\)" {output_positives_fa} | sed 's/)//g' > cofold_temp_mfe_pos && """+
        f"""paste cofold_temp_id_pos cofold_temp_mfe_pos > {output_positives_tsv}""", shell=True)

sp.call(f"""grep -E ">" {output_negatives_fa} | sed 's/>//g' > cofold_temp_id_neg && """+
        f"""grep -oE "\-*[0-9]+\.[0-9]*\)" {output_negatives_fa} | sed 's/)//g' > cofold_temp_mfe_neg && """+
        f"""paste cofold_temp_id_neg cofold_temp_mfe_neg > {output_negatives_tsv}""", shell=True)

sp.call('rm *.ps terminal_temp_* cofold_temp*', shell=True)