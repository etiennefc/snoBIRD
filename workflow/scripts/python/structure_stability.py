#!/usr/bin/python3
import pandas as pd
import subprocess as sp 

fixed_length = snakemake.wildcards.fixed_length
output_positives_fa = snakemake.output.positives_fa
output_negatives_fa = snakemake.output.negatives_fa
output_positives_tsv = snakemake.output.positives_tsv
output_negatives_tsv = snakemake.output.negatives_tsv

# Load dfs and create box location dicts
positives_df = pd.read_csv(snakemake.input.positives_df, sep='\t')
negatives_df = pd.concat([pd.read_csv(path, sep='\t') for path in snakemake.input.negatives_df])
box_pos = pd.read_csv(snakemake.input.box_score_positives, sep='\t')
C_box_pos_dict = dict(zip(box_pos.gene_id, box_pos.C_start))
D_box_pos_dict = dict(zip(box_pos.gene_id, box_pos.D_end))

box_rest = pd.read_csv(snakemake.input.box_score_negatives, sep='\t')
C_box_rest_dict = dict(zip(box_rest.gene_id, box_rest.C_start))
D_box_rest_dict = dict(zip(box_rest.gene_id, box_rest.D_end))

# Function to set a limit of sequence extension 
# (defines the start and end of snoRNA in a window)
def lim_ext(sequence, c_start_, d_end_, extension):
    l = len(sequence)
    if c_start == 0:  # no C box was found
        sno_start = 16  # minimal start of sno in windows due to 15 extended nt on each side
    else:
        sno_start = c_start_ - extension
    if d_end_ == 0:  # no D box was found
        sno_end = l - 15
    else:
        sno_end = d_end_ + extension
    return sno_start, sno_end



# Extract sno sequence from window based on the boxes
with open('temp_structure_pos.fa', 'w') as temp:
    for id_, seq in dict(zip(positives_df.gene_id, positives_df[f'extended_{fixed_length}nt_sequence'])).items():
        c_start, d_end = C_box_pos_dict[id_], D_box_pos_dict[id_]

        # If possible extend 5 nt before c_start (and 5 after d_end) to get the complete snoRNA
        sno_start, sno_end = lim_ext(seq, c_start, d_end, 5)
        sno_seq = seq[sno_start-1:sno_end]
        temp.write(f'>{id_}\n{sno_seq}\n')

with open('temp_structure_neg.fa', 'w') as temp2:
    for id_, seq in dict(zip(negatives_df.gene_id, negatives_df[f'extended_{fixed_length}nt_sequence'])).items():
        c_start, d_end = C_box_rest_dict[id_], D_box_rest_dict[id_]

        # If possible extend 5 nt before c_start (and 5 after d_end) to get the complete snoRNA
        sno_start, sno_end = lim_ext(seq, c_start, d_end, 5)
        sno_seq = seq[sno_start-1:sno_end]
        temp2.write(f'>{id_}\n{sno_seq}\n')
            
            
# Predict the secondary structure stability of the snoRNA (it converts by default to RNA)
sp.call(f'RNAfold --infile=temp_structure_pos.fa --outfile={output_positives_fa} && mv data_references_structure_* {output_positives_fa}', shell=True)
sp.call(f'RNAfold --infile=temp_structure_neg.fa --outfile={output_negatives_fa} && mv data_references_structure_* {output_negatives_fa}', shell=True)

# Extract the MFE of the predicted structure
sp.call(f"""grep -E ">" {output_positives_fa} | sed 's/>//g' > temp_id_pos && """+
        f"""grep -oE "\-*[0-9]+\.[0-9]*\)" {output_positives_fa} | sed 's/)//g' > temp_mfe_pos && """+
        f"""paste temp_id_pos temp_mfe_pos > {output_positives_tsv}""", shell=True)

sp.call(f"""grep -E ">" {output_negatives_fa} | sed 's/>//g' > temp_id_neg && """+
        f"""grep -oE "\-*[0-9]+\.[0-9]*\)" {output_negatives_fa} | sed 's/)//g' > temp_mfe_neg && """+
        f"""paste temp_id_neg temp_mfe_neg > {output_negatives_tsv}""", shell=True)

sp.call('rm temp_structure* *.ps temp_id_* temp_mfe_*', shell=True)