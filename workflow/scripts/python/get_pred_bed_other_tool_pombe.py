#!/usr/bin/python3
import pandas as pd
import subprocess as sp
import numpy as np 
from Bio import SeqIO 


bed_cols = ['seqname', 'start', 'end', 'block_id', 'score', 
            'strand', 'block_length']

infernal_pred = snakemake.input.infernal 
snoreport_pred_pos = snakemake.input.snoreport_pos 
snoreport_pred_neg = snakemake.input.snoreport_neg 
snoscan_pred = snakemake.input.snoscan
genome = snakemake.input.genome_fa
sno_rfam = pd.read_csv(snakemake.input.sno_rfam, sep='\t')
cd_rfam  = sno_rfam[sno_rfam['snoRNA_type'] == 'C/D']

# Filter infernal predictions
sp.call(f"awk -v OFS='\t' '/nucleolar/ {{print $3,$8,$9,$2,$16,$10,$9-$8}}' {infernal_pred} > infernal_pombe_temp.tsv", shell=True)
infernal_sno = pd.read_csv('infernal_pombe_temp.tsv', sep='\t', names=bed_cols)
cd_infernal = infernal_sno[infernal_sno['block_id'].isin(cd_rfam.Rfam_family_id)].sort_values(by=['seqname', 'start']).reset_index(drop=True)
cd_infernal['block_length'] = np.abs(cd_infernal['block_length'])

# Switch start and end for predictions on the - strand
rows = []
for row in cd_infernal.iterrows():
    if row[1].strand == '-':
        rows.append([row[1].seqname, row[1].end, row[1].start, row[1].block_id, 
                row[1].score, row[1].strand, row[1].block_length])
    else:
        rows.append(list(row[1]))
cd_infernal = pd.DataFrame(rows, columns=bed_cols)
cd_infernal.to_csv(snakemake.output.bed_infernal, sep='\t', index=False, header=None)


# Filter snoscan predictions
i = 0
snoscan_rows = []
with open(snoscan_pred, 'r') as f:
    for line in f:
        if line.startswith('>>'):
            i += 1
            chr_ = line.split(' ')[1]
            score = float(line.split(' ')[3])
            st, nd = line.split(' ')[5].strip('()').split('-')
            if int(st) < int(nd):  # + strand prediction
                start, end, strand = int(st), int(nd), '+'
            else:
                start, end, strand = int(nd), int(st), '-'
            length = end - start
            pred_id = f'snoscan_{i}_{strand}'
            snoscan_rows.append([chr_, start, end, pred_id, score, strand, length])
            
snoscan = pd.DataFrame(snoscan_rows, columns=bed_cols).sort_values(by=['seqname', 'start', 'score'], 
                        ascending = [True, True, False]).reset_index(drop=True)
snoscan.to_csv(snakemake.output.bed_snoscan, sep='\t', index=False, header=None)

# For the predictions which have multiple targets, keep only the most probable one
snoscan_filtered = snoscan.drop_duplicates(subset=['seqname', 'start', 'end', 'strand'])
snoscan_filtered.to_csv(snakemake.output.bed_snoscan_filtered_target, sep='\t', index=False, header=None)

# Get chr size for the prediction on the - strand for snoreport2
chr_dict = {record.id: len(str(record.seq)) 
                for record in SeqIO.parse(genome, "fasta")}
# Filter snoreport predictions on + then - strand
j, snoreport_rows = 0, []
with open(snoreport_pred_pos, 'r') as f:
    strand = '+'
    for line in f:
        if '>' in line:
            j += 1
            chr_ = line.split('_')[0].strip('>')
            start = int(line.split(' ')[3])
            end = int(line.split(' ')[4])
            score = float(line.split(' ')[9])
            length = end - start 
            pred_id = f'snoreport_{j}_{strand}'
            snoreport_rows.append([chr_, start, end, pred_id, score, strand, length])

with open(snoreport_pred_neg, 'r') as f:
    strand = '-'
    for line in f:
        if '>' in line:
            j += 1
            chr_ = line.split('_')[0].strip('>')
            chr_len = chr_dict[chr_]
            start = chr_len - int(line.split(' ')[4])
            end = chr_len - int(line.split(' ')[3])
            score = float(line.split(' ')[9])
            length = end - start 
            pred_id = f'snoreport_{j}_{strand}'
            snoreport_rows.append([chr_, start, end, pred_id, score, strand, length])



snoreport = pd.DataFrame(snoreport_rows, columns=bed_cols).sort_values(by=['seqname', 'start', 'score'], 
                        ascending = [True, True, False]).reset_index(drop=True)
snoreport.to_csv(snakemake.output.bed_snoreport, sep='\t', index=False, header=None)

sp.call('rm infernal_pombe_temp.tsv', shell=True)

