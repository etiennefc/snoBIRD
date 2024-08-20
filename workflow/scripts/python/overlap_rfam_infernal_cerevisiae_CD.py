#!/usr/bin/python3
import pandas as pd
import re
import subprocess as sp 
from pybedtools import BedTool

preds = snakemake.input.infernal_rfam
output = snakemake.output.df
haca_yeast = pd.read_csv(snakemake.input.haca_yeast, sep='\t')
haca_yeast = list(haca_yeast[haca_yeast['sno_type'] == 'H/ACA'].gene_id)

cd_sno = pd.read_csv(snakemake.input.cd_sno, sep='\t')
cd_yeast = cd_sno[cd_sno['species_name'] == 'S_cerevisiae']
cd_yeast['score'] = '.'
cd_yeast = cd_yeast[['chr', 'start', 'end', 'gene_id', 'score', 'strand']]
cd_yeast.to_csv('temp_cd_yeastt.bed', sep='\t', index=False, header=False)
cd_yeast_bed = BedTool('temp_cd_yeastt.bed')

positives_, positives_dict = [], {}
with open(preds, 'r') as f:
    i = 0
    for line in f:
        # Select only lines with C/D snoRNAs for which a rfam family was identified
        if ('#' and 'SNORA' not in line) & ('Small nucleolar RNA' in line):
            gene_id = re.split('Small nucleolar RNA ', line)[1].strip('\n')
            
            chr_ = re.split('RF[0-9]{5}   ', line)[1].split(' ')[0]
            if '+' in line:
                strand = '+'
                start2 = line.split(strand)[0].replace(' ', '%%')
                start = re.split('%%*', start2)[-3]
                end = re.split('%%*', start2)[-2]
            else:
                strand = '-'
                start2 = line.split(strand)[1].replace(' ', '%%')
                end = re.split('%%*', start2)[-3]
                start = re.split('%%*', start2)[-2]
            if (gene_id not in positives_dict.keys()) & (gene_id not in haca_yeast):
                positives_dict['gene_id'] = 'predicted_CD'
                positives_.append([chr_, start, end, gene_id, '.', strand])

# Create bed of CD predictions
bed_df = pd.DataFrame(positives_, columns=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])
bed_df.to_csv('temp_infernal_cd.bed', sep='\t', index=False, header=False)
bed_pred = BedTool('temp_infernal_cd.bed')

# Intersect predictions with annotated CD
inter = bed_pred.intersect(cd_yeast_bed, s=True, F=0.5).to_dataframe()
TP = len(inter)
FN = len(cd_yeast) - TP
FP = len(bed_df) - TP
precision = TP/(TP + FP)
recall = TP/(TP+FN)
fscore = 2*precision*recall/(precision + recall)
final_df = pd.DataFrame([['infernal_rfam', len(cd_yeast), TP, FN, FP, precision, recall, fscore]], 
            columns=['predictor', 'Real_nb_CD', 'predicted_nb_CD', 'FN', 'FP', 'precision', 'recall', 'fscore'])

final_df.to_csv(snakemake.output.df, sep='\t', index=False)

sp.call('rm temp_infernal_cd.bed temp_cd_yeastt.bed', shell=True)
