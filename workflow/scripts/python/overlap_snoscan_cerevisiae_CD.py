#!/usr/bin/python3
import pandas as pd
import re
import subprocess as sp 
from pybedtools import BedTool

preds = snakemake.input.snoscan
output = snakemake.output.df

cd_sno = pd.read_csv(snakemake.input.cd_sno, sep='\t')
cd_yeast = cd_sno[cd_sno['species_name'] == 'S_cerevisiae']
cd_yeast['score'] = '.'
cd_yeast = cd_yeast[['chr', 'start', 'end', 'gene_id', 'score', 'strand']]
cd_yeast.to_csv('temp_cd_yeasttt.bed', sep='\t', index=False, header=False)
cd_yeast_bed = BedTool('temp_cd_yeasttt.bed')

positives_ = []
# Extract which examples were predicted as snoRNAs by snoscan
with open(preds, 'r') as f:
    for i, line in enumerate(f):
        if line.startswith('>>'):
            # Search for numbers in parentheses separated by a '-'
            start, end = re.search(r'\(\d+-\d+\)', line).group().strip(')()').split('-')
            if int(start) > int(end):
                strand = '-'
                start2 = int(end)
                end2 = int(start)
            else:
                strand = '+'
                start2 = int(start)
                end2 = int(end)
            chr_ = line.split('  ')[0].strip('> ')
            gene_id = f'snoscan_{i+1}'
            positives_.append([chr_, start2, end2, gene_id, '.', strand])


# Create bed of CD predictions
bed_df = pd.DataFrame(positives_, columns=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])


# Drop duplicate predictions (though with different scores)
bed_df = bed_df.drop_duplicates(subset=['chr', 'start', 'end'])
bed_df = bed_df.sort_values(by=['chr', 'start'])


# Merge overlapping predictions
bed_df.to_csv('temp_snoscan_cd.bed', sep='\t', index=False, header=False)
bed_pred = BedTool('temp_snoscan_cd.bed')
bed_pred = bed_pred.merge(s=True, c=[4, 5, 6], o=['distinct', 'distinct', 'distinct'])

# Intersect predictions with annotated CD (at least 50% of the C/D must be overlapping)
inter = bed_pred.intersect(cd_yeast_bed, s=True, F=0.5).to_dataframe()
TP = len(inter)
FN = len(cd_yeast) - TP
FP = len(bed_df) - TP
precision = TP/(TP + FP)
recall = TP/(TP+FN)
fscore = 2*precision*recall/(precision + recall)
final_df = pd.DataFrame([['snoscan', len(cd_yeast), TP, FN, FP, precision, recall, fscore]], 
            columns=['predictor', 'Real_nb_CD', 'predicted_nb_CD', 'FN', 'FP', 'precision', 'recall', 'fscore'])

final_df.to_csv(snakemake.output.df, sep='\t', index=False)

sp.call('rm temp_snoscan_cd.bed temp_cd_yeasttt.bed', shell=True)
