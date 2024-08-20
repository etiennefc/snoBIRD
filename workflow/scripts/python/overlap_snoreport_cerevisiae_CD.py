#!/usr/bin/python3
import pandas as pd
import re
import subprocess as sp 
from pybedtools import BedTool

preds_pos = snakemake.input.snoreport_pos
preds_neg = snakemake.input.snoreport_neg
chr_size = pd.read_csv(snakemake.input.chr_size, sep='\t', names=['chr', 'size'])
size_dict = dict(zip(chr_size['chr'], chr_size['size']))

cd_sno = pd.read_csv(snakemake.input.cd_sno, sep='\t')
cd_yeast = cd_sno[cd_sno['species_name'] == 'S_cerevisiae']
cd_yeast['score'] = '.'
cd_yeast = cd_yeast[['chr', 'start', 'end', 'gene_id', 'score', 'strand']]
cd_yeast.to_csv('temp_cd_yeastttt.bed', sep='\t', index=False, header=False)
cd_yeast_bed = BedTool('temp_cd_yeastttt.bed')

positives_ = []
# Extract which examples were predicted as snoRNAs by snoreport on + strand
with open(preds_pos, 'r') as f:
    strand = '+'
    for line in f:
        if '>' in line:
            gene_id = line.split(' ')[0].replace('>', '')
            print(gene_id)
            chr_ = gene_id.split('_')[0]
            print(line)
            start = int(line.split(' ')[3])
            end = int(line.split(' ')[4])
            print(chr_, start, end, gene_id, strand)
            positives_.append([chr_, start, end, gene_id, '.', strand])

# Extract which examples were predicted as snoRNAs by snoreport on - strand
with open(preds_neg, 'r') as f:
    strand = '-'
    for line in f:
        if '>' in line:
            gene_id = line.split(' ')[0].replace('>', '')
            chr_ = gene_id.split('_')[0]
            c_size = size_dict[chr_]
            # Convert positions on negative strand with regards to the + strand (left to right)
            fake_start = int(line.split(' ')[3])
            fake_end = int(line.split(' ')[4])
            start = c_size - fake_end
            end = c_size - fake_start
            positives_.append([chr_, start, end, gene_id, '.', strand])



# Create bed of CD predictions
bed_df = pd.DataFrame(positives_, columns=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])


# Drop duplicate predictions (though with different scores)
bed_df = bed_df.drop_duplicates(subset=['chr', 'start', 'end'])
bed_df = bed_df.sort_values(by=['chr', 'start'])


# Merge overlapping predictions
bed_df.to_csv('temp_snoreport_cd.bed', sep='\t', index=False, header=False)
bed_pred = BedTool('temp_snoreport_cd.bed')
bed_pred = bed_pred.merge(s=True, c=[4, 5, 6], o=['distinct', 'distinct', 'distinct'])

# Intersect predictions with annotated CD (at least 50% of the C/D must be overlapping)
inter = bed_pred.intersect(cd_yeast_bed, s=True, F=0.5).to_dataframe()
print(inter)
TP = len(inter)
FN = len(cd_yeast) - TP
FP = len(bed_df) - TP
precision = TP/(TP + FP)
recall = TP/(TP+FN)
fscore = 2*precision*recall/(precision + recall)
final_df = pd.DataFrame([['snoreport', len(cd_yeast), TP, FN, FP, precision, recall, fscore]], 
            columns=['predictor', 'Real_nb_CD', 'predicted_nb_CD', 'FN', 'FP', 'precision', 'recall', 'fscore'])

final_df.to_csv(snakemake.output.df, sep='\t', index=False)

print(final_df)

sp.call('rm temp_snoreport_cd.bed temp_cd_yeastttt.bed', shell=True)
