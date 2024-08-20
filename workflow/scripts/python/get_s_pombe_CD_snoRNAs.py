#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp
import re
from gtfparse import read_gtf

bed_cols = ['seqname', 'start', 'end', 'gene_id', 'score', 
            'strand', 'feature', 'snoRNA_type']

# Create gtf of only the genes in S. pombe
cmd1 = f"""awk 'NR>5 && $3=="gene"' {snakemake.input.gtf} > temp_s_pombe.gtf"""
sp.call(cmd1, shell=True)

# Filter that gtf to keep only snoRNAs found in Pombase (this removes duplicates introduced in the Ensembl annotation)
cmd2 = f"""for i in $(grep "snoRNA gene" {snakemake.input.s_pombe_genes} | awk '{{print $1}}'); """+\
        """do grep $i temp_s_pombe.gtf >> sno_s_pombe.gtf; done"""
sp.call(cmd2, shell=True)
sno_gtf = read_gtf('sno_s_pombe.gtf')


# Find if these snoRNAs are C/D or H/ACA
sno_type_df = pd.read_csv(snakemake.input.sno_type_pombe, sep='\t')
sno_type_df['gene_id'] = sno_type_df['gene_id_sno']


# Create bed out of the sno gtf and keep only C/D
cd_gtf = sno_gtf.merge(sno_type_df[['gene_id', 'snoRNA_type']], how='left', on='gene_id')
cd_gtf = cd_gtf[cd_gtf['snoRNA_type'] == 'C/D']
cd_gtf[bed_cols].to_csv('tempo_bed_s_pombe.bed', header=False, sep='\t', index=False)
sp.call(f"""sort -k1,1 -k2,2n tempo_bed_s_pombe.bed | awk -v OFS="\t" '{{print $1,$2,$3,$4,".",$5,$6,$7}}' > {snakemake.output.bed}""", shell=True)

sno_bed = BedTool(snakemake.output.bed)

# Create fasta of the sequences of these sno
fasta = sno_bed.sequence(fi=snakemake.input.genome, nameOnly=True, s=True)
seq_dict = {}
with open(fasta.seqfn, 'r') as fasta_file:
    for line in fasta_file:
        if '>' in line:
            sno_name = line.strip('\n').replace('>', '').replace('(+)', '').replace('(-)', '')
        else:
            seq = line.strip('\n')
            seq_dict[sno_name] = seq
with open(snakemake.output.fa, 'w') as f:
    for k,v in seq_dict.items():
        f.write(f'>{k}\n{v}\n')



sp.call('rm tempo_bed_s_pombe.bed temp_s_pombe.gtf sno_s_pombe.gtf', shell=True)
