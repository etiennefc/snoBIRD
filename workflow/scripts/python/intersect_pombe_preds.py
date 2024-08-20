#!/usr/bin/python3
import pandas as pd
import subprocess as sp
from pybedtools import BedTool


bed_cols = ['chr', 'start', 'end', 'gene_id', 'score', 
            'strand', 'feature', 'sno_type']
cols = bed_cols + ['chr_pred', 'start_pred', 'end_pred', 'gene_id_pred', 'score_pred', 
            'strand_pred', 'len']

infernal_pred = BedTool(snakemake.input.infernal) 
snoreport_pred = BedTool(snakemake.input.snoreport) 
snoscan_pred = BedTool(snakemake.input.snoscan)
snoBIRD_pred = BedTool(snakemake.input.snoBIRD)
cd_pombe = pd.read_csv(snakemake.input.cd_pombe, sep='\t', names=bed_cols)
fixed_length = snakemake.params.fixed_length


# Filter to keep only CD of length of fixed_length nt
cd_pombe = cd_pombe[cd_pombe['end'] - cd_pombe['start'] <= fixed_length - 15 * 2].reset_index(drop=True)
cd_pombe[['start', 'end']] = cd_pombe[['start', 'end']].astype(int)
cd_pombe.to_csv('s_pombe_cd_temp.bed', sep='\t', index=False, header=False)
cd_pombe_bed = BedTool('s_pombe_cd_temp.bed')


# Intersect cd_pombe with the predictions from the other tools
infernal_intersect = cd_pombe_bed.intersect(infernal_pred, s=True, f=0.95, wb=True, wa=True).to_dataframe(names=cols)
infernal_intersect = infernal_intersect.drop_duplicates(subset=['start', 'end', 'gene_id']).reset_index(drop=True)
infernal_intersect.to_csv(snakemake.output.bed_infernal, sep='\t', index=False)

snoreport_intersect = cd_pombe_bed.intersect(snoreport_pred, s=True, f=0.95, wb=True, wa=True).to_dataframe(names=cols)
snoreport_intersect = snoreport_intersect.drop_duplicates(subset=['start', 'end', 'gene_id']).reset_index(drop=True)
snoreport_intersect.to_csv(snakemake.output.bed_snoreport, sep='\t', index=False)

snoscan_intersect = cd_pombe_bed.intersect(snoscan_pred, s=True, f=0.95, wb=True, wa=True).to_dataframe(names=cols)
snoscan_intersect = snoscan_intersect.drop_duplicates(subset=['start', 'end', 'gene_id']).reset_index(drop=True)
snoscan_intersect.to_csv(snakemake.output.bed_snoscan, sep='\t', index=False)

snoBIRD_intersect = cd_pombe_bed.intersect(snoBIRD_pred, s=True, f=0.95, wb=True, wa=True).to_dataframe(names=cols)
snoBIRD_intersect = snoBIRD_intersect.drop_duplicates(subset=['start', 'end', 'gene_id']).reset_index(drop=True)
snoBIRD_intersect.to_csv(snakemake.output.bed_snoBIRD, sep='\t', index=False)

sp.call('rm s_pombe_cd_temp.bed', shell=True)