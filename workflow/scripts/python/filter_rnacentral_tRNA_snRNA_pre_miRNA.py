#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp
import re 

output = snakemake.output.df
rnacentral_df_paths = snakemake.input.beds
genomes = snakemake.input.genomes
chr_sizes = snakemake.input.chr_sizes
extension = snakemake.params.extension
gff_gallus = pd.read_csv(snakemake.input.gallus_gallus_gff, sep='\t', 
                names=['chr', 'source', 'feature', 'start', 'end', 
                        'dot', 'strand', 'dot2', 'attributes'], skiprows=1)

# Filter and drop duplicates in gff_gallus (this will be used to convert old to new genomic coordinates)
# RNAcentral uses an old version of the gallus gallus genome, which is why we need to convert                       
gff_gallus = gff_gallus[(gff_gallus['feature'] == 'transcript') & 
                (gff_gallus['attributes'].str.contains('tRNA|snRNA|pre_miRNA'))]
gff_gallus = gff_gallus.drop_duplicates(subset=['chr', 'start', 'end', 'strand'])
gff_gallus[['name', 'type', 'db', 'id', 'source', 'identity']] = gff_gallus['attributes'].str.split(';', expand=True)
gff_gallus['id'] = gff_gallus['id'].str.replace('ID=', '')
gff_gallus['id'] = gff_gallus['id'].str.replace('(.*)\.[0-9]*', r'\1', regex=True)
gff_gallus = gff_gallus.drop_duplicates(subset=['id'])

chr_dictio, start_dictio, end_dictio, strand_dictio = {}, {}, {}, {}
for att in gff_gallus.attributes:
    gene_id = att.split(';')[0].split('=')[-1]
    chr = gff_gallus[gff_gallus['attributes'].str.contains(gene_id)]['chr'].values[0]
    start = gff_gallus[gff_gallus['attributes'].str.contains(gene_id)]['start'].values[0]
    end = gff_gallus[gff_gallus['attributes'].str.contains(gene_id)]['end'].values[0]
    strand = gff_gallus[gff_gallus['attributes'].str.contains(gene_id)]['strand'].values[0]
    chr_dictio[gene_id], start_dictio[gene_id] = chr, start
    end_dictio[gene_id], strand_dictio[gene_id] = end, strand
     

# Load beds as df
dfs = []
for path in list(pd.unique(rnacentral_df_paths)):
    species_name = path.split('/')[-1].split('.')[0]
    genome = [path for path in genomes if species_name in path][0]
    chr_size = [path for path in chr_sizes if species_name in path][0]
    bed_df = pd.read_csv(path, sep='\t', names=['chr', 'start', 'end', 
                        'rnacentral_id', 'frame', 'strand', 'start2', 
                        'end2', 'other', 'start3', 'end3', 'score', 
                        'dot', 'gene_biotype', 'gene_source'])
    
    # Drop duplicate entries
    bed_df = bed_df.drop_duplicates(subset=['chr', 'start', 'end', 'strand'])
    bed_df = bed_df[['chr', 'start', 'end', 'rnacentral_id', 'dot', 
                    'strand', 'score', 'gene_source', 'gene_biotype', 'frame']]
    # Convert location of ncRNA for G. gallus to latest genome version
    if species_name == 'gallus_gallus':
        bed_df['chr'] = bed_df['rnacentral_id'].map(chr_dictio)
        bed_df['start'] = bed_df['rnacentral_id'].map(start_dictio)
        bed_df['end'] = bed_df['rnacentral_id'].map(end_dictio)
        bed_df['strand'] = bed_df['rnacentral_id'].map(strand_dictio)
        bed_df = bed_df[bed_df['gene_biotype'].isin(['snRNA', 'tRNA', 'pre_miRNA'])].dropna(subset=['start', 'end'])
        bed_df[['start', 'end']] = bed_df[['start', 'end']].astype(int)
        
    # Modify chr column (w/r to Mt DNA and remove 'chr' in front of chr number)
    if species_name in ['macaca_mulatta', 'ornithorhynchus_anatinus', 'homo_sapiens', 'mus_musculus']:
        bed_df['chr'] = bed_df['chr'].str.replace('chrM', 'chrMT')
    if species_name != 'gallus_gallus':
        bed_df['chr'] = bed_df['chr'].str.replace('chr', '')
    # Merge overlapping entries
    bed_df['attributes'] = bed_df['gene_biotype']
    bed_df['dot2'] = '.'
    bed_df = bed_df[['chr', 'start', 'end', 'rnacentral_id', 'dot2', 
                    'strand', 'score', 'gene_source', 'gene_biotype']]
    bed_df.to_csv('temp_rnacentral.bed', sep='\t', index=False, header=False)
    sp.call('sort -k1,1 -k2,2n temp_rnacentral.bed | grep -vE "MTL14|QNVO|CHR_HG|CHR_HS|HSCHR|MTU069"  > temp_rnacentral.sorted.bed', shell=True)  # remove weird scaffolds
    temp_bed = BedTool('temp_rnacentral.sorted.bed')
    merged_bed_real = temp_bed.merge(s=True, c=[4, 5, 6, 8, 9], o='distinct').saveas('temp_rnacentral.sorted.merged.bed')
    # Get sequence of the ncRNAs into dict
    fasta = merged_bed_real.sequence(fi=genome, nameOnly=True, s=True)
    dictio_seq = {}
    with open(fasta.seqfn, 'r') as fasta_file:
        for line in fasta_file:
            if '>' in line:
                ncRNA_id = line.replace('>', '').replace('\n', '').replace(',', '_')
                ncRNA_id = re.sub("\(.*\)", "", ncRNA_id)
            else:
                seq = line.replace('\n', '')
                dictio_seq[ncRNA_id] = seq
    # Remove overlapping entries (ex: URS0000630DA2_9258,URS000231C311_9258)
    merged_bed = pd.read_csv('temp_rnacentral.sorted.merged.bed', sep='\t', 
                            names=['chr', 'start', 'end', 'rnacentral_id', 'dot',
                                    'strand', 'gene_source', 'gene_biotype'])
    merged_bed = merged_bed[~merged_bed.rnacentral_id.str.contains(',')]
    
    # Add sequence column and species
    merged_bed['species'] = species_name
    merged_bed['sequence'] = merged_bed['rnacentral_id'].map(dictio_seq)

    # Add extended sequence column (15 flanking nt up/downstream of ncRNA)
    extended_ncRNA = merged_bed_real.slop(g=chr_size, r=extension, l=extension)
    extended_sequences_fa = extended_ncRNA.sequence(fi=genome, nameOnly=True, s=True)
    ext_seq_dict = {}
    with open(extended_sequences_fa.seqfn, 'r') as fa:
        for line in fa:
            if '>' in line:
                ncRNA_name = line.strip('\n').replace('>', '').replace('(+)', '').replace('(-)', '')
                ncRNA_name = re.sub("\(.*\)", "", ncRNA_name)
            else:
                seq = line.strip('\n')
                ext_seq_dict[ncRNA_name] = seq

    merged_bed['extended_sequence'] = merged_bed['rnacentral_id'].map(ext_seq_dict)
    
    # Keep only tRNAs, snRNAs and pre_miRNA
    merged_bed = merged_bed[merged_bed.gene_biotype.isin(['tRNA', 'snRNA', 'pre_miRNA'])]
    
    dfs.append(merged_bed)

    
# Concat all dfs and drop duplicate based on sequence
final_df = pd.concat(dfs)
final_df = final_df.drop_duplicates(subset=['sequence'])
final_df = final_df.drop_duplicates(subset=['extended_sequence'])


final_df.to_csv(output, sep='\t', index=False)

sp.call('rm temp_rnacentral.*', shell=True)

