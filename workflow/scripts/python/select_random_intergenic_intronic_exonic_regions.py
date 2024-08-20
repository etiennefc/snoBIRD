#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp
import numpy as np

species = str(snakemake.wildcards.species)
rs = snakemake.params.random_state
np.random.seed(rs) # set a fixed reproducible randomness
positives = pd.read_csv(snakemake.input.expressed_cd_all_sets, sep='\t')
intergenic_regions = pd.read_csv(snakemake.input.intergenic_regions, 
                        sep='\t', names=['chr', 'start', 'end'])
intronic_regions = pd.read_csv(snakemake.input.intronic_regions, 
                        sep='\t', names=['chr', 'start', 'end'])
exonic_regions = pd.read_csv(snakemake.input.exonic_regions, 
                        sep='\t', names=['chr', 'start', 'end', 'gene_id', 'score', 
                                        'strand', 'feature', 'gene_biotype'])
exonic_regions = exonic_regions[['chr', 'start', 'end', 'strand']]
bed_cols = ['chr', 'start', 'end', 'gene_id', 'score', 'strand', 'gene_biotype']

output_intergenic = snakemake.output.random_intergenic_regions
output_intronic = snakemake.output.random_intronic_regions
output_exonic = snakemake.output.random_exonic_regions


# Get length of all positive examples (expressed C/D snoRNAs with flanking 15 nt) #Extended sequence!
seqs = list(positives.extended_sequence)
lengths = [len(seq) for seq in seqs]
positives['len'] = positives.end - positives.start

# The intergenic/intronic regions are given without strand, so
# we need to define a random choice between + or - for these
strand_dict = {0: '+', 1: '-'}
intergenic_strand = list(np.random.choice([0, 1], 
                        size=len(intergenic_regions), p=[0.5, 0.5]))
intronic_strand = list(np.random.choice([0, 1], 
                        size=len(intronic_regions), p=[0.5, 0.5]))
intronic_regions['strand'] = intronic_strand
intronic_regions['strand'] = intronic_regions['strand'].map(strand_dict)
intergenic_regions['strand'] = intergenic_strand
intergenic_regions['strand'] = intergenic_regions['strand'].map(strand_dict)


def select_regions(df, location, output):
    """ Select n intergenic/intronic regions of size s in the pool of regions, 
        where n is the number of positive examples and s is a size in the 
        distribution of positive example lengths."""
    # Select only regions that are at least as large as the longest positive
    df['len'] = df['end'] - df['start']
    df = df[df['len'] >= max(lengths) + 2] # regions are least greater than 1 nt than sno length

    # Select randomly n regions of size s (where n=nb of positives 
    # and s is part of the size distribution of n)
    random_df = df.sample(n=len(positives), random_state=rs)
    random_df = random_df.reset_index(drop=True)
    rows = []
    for i, length in enumerate(lengths):
        row_dict = dict(random_df.iloc[i, :])
        # Get a random start/end in the given region that respects the overall sno length
        random_start_index = int(np.random.choice(range(0, row_dict['len'] - length)))
        random_start = row_dict['start'] + random_start_index
        random_end = row_dict['start'] + random_start_index + length
        region_id = f'{species}_{location}_region_{i}'
        row = [row_dict['chr'], random_start, random_end, region_id, '.', 
                row_dict['strand'], f'random_{location}_region']
        rows.append(row)

    selected_region_df = pd.DataFrame(rows, columns=bed_cols)
    selected_region_df.to_csv(f'{location}_regions_{species}_temp.bed', 
                                        sep='\t', index=False, header=False)

    # Get sequence of selected regions and add it to df
    selected_region_bed = BedTool(f'{location}_regions_{species}_temp.bed')
    fasta = selected_region_bed.sequence(fi=snakemake.input.genome, nameOnly=True, s=True)
    seq_dict = {}
    with open(fasta.seqfn, 'r') as fasta_file:
        for line in fasta_file:
            if '>' in line:
                name = line.strip('\n').replace('>', '').replace('(+)', '').replace('(-)', '')
            else:
                seq = line.strip('\n')
                seq_dict[name] = seq

    selected_region_df['sequence'] = selected_region_df['gene_id'].map(seq_dict)
    
    selected_region_df.to_csv(output, sep='\t', index=False)
    sp.call(f'rm {location}_regions_{species}_temp.bed', shell=True)

# Select intergenic regions
select_regions(intergenic_regions, 'intergenic', output_intergenic)

# Select exonic regions
select_regions(exonic_regions, 'exonic', output_exonic)

# Select intronic regions for species with a high number of large introns
if species not in ['giardia_lamblia', 'leishmania_major', 
                    'saccharomyces_cerevisiae', 
                    'dictyostelium_discoideum']:
    select_regions(intronic_regions, 'intronic', output_intronic)
# There is a limited number of intronic regions in G. lamblia, L. major, S.cerevisiae and D. discoideum
# (not enough to have equal number compared to expressed CD)
else:
    df = pd.DataFrame(columns=['Not enough large introns to select introns'])
    df.to_csv(output_intronic, sep='\t', index=False)







