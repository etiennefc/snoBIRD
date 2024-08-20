#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp
from gtfparse import read_gtf

species = str(snakemake.wildcards.species)
gtf = read_gtf(snakemake.input.gtf)
species_dict = snakemake.params.species_dict
short_species_name = [k for k,v in species_dict.items() if v==species][0]
chr_size_file = [path for path in snakemake.input.genome_size 
                if species in path][0]
cd_path = snakemake.input.expressed_cd_all_sets
ncRNA_path = snakemake.input.ncRNA
human_pseudo_path = snakemake.input.human_pseudogene
mouse_pseudo_path = snakemake.input.mouse_pseudogene
droso_pseudo_path = snakemake.input.droso_pseudogene

if species == 'tetrahymena_thermophila':
    haca_path = None  # There are no annotated H/ACA snoRNAs in that species
else:
    haca_path = [path for path in snakemake.input.haca if species in path][0]
intronic_output = snakemake.output.intronic_regions_bed
intergenic_output = snakemake.output.intergenic_regions_bed
exonic_output = snakemake.output.exonic_regions_bed

bed_cols = ['seqname', 'start', 'end', 'gene_id', 'score', 
            'strand', 'feature', 'gene_biotype']
df_cols = ['chr', 'start', 'end', 'gene_id', 
            'score', 'strand', 'gene_biotype']

# Load cd df
cd_df = pd.read_csv(str(cd_path), sep='\t')
cd_df['score'] = '.'
cd_df = cd_df[cd_df['species_name'] == short_species_name]
cd_df['gene_biotype'] = 'C/D'
cd_df = cd_df[df_cols]

# Load ncRNA df
ncRNA_df = pd.read_csv(ncRNA_path, sep='\t')
ncRNA_df = ncRNA_df[ncRNA_df['species'] == species]
ncRNA_df[['score', 'gene_biotype']] = '.', 'ncRNA'
ncRNA_df = ncRNA_df.rename(columns={'rnacentral_id': 'gene_id'})
ncRNA_df = ncRNA_df[df_cols]

# Load H/ACA df
if species == 'tetrahymena_thermophila':
    haca_df = pd.DataFrame(columns=df_cols)  # create empty df because no annotated in H/ACA
else:
    haca_df = pd.read_csv(haca_path, sep='\t')
    haca_df[['score', 'gene_biotype']] = '.', 'H/ACA'
    haca_df = haca_df.rename(columns={'gene_name': 'gene_id'})
    haca_df = haca_df[df_cols]

# Load C/D pseudogene dfs
if species == 'homo_sapiens':
    human_pseudo = pd.read_csv(human_pseudo_path[0], sep='\t')
    human_pseudo['score'] = '.'
    human_pseudo = human_pseudo[df_cols]
    all_ncRNA_df = pd.concat([cd_df, ncRNA_df, haca_df, human_pseudo])
elif species == 'mus_musculus':
    mouse_pseudo = pd.read_csv(mouse_pseudo_path[0], sep='\t')
    mouse_pseudo['score'] = '.'
    mouse_pseudo = mouse_pseudo[df_cols]
    all_ncRNA_df = pd.concat([cd_df, ncRNA_df, haca_df, mouse_pseudo])
elif species == 'drosophila_melanogaster':
    droso_pseudo = pd.read_csv(droso_pseudo_path[0], sep='\t')
    droso_pseudo['score'] = '.'
    droso_pseudo = droso_pseudo[df_cols]
    all_ncRNA_df = pd.concat([cd_df, ncRNA_df, haca_df, droso_pseudo])
else:
    all_ncRNA_df = pd.concat([cd_df, ncRNA_df, haca_df])

def get_bed_ncRNA(a_ncRNA_df, species_name):
    """ Convert a ncRNA df into a temp bed"""

    a_ncRNA_df.to_csv(f'all_ncRNA_{species_name}.bed', sep='\t', index=False, header=False)
    sp.call(f'sort -k1,1 -k2,2n all_ncRNA_{species_name}.bed > all_ncRNA_{species_name}.sorted.bed', shell=True)
    all_ncRNA_bed = BedTool(f'all_ncRNA_{species_name}.sorted.bed')

    return all_ncRNA_bed


def get_intergenic_regions(gtf_df, species_name, chr_size, ncRNA_bed, output):
    """ Get intergenic regions in a genome by taking
        the bedtools complement of all genes (i.e. all the 
        regions not included in these genes)."""
    # Create bed of all genes
    gene_gtf = gtf_df[gtf_df['feature'] == 'gene']
    gene_gtf[bed_cols].to_csv(f'all_genes_{species_name}.bed', 
                        sep='\t', index=False, header=False)
    sp.call(f'sort -k1,1 -k2,2n all_genes_{species_name}.bed > all_genes_{species_name}_sorted.bed', 
            shell=True)
    gene_bed = BedTool(f'all_genes_{species_name}_sorted.bed')

    # Create temporary sorted chr size file
    sp.call(f'sort -k1,1 -k2,2n {chr_size} > temp_chr_size{species_name}.tsv', shell=True)

    # Get the complement of these genes
    complement = gene_bed.complement(g=f'temp_chr_size{species_name}.tsv')

    # Subtract to make sure all ncRNAs obtained before are not included in these regions
    final_intergenic_regions = complement.subtract(b=ncRNA_bed).saveas(output)

    return final_intergenic_regions


def get_intronic_regions(gtf_df, species_name, chr_size, intergenic_bed, ncRNA_bed, output):
    """ Get the intronic regions in biotypes of interest (intronic 
        regions that don't overlap with embedded genes)."""
    # Select intronic and intergenic regions using the exons of all genes as boundaries
    gene_gtf = gtf_df[gtf_df['feature'] == 'exon']
    gene_gtf[bed_cols].to_csv(f'all_exons_{species_name}.bed', 
                        sep='\t', index=False, header=False)
    sp.call(f'sort -k1,1 -k2,2n all_exons_{species_name}.bed > all_exons_{species_name}_sorted.bed', 
            shell=True)
    gene_bed = BedTool(f'all_exons_{species_name}_sorted.bed')

    # Create temporary sorted chr size file
    sp.call(f'sort -k1,1 -k2,2n {chr_size} > temp2_chr_size{species_name}.tsv', shell=True)

    # Get the complement of these exons (i.e intron and intergenic regions)
    complement = gene_bed.complement(g=f'temp2_chr_size{species_name}.tsv')

    # Subtract intergenic regions to the complement to get intronic regions only
    intronic_regions = complement.subtract(b=intergenic_bed)

    # Subtract to make sure all ncRNAs obtained before are not included in these regions
    final_intronic_regions = intronic_regions.subtract(b=ncRNA_bed).saveas(output)


    sp.call(f'rm all_genes_{species_name}*.bed temp_chr_size{species_name}.tsv', shell=True)
    sp.call(f'rm all_exons_{species_name}*.bed temp2_chr_size{species_name}.tsv', shell=True)


def get_exonic_regions(gtf_df, species_name, ncRNA_bed, output):
    """ Get the exonic regions in biotypes of interest (exonic regions that don't 
        overlap with ncRNAs (negative examples))."""
    # Select the exons of all genes that are not snoRNA/snRNA/tRNA
    gene_gtf = gtf_df[(gtf_df['feature'] == 'exon') & (~gtf_df['gene_biotype'].isin(['snoRNA', 'snRNA', 'tRNA']))]

    gene_gtf[bed_cols].to_csv(f'all_exons_real_{species_name}.bed', 
                        sep='\t', index=False, header=False)
    sp.call(f'sort -k1,1 -k2,2n all_exons_real_{species_name}.bed > all_exons_real_{species_name}_sorted.bed', 
            shell=True)
    gene_bed = BedTool(f'all_exons_real_{species_name}_sorted.bed')


    # Subtract to make sure all ncRNAs obtained before are not included in these regions
    final_exonic_regions = gene_bed.subtract(b=ncRNA_bed).saveas(output)

    sp.call(f'rm all_exons_real_{species_name}*.bed', shell=True)


def main():
    ncRNA_bed_final = get_bed_ncRNA(all_ncRNA_df, species)
    intergenic_bed = get_intergenic_regions(gtf, species, chr_size_file, 
                                ncRNA_bed_final, intergenic_output)
    get_intronic_regions(gtf, species, chr_size_file, intergenic_bed, 
                                ncRNA_bed_final, intronic_output)
    get_exonic_regions(gtf, species, ncRNA_bed_final, exonic_output)
    sp.call(f'rm all_ncRNA_{species}.*bed', shell=True)

main()