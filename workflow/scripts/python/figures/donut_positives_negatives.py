#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
import collections as coll

#positives = pd.read_csv(snakemake.input.positives, sep='\t')
positives = pd.read_csv('~/Desktop/Etienne/cd_predictor/workflow/data/references/data_augmentation/cd_rfam_filtered_all_sno_pseudo_fixed_length_194nt.tsv', sep='\t')
#negatives = pd.concat([pd.read_csv(path, sep='\t') for path in snakemake.input.negatives])
negatives = pd.read_csv('~/Desktop/Etienne/cd_predictor/workflow/data/references/negatives/data_augmentation/all_negatives_wo_pseudo_1_ratio_fixed_length_194nt.tsv', sep='\t')
pie_species = snakemake.output.pie_species
pie_neg_type = snakemake.output.pie_neg_type
sp_colors = snakemake.params.species_colors
ty_colors = snakemake.params.biotype_colors

# Convert species short name to long name
sp_name = snakemake.params.species_name
positives['species_name'] = positives['species_name'].replace(sp_name)

# Remove snoRNA pseudogenes from negatives
negatives = negatives[negatives['gene_biotype'] != 'snoRNA_pseudogene']

# Count the number of examples that are part of a given species/gene_biotype
species = ['homo_sapiens', 'mus_musculus', 'drosophila_melanogaster', 'ornithorhynchus_anatinus', 
            'caenorhabditis_elegans', 'macaca_mulatta', 'gallus_gallus', 
            'tetrahymena_thermophila', 'dictyostelium_discoideum', 'giardia_lamblia', 'leishmania_major', 
            'neurospora_crassa', 'saccharomyces_cerevisiae', 'candida_albicans', 'aspergillus_fumigatus',
            'arabidopsis_thaliana', 'oryza_sativa', 'ostreococcus_tauri']
biotype = ['random_exonic_region', 'random_intronic_region', 'random_intergenic_region',
            'tRNA', 'HACA_snoRNA', 'snRNA', 'pre_miRNA', 'shuffled_expressed_CD_snoRNA']
print(coll.Counter(positives.species_name))
print(coll.Counter(negatives.gene_biotype))
species_colors = [sp_colors[sp] for sp in species]
biotype_colors = [ty_colors[biot] for biot in biotype]
species_count, neg_type_count = [], []

for sp in species:
    nb = len(positives[positives['species_name'] == sp])
    species_count.append(nb)
for ty in biotype:
    nb2 = len(negatives[negatives['gene_biotype'] == ty])
    neg_type_count.append(nb2)

# Create pie chart for species distribution in FP
ft.donut(species_count, species, species_colors, '', 
                '', '', pie_species)
# Create pie chart for gene_biotype distribution in FP
ft.donut(neg_type_count, biotype, biotype_colors, '', 
                '', '', pie_neg_type)