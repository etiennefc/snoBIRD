""" Common functions used across the workflow"""

wildcard_constraints:
    simple_models = "({})".format("|".join(config["simple_models"]))

def get_species_genome(wildcards):
    # Get the fasta of the genome of a given species
    # The wildcard here is sno_fasta
    species_name = str(wildcards).split('_cd_')[0]
    species_name = species_name.split(' ')[0]
    protists = ["leishmania_major", "dictyostelium_discoideum", "giardia_lamblia"]
    fungi = ['saccharomyces_cerevisiae', 'schizosaccharomyces_pombe', 
            'aspergillus_fumigatus', 'neurospora_crassa', 'candida_albicans']
    if species_name == 'tetrahymena_thermophila':
        path = rules.download_tetrahymena_genome.output.genome
    elif species_name == 'ostreococcus_tauri':
        path = rules.download_o_tauri_genome.output.genome
    elif species_name in protists:
        path = expand(rules.download_other_protist_genome.output.genome, species=species_name)
    elif species_name in ['homo_sapiens', 'mus_musculus']:
        path = rules.download_mammal_genome.output.genome
    elif species_name in fungi:
        path = expand(rules.download_yeast_genome.output.genome, species=species_name)
    else:
        path = expand(rules.download_genome.output.genome, species=species_name)
    return path
    print(path)

def get_all_genomes(dir):
    # Get all the fasta files of genomes in a given directory
    return glob.glob(dir)

def get_chr_size(species):
    # Get the size of all chr of a given species
    sp = str(species).split(' ')[0]
    return glob.glob(f'data/references/chr_size/{sp}*_chr_size.tsv')[0]

def get_species_gtf(species):
    # Get the gtf of the genome of a given species
    species = str(species)
    protists = ["leishmania_major", "dictyostelium_discoideum", "giardia_lamblia"]
    animals = ['macaca_mulatta', 'ornithorhynchus_anatinus', 'gallus_gallus', 
                'caenorhabditis_elegans', 'drosophila_melanogaster']
    plants = ['arabidopsis_thaliana', 'oryza_sativa']
    fungi = ['saccharomyces_cerevisiae', 'schizosaccharomyces_pombe', 
            'aspergillus_fumigatus', 'neurospora_crassa', 'candida_albicans']
    if species in fungi:
        path = rules.download_yeast_gtf.output.gtf
    elif species == 'homo_sapiens':
        path = rules.download_human_gtf.output.gtf
    elif species == 'mus_musculus':
        path = rules.download_mouse_gtf.output.gtf
    elif species == 'tetrahymena_thermophila':
        path = rules.download_tetrahymena_gtf.output.gtf
    elif species in protists:
        path = rules.download_protist_gtf.output.gtf
    elif species in plants:
        path = rules.download_plant_gtf.output.gtf
    elif species in animals:
        path = rules.download_animal_gtf.output.gtf
    return path

def join_list(l, subl, remove=False):
    # From a list l, return a string of all items in subl joined by '|'
    small_list = [a for a in l if a in subl]
    # If we want to remove (instead of only keeping) items of subl from l
    if remove==True:
        small_list = [a for a in l if a not in subl]
    return "{}".format("|".join(small_list))

