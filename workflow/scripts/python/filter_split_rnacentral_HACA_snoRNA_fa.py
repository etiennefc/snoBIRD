#!/usr/bin/python3
import pandas as pd

input_fa = snakemake.input.fa 
outputs = snakemake.output.fa 

species = [sp.split('/')[-1].split('_HACA')[0] for sp in outputs]


# Get rnacentral and sequence of H/ACA in dict per species
d = {sp: {} for sp in species}
with open(input_fa, 'r') as f:
    for line in f:
        if line.startswith('>'):
            strain = ''
            rnacentral_id, other = line.split(' ', maxsplit=1)
            rnacentral_id = rnacentral_id.replace('>', '')
            species_genus, species_species = other.split(' ')[0:2]
            species_name = f'{species_genus}_{species_species}'.lower()
            if species_name == 'saccharomyces_cerevisiae':
                if 'AWRI796' not in other:  # keep only this yeast strain, otherwise too many duplicates
                    strain = 'bad_yeast_strain'
                    continue
            if species_name not in d.keys():  # exclude species other than our species of interest that are in the rnacentral fasta
                continue
            if rnacentral_id not in d[species_name].keys():
                d[species_name][rnacentral_id] = ''
        else:
            if (species_name not in d.keys()) | (strain == 'bad_yeast_strain'):  # skip sequence lines if it's not our species of interest
                continue  
            seq = line.replace('\n', '')
            d[species_name][rnacentral_id] += seq

# Remove duplicate and generate 1 H/ACA fasta per species
for sp in species:
    dictio = d[sp]
    df = pd.DataFrame(dictio.items(), columns=['id', 'seq'])
    df = df.drop_duplicates(subset=['seq'])
    final_dictio = dict(zip(df.id, df.seq))
    path = [output for output in outputs if sp in output][0]
    with open(path, 'w') as o:
        for k, v in final_dictio.items():
            o.write(f'>{k}\n{v}\n')



