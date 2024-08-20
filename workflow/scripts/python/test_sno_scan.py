#!/usr/bin/python3
import glob
import subprocess as sp
import pandas as pd

input_rDNA = snakemake.input.target_rDNA
test_set = pd.read_csv(snakemake.input.test_set, sep='\t')
output = snakemake.output.predictions

# Split rDNA of all species in separate fasta files
sp.call("csplit -sz -f temp_rDNA_ "+input_rDNA+" '/>/' '{*}'", shell=True)

for species in pd.unique(test_set.species_name):
    print(species)
    species_rDNA = []
    for fa in glob.glob('temp_rDNA_*'):
        with open(fa, 'r') as f:
            for i, line in enumerate(f):
                if species in line:
                    species_rDNA.append(fa)
    # For a given species, concat all of its rDNA fa into 1 fa (ex: 5.8S, 18S and 28S)
    sp.call('cat '+' '.join(species_rDNA)+' > temp_rDNA_'+species+'.fa', shell=True)
    
    # Select examples in test set of given species and create fasta of these examples
    sp.call("""awk 'NR>1 && $3 == \""""+species+"""\" {print ">"$1" "$2" "$3"\\n"$4}' """+snakemake.input.test_set+""" > testt_"""+species, shell=True)

    # Run snoscan on each species 
    sp.call("snoscan temp_rDNA_"+species+".fa testt_"+species+" -o output__"+species, shell=True)

# Concat all predictions into 1 output
sp.call('cat output__* > '+output, shell=True)

sp.call('rm temp_rDNA_* testt_* output__*', shell=True)