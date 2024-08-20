#!/usr/bin/python3
import pandas as pd
import subprocess as sp

cols = ['gene_id', 'chr', 'strand', 'start', 'end', 'sequence', 'extended_sequence', 'species_name']
sno_literature = pd.read_csv(snakemake.input.sno_literature, sep='\t')
sno_literature = sno_literature[cols]

# Filter out Drosophila C/D found in Huang et al. 2005 to keep only the 
# Drosophila CD found in Sklias et al. 2024 (remove duplicates)
sno_literature = sno_literature[sno_literature['species_name'] != 'D_melanogaster']

# Load Drosophila expressed C/D (in TGIRT-Seq)
droso = pd.read_csv(snakemake.input.sno_tgirt_droso[0], sep='\t')
droso = droso[cols]

# Load human, mouse and S. cerevisiae dfs (C/D expressed in TGIRT-Seq)
sno_dfs = [sno_literature, droso]
for path in snakemake.input.sno_tgirt:
    df = pd.read_csv(path, sep='\t')
    sp_name = path.split('/')[-1].split('_expressed_')[0]
    df['species_name'] = sp_name
    df = df[cols]
    sno_dfs.append(df)

# Concat all cd dfs together
all_cd_df = pd.concat(sno_dfs)
all_cd_df.to_csv(snakemake.output.df, sep='\t', index=False)

# Create fasta of all expressed C/D sno sequences 
d = dict(zip(all_cd_df.gene_id, all_cd_df.sequence))

with open(snakemake.output.fa, 'w') as f:
    for id, sequence in d.items():
        f.write(f'>{id}\n')
        f.write(f'{sequence}\n')

