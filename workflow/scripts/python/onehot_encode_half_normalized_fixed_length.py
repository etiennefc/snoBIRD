#!/usr/bin/python3
import pandas as pd


tuning = pd.read_csv(snakemake.input.tuning, sep='\t')
training = pd.read_csv(snakemake.input.training, sep='\t')
test = pd.read_csv(snakemake.input.test, sep='\t')
fixed_length = snakemake.wildcards.fixed_length
one_hot_output = [path for path in snakemake.output 
                    if 'one_hot' in path]
target_output = [path for path in snakemake.output 
                    if 'target' in path]

# One-hot encode the extended sequences (5 possibilities: A/T/C/G/N)
def encode(seq_, ref):
    ind = [ref.find(nt) for nt in seq_]  # 0-indexed position of nt with regards 
                                         # to ref string of ATCGN
    encoded = []
    for i in ind:
        l = [0] * 5  # create [0,0,0,0,0]
        l[i] = 1 # Add the 1 at the right nt position
        encoded.append(l)
    return [nb for subl in encoded for nb in subl]


ref_nt = 'ATCGN'
new_cols = [nt_+'_'+str(nb) 
            for nb in range(1, int(fixed_length)+1) 
            for nt_ in ref_nt]  # all the new one-hot encoded column names (A_1, T_1, ..., G_211, N_211)

for i, df in enumerate([tuning, training, test]):
    seq_dict = dict(zip(df.gene_id, df[f'extended_{fixed_length}nt_sequence']))
    encoded_seqs = []
    # One-hot encode the sequence of each example
    for gene_id, seq in seq_dict.items():
        encoded_seqs.append(encode(seq, ref_nt))
    temp_df = pd.DataFrame(encoded_seqs, columns=new_cols)
    df = pd.concat([df, temp_df], axis=1)  # Add new cols horizontally

    # Create target df and create three classes to predict (0: other, 1: CD_snoRNA_pseudogene, 2: expressed_CD_snoRNA)
    target_df = df.copy()
    target_df['target'] = target_df['target'].replace({'other': 0, 'CD_snoRNA_pseudogene':1, 'expressed_CD_snoRNA': 2})
    target_df = target_df[['gene_id', 'target']]
    target_df.to_csv(target_output[i], sep='\t', index=False)


    # Scale/normalize the added features only using standardization
    df_num = df.copy()
    cols = ['box_score', 'structure_mfe', 'terminal_stem_mfe', 'length']
    df_num = df_num[cols]
    num_cols = list(df_num.columns)
    col_values = [df.drop(columns=cols + ['gene_biotype', 'species_name', f'extended_{fixed_length}nt_sequence', 'target'])]
    for col in num_cols:
        mean = df_num[col].mean()
        std = df_num[col].std()
        if std != 0:
            col_values.append(pd.DataFrame((df[col] - mean) / std).rename(
                                columns={col: f'{col}_norm'}))
        else:  # to deal with column that has all the same value, thus a std=0
            col_values.append(pd.DataFrame(df[col]).rename(  # we don't scale, but these values will either be all 0 or all 1
                            columns={col: f'{col}_norm'}))  
    scaled_df = pd.concat(col_values, axis=1)
    scaled_df.to_csv(one_hot_output[i], index=False, sep='\t')  # save one-hot-encoded scaled output



