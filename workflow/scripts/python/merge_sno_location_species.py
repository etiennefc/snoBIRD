#!/usr/bin/python3
import pandas as pd

df_paths = list(snakemake.input.dfs)
dfs = []
for path in df_paths:
    df = pd.read_csv(path, sep='\t')
    dfs.append(df)

final_df = pd.concat(dfs)
final_df.to_csv(snakemake.output.df, sep='\t', index=False)
