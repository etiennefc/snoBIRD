#!/usr/bin/python3
import pandas as pd


tuning = pd.read_csv(snakemake.input.tuning, sep='\t')
training = pd.read_csv(snakemake.input.training, sep='\t')
test = pd.read_csv(snakemake.input.test, sep='\t')

print(tuning, training, test)