#!/usr/bin/python3
import time
st_time = time.time()
import glob
import subprocess as sp
import pandas as pd


input_rDNA = snakemake.input.target_rDNA
test_fa = snakemake.input.test_fa
output = snakemake.output.predictions
time_output = snakemake.output.time_predictions

# Run snoscan
sp.call(f"snoscan {input_rDNA} {test_fa} -o {output}", shell=True)
end_time = time.time()
sp.call(f"echo {end_time - st_time}s Ca22chrM_C_albicans_SC5314 > {time_output}", shell=True)
