#!/usr/bin/env python3
import os
import sys
from math import ceil
from snakemake.utils import read_job_properties
import subprocess as sp
import json

with open("../config/config.json", 'r') as f:
    data = json.load(f)
#sp.call(f"echo {data}", shell=True)


# Define function to get the time needed for a given chr_size
def time_limit(chr_size_, gpu='A100'):
    if gpu == 'H100':
        rate = 1250000
    elif gpu == 'A100':
        rate = 550000 #~550KB/h
    elif gpu == 'V100':
        rate = 400000
    elif gpu == 'P100':
        rate = 200000
    if chr_size <= rate:
        time = "0-1:00:00"
    else:
        hours = ceil(chr_size_/rate)
        time = f"0-{hours}:00:00"
    return time


jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)
#print(config)
# Adapt time limit depending on the GPU generation
#if job_properties["rule"] in ["genome_prediction", "shap_snoBIRD", 
#                            "sno_pseudo_prediction"]:
#
#    # Update the cluster properties with the new time limit
#    inputs = job_properties.get("input", [])
#    chr_path = [i for i in inputs if i.endswith('.fa')][0]
#    chr_size = os.path.getsize(chr_path) # in bytes
#    job_properties['cluster']['time'] = time_limit(chr_size, gpu='V100')

cmdline = "sbatch "
for param, val in job_properties['cluster'].items():
    cmdline += "--{param} {val} ".format(param=param, val=val)

# Set up dependencies
dependencies = set(sys.argv[1:-1])
if dependencies:
    cmdline += " --dependency=afterok:{} ".format(":".join(dependencies))

# Adding the actual job
cmdline += jobscript

# remove the leading and trailing white space for the submitted jobid
cmdline += r" | awk '{print substr($NF, 0, length($NF))}'"

sys.stdout.write(cmdline)

os.system(cmdline)
