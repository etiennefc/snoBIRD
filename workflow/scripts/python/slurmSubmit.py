#!/usr/bin/env python3
import os
import sys
from math import ceil
from snakemake.utils import read_job_properties
import subprocess as sp



# Define function to get the time needed for a given chr_size
def time_limit(chr_size_, step_size, gpu='A100'):
    # the bigger the chr_size_, the longer it will take to run
    # the bigger the step_size is, the faster it will take to run
    if gpu == 'H100':
        rate = 1250000 * step_size
    elif gpu == 'A100':
        rate = 550000 * step_size #~550KB/h
    elif gpu == 'V100':
        rate = 400000 * step_size
    elif gpu == 'P100':
        rate = 200000 * step_size
    if chr_size_ <= rate:
        time = "0-1:00:00"
    else:
        hours = ceil(chr_size_/rate)
        time = f"0-{hours}:00:00"
    return time


jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

# Adapt time limit depending on the GPU generation
if job_properties["rule"] in ["genome_prediction"]:
    gpu_type = job_properties['params']['gpu']
    step_size = job_properties['params']['step_size']
    # Update the cluster properties with the new time limit
    if gpu_type != 'Unknown':
        if job_properties['params']['chunks'] == True:
            # get chunk size in bytes
            chunk_size = int(job_properties['params']['chunk_size'])
            job_properties['cluster']['time'] = time_limit(chunk_size, 
                                                    step_size, gpu=gpu_type)
            
        else:
            chr_size_dict = job_properties['params']['chr_dict']
            chr_wildcard = job_properties['wildcards']['chr_']
            chr_wildcard_size = chr_size_dict[chr_wildcard]
            job_properties['cluster']['time'] = time_limit(chr_wildcard_size, 
                                                    step_size, gpu=gpu_type)



cmdline = "sbatch "
cmd_line = "sbatch --acount def-scottmic "
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
