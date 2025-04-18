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
        rate = 2500000 * step_size
    elif gpu == 'A100':
        rate = 1650000 * step_size 
    elif gpu == 'V100':
        rate = 450000 * step_size
    elif gpu == 'P100':
        rate = 200000 * step_size
    if chr_size_ <= rate:
        time = "0-1:00:00"
    else:
        hours = ceil(chr_size_/rate)
        time = f"0-{hours}:00:00"
    return time

# Define function to get the time needed for a given number of entries in a 
# input bed
def time_limit_bed(bed_entries, gpu='A100'):
    # the bigger the number of bed_entries, the longer it will take to run
    if gpu == 'H100':
        rate = 2500000  # predicted window per hour
    elif gpu == 'A100':
        rate = 1500000 
    elif gpu == 'V100':
        rate = 450000
    elif gpu == 'P100':
        rate = 200000
    if int(bed_entries) <= rate:
        time = "0-1:00:00"
    else:
        hours = ceil(int(bed_entries)/rate)
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
        # get chunk/chr size in bytes
        chunk_chr_size = int(job_properties['params']['real_chunk_chr_size'])
        job_properties['cluster']['time'] = time_limit(chunk_chr_size, 
                                                    step_size, gpu=gpu_type)

if job_properties["rule"] in ["predict_and_filter_bed_windows"]:
    gpu_type = job_properties['params']['gpu']
    nb_entries = job_properties['params']['bed_entries']
    # Update the cluster properties with the new time limit
    if gpu_type != 'Unknown':
        job_properties['cluster']['time'] = time_limit_bed(nb_entries, 
                                                gpu=gpu_type)


# No GPU needed for merge_filter_windows if step_size = 1
if job_properties["rule"] in ["merge_filter_windows"]:
    step_size = job_properties['params']['step_size']
    if step_size == 1:
        job_properties['cluster'].pop('gpus-per-node', None)

cmdline = "sbatch "
#cmdline = "sbatch --account=[def-your_account] "

# Change slurm log name for each job to include the job name
if job_properties["rule"] == 'shap_snoBIRD':
    job_name = 'shap_snoBIRD'
else:
    job_name = job_properties['cluster']['job-name']
# logs directory is already created by snoBIRD.py script
cmdline = cmdline + f"--output=logs/slurmLog-{job_name}-%j.out "
cmdline = cmdline + f"--error=logs/slurmLog-{job_name}-%j.out "

# Add cluster params value
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
