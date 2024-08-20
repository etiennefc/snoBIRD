#!/usr/bin/env python3

"""
Adapted from https://bitbucket.org/snakemake/snakemake/issues/28/clustering-jobs-with-snakemake
Launch with :
snakemake -j 99 --use-conda --cluster-config cluster.json --immediate-submit --notemp --cluster 'python3 slurmSubmit.py {dependencies}'
"""
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

cmdline = "sbatch "
for param, val in job_properties['cluster'].items():
    cmdline += "--{param} {val} ".format(param=param, val=val)

# Check if the rule name is 'genome_windows_separate_chrom'
#if job_properties["rule"] == "genome_windows_separate_chrom":
#
#    # Update the cluster properties with the new time limit
#    inputs = job_properties.get("input", [])
#    chr_path = [i for i in inputs if i.endswith('.fa') | i.endswith('.fasta')][0]
#    chr_size = os.path.getsize(chr_path) # in bytes
#    job_properties['cluster']['time'] = time_limit(chr_size)


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
