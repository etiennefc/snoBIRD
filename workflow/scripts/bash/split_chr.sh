#!/bin/bash

# Initialize variables
python_script=$1
input_fasta=$2
output_dir=$3
chunks=$4
chunk_size=$5
cluster_env=$6

# Load module
module load python

# Copy and activate virtualenv on computing node
cp $cluster_env $SLURM_TMPDIR/
tar -xzf $SLURM_TMPDIR/snoBIRD_env.tar.gz -C $SLURM_TMPDIR/
source $SLURM_TMPDIR/snoBIRD_env/bin/activate

# Run python script
python3 $python_script $input_fasta $output_dir $chunks $chunk_size 

echo "Fasta split completed!"
