#!/bin/bash

# Initialize variables
model=$1
genome=$2
pretrained_model=$3
tokenizer=$4
output=$5
fixed_length=$6
step_size=$7
strand=$8
python_script=$9
batch_size=${10}
num_labels=${11}
cluster_env=${12}
profile=${13}

# Load module
module load python

# Copy and activate virtualenv on computing node
cp $cluster_env $SLURM_TMPDIR/
tar -xzf $SLURM_TMPDIR/snoBIRD_env.tar.gz -C $SLURM_TMPDIR/
source $SLURM_TMPDIR/snoBIRD_env/bin/activate

# Copy all data on $SLURM_TMPDIR
cp $python_script $SLURM_TMPDIR/genome_prediction.py
cp scripts/python/utils.py $SLURM_TMPDIR/utils.py
cp $model $SLURM_TMPDIR/snoBIRD_first_model.pt
cp --parents $genome $SLURM_TMPDIR/
cp -r $pretrained_model $SLURM_TMPDIR/
cp -r $tokenizer $SLURM_TMPDIR/


# To fix bug happening on V100l (CUDA not available otherwise...)
nvidia-modprobe

python3 $SLURM_TMPDIR/genome_prediction.py \
$SLURM_TMPDIR/snoBIRD_first_model.pt \
$SLURM_TMPDIR/$genome \
$SLURM_TMPDIR/DNA_BERT_6_pretrained_model/ \
$SLURM_TMPDIR/DNA_BERT_6_tokenizer/ \
$output \
$fixed_length \
$step_size \
$strand \
$batch_size \
$num_labels \
$profile

echo "Genome-wide prediction completed for SnoBIRD's first model!"
