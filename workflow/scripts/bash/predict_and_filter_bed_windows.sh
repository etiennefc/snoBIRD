#!/bin/bash

# Initialize variables
model=$1
fasta=$2
bed=$3
pretrained_model=$4
tokenizer=$5
output_centered=$6
fixed_length=$7
python_script=$8
batch_size=$9
num_labels=${10}
prob_thresh=${11}
cluster_env=${12}
profile=${13}
gpu=${14}
final_output=${15}

# Load module
module load python

# Copy and activate virtualenv on computing node
cp $cluster_env $SLURM_TMPDIR/
tar -xzf $SLURM_TMPDIR/snoBIRD_env.tar.gz -C $SLURM_TMPDIR/
source $SLURM_TMPDIR/snoBIRD_env/bin/activate

# Copy all data on $SLURM_TMPDIR
cp $python_script $SLURM_TMPDIR/predict_and_filter_bed_windows.py
cp scripts/python/utils.py $SLURM_TMPDIR/utils.py
cp --parents $fasta $SLURM_TMPDIR/
cp --parents $bed $SLURM_TMPDIR/
cp -r $pretrained_model $SLURM_TMPDIR/
cp -r $tokenizer $SLURM_TMPDIR/
cp $model $SLURM_TMPDIR/snoBIRD_first_model.pt



# To fix bug happening on V100l (CUDA not available otherwise...)
nvidia-modprobe

python3 $SLURM_TMPDIR/predict_and_filter_bed_windows.py \
$SLURM_TMPDIR/$fasta \
$SLURM_TMPDIR/$bed \
$SLURM_TMPDIR/DNA_BERT_6_pretrained_model/ \
$SLURM_TMPDIR/DNA_BERT_6_tokenizer/ \
$SLURM_TMPDIR/snoBIRD_first_model.pt \
$output_centered \
$fixed_length \
$batch_size \
$num_labels \
$prob_thresh \
$profile \
$gpu \
$final_output 

echo "Predict on with SnoBIRD's first model and filter bed windows completed!"
