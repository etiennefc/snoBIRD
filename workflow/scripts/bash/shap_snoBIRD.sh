#!/bin/bash

echo START bash

# Initialize variables
model=$1
df=$2
pretrained_model=$3
tokenizer=$4
output=$5
fixed_length=$6
python_script=$7
batch_size=$8
num_labels=$9
cluster_env=${10}
profile=${11}

## Load module
module load python

# Copy and activate virtualenv on computing node
cp $cluster_env $SLURM_TMPDIR/
tar -xzf $SLURM_TMPDIR/snoBIRD_env.tar.gz -C $SLURM_TMPDIR/
source $SLURM_TMPDIR/snoBIRD_env/bin/activate

# Copy all data on $SLURM_TMPDIR
cp $python_script $SLURM_TMPDIR/shap_snoBIRD.py
cp scripts/python/utils.py $SLURM_TMPDIR/utils.py
cp $model $SLURM_TMPDIR/snoBIRD_first_model.pt
cp $df $SLURM_TMPDIR/filtered_center_positive_windows.bed
cp -r $pretrained_model $SLURM_TMPDIR/
cp -r $tokenizer $SLURM_TMPDIR/

# To fix bug happening on V100l (CUDA not available otherwise...)
nvidia-modprobe

python3 $SLURM_TMPDIR/shap_snoBIRD.py \
$SLURM_TMPDIR/snoBIRD_first_model.pt \
$SLURM_TMPDIR/filtered_center_positive_windows.bed \
$SLURM_TMPDIR/DNA_BERT_6_pretrained_model/ \
$SLURM_TMPDIR/DNA_BERT_6_tokenizer/ \
$output \
$fixed_length \
$batch_size \
$num_labels \
$profile



echo "SHAP values computation completed!"
