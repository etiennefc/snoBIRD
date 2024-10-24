#!/bin/bash

echo START bash

# Initialize variables
model=$1
pretrained_model=$2
tokenizer=$3
preds=$4
output_preds=$5
fixed_length=$6
python_script=$7
batch_size=$8
num_labels=$9
cluster_env=${10}
profile=${11}

# Load module
module load python

# Copy and activate virtualenv on computing node
cp $cluster_env $SLURM_TMPDIR/
tar -xzf $SLURM_TMPDIR/snoBIRD_env.tar.gz -C $SLURM_TMPDIR/
source $SLURM_TMPDIR/snoBIRD_env/bin/activate

# Copy all data on $SLURM_TMPDIR
cp $python_script $SLURM_TMPDIR/sno_pseudo_prediction.py
cp scripts/python/utils.py $SLURM_TMPDIR/utils.py
cp $model $SLURM_TMPDIR/snoBIRD_second_model.pt
cp -r $pretrained_model $SLURM_TMPDIR/
cp -r $tokenizer $SLURM_TMPDIR/
cp --parents $preds $SLURM_TMPDIR/

# To fix bug happening on V100l (CUDA not available otherwise...)
nvidia-modprobe


python3 $SLURM_TMPDIR/sno_pseudo_prediction.py \
$SLURM_TMPDIR/snoBIRD_second_model.pt \
$SLURM_TMPDIR/DNA_BERT_6_pretrained_model/ \
$SLURM_TMPDIR/DNA_BERT_6_tokenizer/ \
$SLURM_TMPDIR/$preds \
$output_preds \
$fixed_length \
$batch_size \
$num_labels \
$profile

echo "Predictions with SnoBIRD second model (expressed vs pseudogene) completed!"
