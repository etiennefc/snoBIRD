#!/bin/bash

echo START bash

# Initialize variables
model=$1
df=$2
pretrained_model=$3
tokenizer=$4
fixed_length=$5
python_script=$6
output=$7
batch_size=$8
num_labels=$9

## Load modules
module load StdEnv/2020
module load python/3.11.2


virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate

## Install packages
pip install torch --no-index
pip install pandas --no-index
pip install scikit_learn --no-index
pip install numpy --no-index
pip install transformers --no-index
pip install biopython --no-index
pip install shap --no-index
echo Activated env


# To fix bug happening on V100l (CUDA not available otherwise...)
nvidia-modprobe


python3 $python_script \
$model \
$df \
$pretrained_model \
$tokenizer \
$fixed_length \
$output \
$batch_size \
$num_labels

echo SHAP values computation completed!
