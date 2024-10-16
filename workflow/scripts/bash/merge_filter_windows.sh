#!/bin/bash

echo START merge and filter

# Initialize variables
pretrained_model=$1
tokenizer=$2
preds=$3
model=$4
fasta_dir=$5
genome=$6
output_preds=$7
output_center_preds=$8
fixed_length=$9
step_size=${10}
chunk_size=${11}
batch_size=${12}
num_labels=${13}
prob_threshold=${14}
min_cons_windows=${15}
python_script=${16}


# Load modules
module load StdEnv/2020
module load python/3.11.2

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate

# install packages
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
$pretrained_model \
$tokenizer \
$preds \
$model \
$fasta_dir \
$genome \
$output_preds \
$output_center_preds \
$fixed_length \
$step_size \
$chunk_size \
$batch_size \
$num_labels \
$prob_threshold \
$min_cons_windows

echo Merge and filter completed!
