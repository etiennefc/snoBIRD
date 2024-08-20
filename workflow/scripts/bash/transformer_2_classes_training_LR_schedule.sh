#!/bin/bash

echo START bash

# Initialize variables
pretrained_model=$1
fold_num=$2
random_state=$3
x_train=$4
y_train=$5
best_hyperparams=$6
output_model=$7
output_loss=$8
output_f1=$9

# Load modules
module load StdEnv/2020
module load python/3.11.2

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate

# install packages
pip install torch==2.0.1 --no-index
pip install pandas==2.0.0 --no-index
pip install scikit_learn==1.3.0 --no-index
pip install numpy==1.24.2 --no-index
pip install transformers==4.31.0 --no-index
echo Activated env

# To fix bug happening on V100l (CUDA not available otherwise...)
nvidia-modprobe

python3 scripts/python/training_2_classes_transformer_LRschedule.py \
$pretrained_model \
$fold_num \
$random_state \
$x_train \
$y_train \
$best_hyperparams \
$output_model \
$output_loss \
$output_f1

echo Training completed
