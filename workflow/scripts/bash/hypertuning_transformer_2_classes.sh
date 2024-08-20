#!/bin/bash

echo START bash

# Initialize variables
#hyperparams_space=$1
random_state=$1
pretrained_model=$2
x_tuning=$3
y_tuning=$4
output_best_hyperparams=$5

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
pip install optuna==3.1.0 --no-index
echo Activated env

# To fix bug happening on V100l (CUDA not available otherwise...)
nvidia-modprobe

python3 scripts/python/hypertuning_transformer_2_classes.py \
$random_state \
$pretrained_model \
$x_tuning \
$y_tuning \
$output_best_hyperparams 

echo Tuning completed!
