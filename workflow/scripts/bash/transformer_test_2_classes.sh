#!/bin/bash

echo START bash

# Initialize variables
pretrained_model=$1
x_test=$2
y_test=$3
best_hyperparams=$4
model=$5
output_metrics=$6
output_preds=$7
python_script=$8

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

python3 $python_script \
$pretrained_model \
$x_test \
$y_test \
$best_hyperparams \
$model \
$output_metrics \
$output_preds

echo Test prediction completed
