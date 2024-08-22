#!/bin/bash

echo START bash

# Initialize variables
model=$1
genome=$2
output=$3
pretrained_model=$4
fixed_length=$5
step_size=$6
strand=$7
python_script=$8

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
#pip list
echo Activated env

# To fix bug happening on V100l (CUDA not available otherwise...)
nvidia-modprobe

python3 $python_script \
$model \
$genome \
$output \
$pretrained_model \
$fixed_length \
$step_size \
$strand \
$python_script

echo Prediction completed!
