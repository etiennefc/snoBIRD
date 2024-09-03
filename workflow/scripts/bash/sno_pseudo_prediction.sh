#!/bin/bash

echo START bash

# Initialize variables
pretrained_model=$1
tokenizer=$2
preds=$3
model=$4
output_preds=$5
fixed_length=$6
python_script=$7

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


python3 $python_script \
$pretrained_model \
$tokenizer \
$preds \
$model \
$fixed_length \
$output_preds

echo Predictions with SnoBIRD's second model (expressed vs pseudogene) completed!
