#!/bin/bash

echo START bash

# Initialize variables
shap_df=$1
pred_df=$2
output=$3
fixed_length=$4
python_script=$5
min_box_dist=$6
flanking_nt=$7
cluster_env=$8

## Load module
module load python

# Copy and activate virtualenv on computing node
cp $cluster_env $SLURM_TMPDIR/
tar -xzf $SLURM_TMPDIR/snoBIRD_env.tar.gz -C $SLURM_TMPDIR/
source $SLURM_TMPDIR/snoBIRD_env/bin/activate

# Copy all data on $SLURM_TMPDIR
cp $python_script $SLURM_TMPDIR/find_sno_limits_shap.py
cp scripts/python/utils.py $SLURM_TMPDIR/utils.py
cp --parents $shap_df $SLURM_TMPDIR/
cp --parents $pred_df $SLURM_TMPDIR/


python3 $SLURM_TMPDIR/find_sno_limits_shap.py \
$SLURM_TMPDIR/$shap_df \
$SLURM_TMPDIR/$pred_df \
$output \
$fixed_length \
$min_box_dist \
$flanking_nt 



echo 'Finished finding snoRNA limits based on SHAP values!'
