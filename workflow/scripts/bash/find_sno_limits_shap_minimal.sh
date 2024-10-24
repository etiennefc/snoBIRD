#!/bin/bash

echo START bash

# Initialize variables
shap_df=$1
pred_df=$2
output=$3
output_type=$4
fixed_length=$5
python_script=$6
min_box_dist=$7
flanking_nt=$8
cluster_env=$9

## Load module
module load python

# Copy and activate virtualenv on computing node
cp $cluster_env $SLURM_TMPDIR/
tar -xzf $SLURM_TMPDIR/snoBIRD_env.tar.gz -C $SLURM_TMPDIR/
source $SLURM_TMPDIR/snoBIRD_env/bin/activate

# Copy all data on $SLURM_TMPDIR
cp $python_script $SLURM_TMPDIR/find_sno_limits_shap_minimal.py
cp scripts/python/utils.py $SLURM_TMPDIR/utils.py
cp --parents $shap_df $SLURM_TMPDIR/
cp --parents $pred_df $SLURM_TMPDIR/


python3 $SLURM_TMPDIR/find_sno_limits_shap_minimal.py \
$SLURM_TMPDIR/$shap_df \
$SLURM_TMPDIR/$pred_df \
$output \
$output_type \
$fixed_length \
$min_box_dist \
$flanking_nt 



echo "SnoBIRD's run is fully completed!"
