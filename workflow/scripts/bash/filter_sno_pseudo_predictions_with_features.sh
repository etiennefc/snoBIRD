#!/bin/bash

# Initialize variables
sno_limits=$1
pred_df=$2
output=$3
output_type=$4
fixed_length=$5
python_script=$6
prob_second_model=$7
box_score_thresh=$8
score_c_thresh=$9
score_d_thresh=${10}
terminal_stem_score_thresh=${11}
normalized_sno_stability_thresh=${12}
cluster_env=${13}

## Load module
module load python

# Copy and activate virtualenv on computing node
cp $cluster_env $SLURM_TMPDIR/
tar -xzf $SLURM_TMPDIR/snoBIRD_env.tar.gz -C $SLURM_TMPDIR/
source $SLURM_TMPDIR/snoBIRD_env/bin/activate

# Copy all data on $SLURM_TMPDIR
cp $python_script $SLURM_TMPDIR/filter_sno_pseudo_predictions_with_features.py
cp scripts/python/utils.py $SLURM_TMPDIR/utils.py
cp --parents $sno_limits $SLURM_TMPDIR/
cp --parents $pred_df $SLURM_TMPDIR/


python3 $SLURM_TMPDIR/filter_sno_pseudo_predictions_with_features.py \
$SLURM_TMPDIR/$sno_limits \
$SLURM_TMPDIR/$pred_df \
$output \
$output_type \
$fixed_length \
$prob_second_model \
$box_score_thresh \
$score_c_thresh \
$score_d_thresh \
$terminal_stem_score_thresh \
$normalized_sno_stability_thresh 



echo "SnoBIRD's run is fully completed!"
