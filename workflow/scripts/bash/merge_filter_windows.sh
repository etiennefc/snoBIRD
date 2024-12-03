#!/bin/bash

# Initialize variables
preds=$1
fasta_dir=$2
fasta=$3
pretrained_model=$4
tokenizer=$5
model=$6
output_filtered=$7
output_centered=$8
fixed_length=$9
step_size=${10}
chunk_size=${11}
batch_size=${12}
num_labels=${13}
prob_thresh=${14}
cons_window=${15}
python_script=${16}
cluster_env=${17}
profile=${18}
gpu=${19}
final_output=${20}

# Load module
module load python

# Copy and activate virtualenv on computing node
cp $cluster_env $SLURM_TMPDIR/
tar -xzf $SLURM_TMPDIR/snoBIRD_env.tar.gz -C $SLURM_TMPDIR/
source $SLURM_TMPDIR/snoBIRD_env/bin/activate

# Copy all data on $SLURM_TMPDIR
cp $python_script $SLURM_TMPDIR/merge_filter_windows.py
cp scripts/python/utils.py $SLURM_TMPDIR/utils.py
cp --parents results/intermediate/predictions/first_model/positive_windows*tsv $SLURM_TMPDIR/
cp -r --parents $fasta_dir $SLURM_TMPDIR/
cp --parents $fasta $SLURM_TMPDIR/
cp -r $pretrained_model $SLURM_TMPDIR/
cp -r $tokenizer $SLURM_TMPDIR/
cp $model $SLURM_TMPDIR/snoBIRD_first_model.pt
cp -r $pretrained_model $SLURM_TMPDIR/
cp -r $tokenizer $SLURM_TMPDIR/


# To fix bug happening on V100l (CUDA not available otherwise...)
# Use only GPU if step_size > 1 (not requested otherwise, because all possible windows will 
# already have been predicted (i.e. no need to use a GPU)
if [ $step_size -gt 1 ]; then nvidia-modprobe; fi

python3 $SLURM_TMPDIR/merge_filter_windows.py \
$preds \
$SLURM_TMPDIR/$fasta_dir \
$SLURM_TMPDIR/$fasta \
$SLURM_TMPDIR/DNA_BERT_6_pretrained_model/ \
$SLURM_TMPDIR/DNA_BERT_6_tokenizer/ \
$SLURM_TMPDIR/snoBIRD_first_model.pt \
$output_filtered \
$output_centered \
$fixed_length \
$step_size \
$chunk_size \
$batch_size \
$num_labels \
$prob_thresh \
$cons_window \
$profile \
$SLURM_TMPDIR/ \
$gpu \
$final_output

echo "Merge and filter windows completed!"
