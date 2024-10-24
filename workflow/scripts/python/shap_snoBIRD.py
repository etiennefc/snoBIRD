#!/usr/bin/python3
import sys
import pandas as pd
import subprocess as sp
import os
import numpy as np
from math import ceil
import time 
import torch
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification 
from transformers import TextClassificationPipeline
import shap
from utils import seq2kmer, kmer_score_per_nt, batch_generator, shap_batch
import warnings
warnings.filterwarnings("ignore")

# Define model path and other variables
model_path = str(sys.argv[1])
pretrained_model = str(sys.argv[3])
tokenizer_path = str(sys.argv[4])
output_df = str(sys.argv[5])
fixed_length = int(sys.argv[6])
batch_size = int(sys.argv[7])
num_labels = int(sys.argv[8])
profile = str(sys.argv[9])
cd_df = pd.read_csv(str(sys.argv[2]), names = ['chr', 'start', 'end', 'gene_id', 
                            'score', 'strand', 'block_id', 
                            f'extended_{fixed_length}nt_sequence'], sep='\t')


# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# Force to not parallelize tokenizing before dataloader (causes forking errors otherwise)
os.environ["TOKENIZERS_PARALLELISM"] = "true"

# Limit the number of threads that torch can spawn with (to avoid core oversubscription)
# i.e. set the number of threads to the number of CPUs requested (not all CPUs physically installed)
# Do only on cluster, not in local
if profile != 'local':
    N_CPUS = os.environ.get("SLURM_CPUS_PER_TASK")
    torch.set_num_threads(int(N_CPUS))
    
    # Allow TF32 on matrix multiplication to speed up computations
    torch.backends.cuda.matmul.allow_tf32 = True
    
    # Allow TF32 when using cuDNN library (GPU-related library usually automatically installed on the cluster)
    torch.backends.cudnn.allow_tf32 = True


# Show packages versions
sp.call(f'echo PANDAS VERSION: {pd.__version__}', shell=True)
sp.call(f'echo TORCH VERSION: {torch.__version__}', shell=True)
sp.call(f'echo NUMPY VERSION: {np.__version__}', shell=True)
sp.call(f'echo TRANSFORMERS VERSION: {transformers.__version__}', shell=True)
sp.call(f'echo IS CUDA AVAILABLE?: {torch.cuda.is_available()}', shell=True)

# Load model and tokenizer 
start_time = time.time()
tokenizer = AutoTokenizer.from_pretrained(tokenizer_path)
model = BertForSequenceClassification.from_pretrained(pretrained_model, 
        num_labels=num_labels)
model.load_state_dict(torch.load(model_path))
model.to(device)
model.classifier.to(device)
model.eval()
end_time = time.time()
sp.call(f'echo LOAD INITIAL MODEL: {end_time-start_time}s', shell=True)


# Load predicted CD sequences   
all_seqs = list(cd_df[f'extended_{fixed_length}nt_sequence'])
all_gene_ids = list(cd_df['gene_id'])
pipe2 = TextClassificationPipeline(model=model, tokenizer=tokenizer, 
                        device=device, batch_size=batch_size)


# Compute SHAP values for each prediction to find the importance of each nt in 
# the window for a given prediction
all_scores = []
index_i = 0
for i, batch in enumerate(batch_generator(all_seqs, batch_size)):
    batch_scores = shap_batch(batch, 
                            all_gene_ids[index_i:index_i + len(batch)], pipe2)
    index_i += len(batch)
    all_scores.extend(batch_scores)

df = pd.DataFrame(all_scores, 
            columns=['gene_id', 'predicted_label', 'probability', 'CLS'] + \
                    [f'SHAP_pos_{i}' for i in range(int(fixed_length))] + \
                    ['SEP'])

df.to_csv(output_df, sep='\t', index=False)
