#!/usr/bin/python3
import sys
import pandas as pd
import subprocess as sp
from Bio import SeqIO
import os
import numpy as np
from math import ceil
import time 
import torch
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification, TextClassificationPipeline
from utils import seq2kmer, batch_generator
import warnings
warnings.filterwarnings("ignore")

model_path = str(sys.argv[1])
genome = str(sys.argv[2])
chr_name = genome.split('/')[-1].replace('.fasta', '').replace(
            '.fa', '').replace('.FASTA', '').replace('.FA', '')
pretrained_model = str(sys.argv[3])
tokenizer_path = str(sys.argv[4])
window_size = int(sys.argv[6])
step_size_defined = int(sys.argv[7])
strand = str(sys.argv[8])
batch_size = int(sys.argv[9])
num_labels = int(sys.argv[10])
profile = str(sys.argv[11])
output = str(sys.argv[5])

sp.call(f'echo {model_path} {genome} {chr_name} {pretrained_model} {tokenizer_path} {window_size} {step_size_defined} {strand} {output}', shell=True)

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
model = BertForSequenceClassification.from_pretrained(pretrained_model, num_labels=num_labels)
model.load_state_dict(torch.load(model_path)) 
model.to(device, non_blocking=True)
model.classifier.to(device, non_blocking=True)
model.eval()
end_time = time.time()
sp.call(f'echo LOAD INITIAL MODEL: {end_time-start_time}s', shell=True)


if strand == 'plus':
    chr_dict = {record.id: str(record.seq) 
                for record in SeqIO.parse(genome, "fasta")}
elif strand == 'minus':
    chr_dict = {record.id: str(record.seq.reverse_complement()) 
                for record in SeqIO.parse(genome, "fasta")}
else: # predict on both strands
    chr_dict = {record.id: str(record.seq) 
                for record in SeqIO.parse(genome, "fasta")}
    chr_dict_neg = {record.id: str(record.seq.reverse_complement()) 
                for record in SeqIO.parse(genome, "fasta")}



total_window = ceil(((len(chr_dict[chr_name]) - int(window_size))/int(step_size_defined)) + 1)
sp.call(f'echo TOTAL WINDOWS {total_window}', shell=True)


# Create a generator which yields one batch of consecutive sequence at a time
def scan_fasta(chr_dict_, chr_seq, window_size, mini_batch_size, step_size=step_size_defined, mem_batch_size=10):
    id_chr = [chr_id for chr_id, s in chr_dict_.items() if s == chr_seq][0]
    len_chr = len(chr_seq)
    if len_chr > window_size:
        # Create big batches of sequences across all the chr length
        for i in range(0, int(len_chr - window_size + 1), int(step_size_defined*mini_batch_size*mem_batch_size)):
            #print(f'{id_chr}: {i/int(len_chr - window_size / step_size + 1)*100}')
            big_batch, starts, ends = [], [], []
            # Get sequence/start/end of each window in big batch
            for j in range(i, i+step_size_defined*mini_batch_size*mem_batch_size, step_size_defined):
                window = chr_seq[j:j+window_size]
                if len(window) == window_size: # all windows except the last of the batch
                    starts.append(j + 1)
                    ends.append(j + 1 + window_size)
                    big_batch.append(window)
                else:  # if last window doesn't end flush with the end of big_batch (i.e. smaller), change last window to fit with window_size length
                    window_diffe = window_size - len(window)
                    window = chr_seq[j-window_diffe:j-window_diffe+window_size]
                    starts.append(j-window_diffe + 1)
                    ends.append(j-window_diffe + 1 + window_size)
                    big_batch.append(window)
                    break
            #print(big_batch)  # last big_batch might be smaller than mini_batch_size*mem_batch_size
            kmer_seqs = [seq2kmer(seq, 6) for seq in big_batch]
            eval_dataset = tokenizer(kmer_seqs, return_tensors='pt', padding=True)
            #eval_dataset = TensorDataset(eval_dataset.input_ids, eval_dataset.attention_mask)
            #eval_dataloader = DataLoader(eval_dataset, batch_size=batch_size)

            # Convert starts/ends lists to tensor
            starts = torch.tensor(starts).to(device, non_blocking=True)
            ends = torch.tensor(ends).to(device, non_blocking=True)

            # where len(kmer_seqs) = mini_batch_size * mem_batch_size 
            yield eval_dataset, starts, ends   

    elif len_chr < window_size:
        ##TO ADAPT for just the snoRNA sequence (to pad to window_length?)
        raise ValueError(f'Length of sequence {id_chr} ({len_chr}) must be >= than length of window size ({window_size})')




df_cols = ['chr', 'strand', 'start', 'end', 'probs']

def predict(chr_dictio, chr_seq, window_size, mini_batch_size, strand_name, mem_batch=10):
    all_positives = []
    n_batch = 0
    for big_batch in scan_fasta(chr_dictio, chr_seq, window_size, mini_batch_size, mem_batch_size=mem_batch):
        ev_preds = []
        sp.call(f'echo EVAL BATCH {n_batch} {chr_name} {round(n_batch/(total_window/(mem_batch*mini_batch_size))*100, 2)}%', shell=True)
        n_batch += 1
        eval_dataset_real = TensorDataset(big_batch[0].input_ids, big_batch[0].attention_mask)
        for i, ev_batch in enumerate(DataLoader(eval_dataset_real, batch_size=batch_size, num_workers=2, pin_memory=True)):
            ev_input_ids, ev_attention_mask = ev_batch
            ev_input_ids = ev_input_ids.to(device, non_blocking=True)
            ev_attention_mask = ev_attention_mask.to(device, non_blocking=True)
            with torch.no_grad():  # nor gradient computation
                s_time = time.time()
                ev_input_ids_np_ = ev_input_ids.cpu().numpy()
                #print(ev_input_ids_np, len(ev_input_ids_np))
                if len(ev_input_ids_np_) != mini_batch_size:  # to account for the last batch at the end of chromosome
                    start_time = time.time()
                    ev_attention_mask_np_ = ev_attention_mask.cpu().numpy()
                    len_diff = mini_batch_size - len(ev_input_ids_np_)
                    sp.call('echo LAST BATCH MISMATCH', shell=True)
                    # Pad it on the right side to be same size as mini_batch_size
                    original_len = len(ev_input_ids_np_)
                    ev_input_ids_np2 = torch.tensor(np.pad(ev_input_ids_np_, ((0, len_diff), (0, 0)), 'constant', constant_values=0)).to(device, non_blocking=True)
                    ev_attention_mask_np2 = torch.tensor(np.pad(ev_attention_mask_np_, ((0, len_diff), (0, 0)), 'constant', constant_values=0)).to(device, non_blocking=True)
                    #sp.call(f'echo {type(ev_input_ids_np2)} {type(ev_attention_mask_np2)}', shell=True)
                    #sp.call(f'echo {ev_input_ids_np2.size()} {ev_attention_mask_np2.size()}', shell=True)
                    #sp.call(f'echo {ev_input_ids_np2}', shell=True)
                    #sp.call(f'echo {ev_attention_mask_np2}', shell=True)
                    outputs = model(ev_input_ids_np2, attention_mask=ev_attention_mask_np2)
                    # Get only the prediction for the relevant examples, not the padding
                    outputs2 = outputs[0][0:original_len]
                    probabilities = torch.softmax(outputs2, dim=1).to(device, non_blocking=True)
                    pred_labels = torch.argmax(probabilities, dim=1).to(device, non_blocking=True)
                    end_time = time.time()
                    sp.call(f'echo batch mismatch whole time: {end_time -start_time}s', shell=True)
                    if 1 in pred_labels:  # i.e. positive predictions
                        positives_index = (pred_labels == 1).nonzero().squeeze().to(device, non_blocking=True)
                        starts = big_batch[1][i*mini_batch_size:i*mini_batch_size+mini_batch_size].to(device, non_blocking=True)
                        pos_starts = starts[positives_index].to(device, non_blocking=True)
                        if pos_starts.dim() == 0:
                            pos_starts = pos_starts.unsqueeze(0)
                        ends = big_batch[2][i*mini_batch_size:i*mini_batch_size+mini_batch_size].to(device, non_blocking=True)
                        pos_ends = ends[positives_index].to(device, non_blocking=True)
                        if pos_ends.dim() == 0:
                            pos_ends = pos_ends.unsqueeze(0)

                        if probabilities[positives_index].dim() == 1:
                            probs = probabilities[positives_index].unsqueeze(0)[:, -1]
                        else:
                            probs = probabilities[positives_index][:, -1]
                            
                        l = []
                        for k, st in enumerate(pos_starts.tolist()):
                            l.append([chr_name, strand_name, st, pos_ends[k].item(), probs[k].item()])
                        all_positives.extend(l)
                else:  # for all other normal-sized batches
                    start_time = time.time()
                    outputs = model(ev_input_ids, attention_mask=ev_attention_mask)
                    end_time = time.time()
                    #print(f'run: {end_time -start_time}s')
                    #print(outputs)
                    start_time = time.time()
                    probabilities = torch.softmax(outputs.logits, dim=1).to(device, non_blocking=True)
                    pred_labels = torch.argmax(probabilities, dim=1).to(device, non_blocking=True)
                    end_time = time.time()
                    #print(f'conv: {end_time -start_time}s')
                    if 1 in pred_labels:  # i.e. positive predictions
                        #print(pred_labels)
                        positives_index = (pred_labels == 1).nonzero().squeeze().to(device, non_blocking=True)
                        #print(positives_index)
                        starts = big_batch[1][i*mini_batch_size:i*mini_batch_size+mini_batch_size].to(device, non_blocking=True)
                        #print(starts)
                        pos_starts = starts[positives_index].to(device, non_blocking=True)
                        if pos_starts.dim() == 0:
                            pos_starts = pos_starts.unsqueeze(0)
                        #print(pos_starts)
                        #print('STARTSSSSS', pos_starts, len(pos_starts))
                        ends = big_batch[2][i*mini_batch_size:i*mini_batch_size+mini_batch_size].to(device, non_blocking=True)
                        pos_ends = ends[positives_index].to(device, non_blocking=True)
                        if pos_ends.dim() == 0:
                            pos_ends = pos_ends.unsqueeze(0)
                        #print(pos_ends)
                        #print('ENDSSSS', pos_ends, len(pos_ends))
                        if probabilities[positives_index].dim() == 1:
                            probs = probabilities[positives_index].unsqueeze(0)[:, -1]
                        else:
                            probs = probabilities[positives_index][:, -1]
                        l = []
                        for k, st in enumerate(pos_starts.tolist()):
                            l.append([chr_name, strand_name, st, pos_ends[k].item(), probs[k].item()])
                        all_positives.extend(l)
                n_time = time.time()
                sp.call(f'echo pred time: {n_time -s_time}s', shell=True)
    prediction_df = pd.DataFrame(all_positives, columns=df_cols)
    #print(prediction_df)
    return prediction_df


if strand != 'both':
    strand_symbol = {'plus': '+', 'minus': '-'}
    s = strand_symbol[strand]
    start_time = time.time()
    seq_ = chr_dict[chr_name]
    sp.call(f'echo PREDICT ON {s} STRAND {chr_name}', shell=True)
    if len(seq_) <= window_size:
        kmer_seq = [seq2kmer(seq_, 6)]
        pipe = TextClassificationPipeline(model=model, tokenizer=tokenizer, device=device, batch_size=batch_size)
        tokenizer_kwargs = {'padding': True}
        pred_batch = pipe(kmer_seq, **tokenizer_kwargs)
        pred_label = [pred['label'] for pred in pred_batch if pred['label'] == 'LABEL_1']
        probability = [pred['score'] for pred in pred_batch][0]
        if len(pred_label) > 0:
            results_df = pd.DataFrame([[chr_name, s, 0, window_size, probability]], columns=df_cols)
        else:
            results_df = pd.DataFrame(columns=df_cols)
    else:
        results_df = predict(chr_dict, seq_, window_size, batch_size, s, mem_batch=10)
        results_df = results_df.sort_values(by=['start', 'end'])
    if strand == 'minus':
        # Correct for the actual start and ends of snoRNAs based on the first nt not the last
        results_df['start'] = len(seq_) - results_df['start'] + 1
        results_df['end'] = len(seq_) - results_df['end'] + 1
        # Switch start and end because it is the opposite on the - strand
        results_df = results_df.rename(columns={'start': 'end', 'end': 'start'})
        results_df = results_df[df_cols].sort_values(by=['start', 'end']).reset_index(drop=True)
    end_time = time.time()
    sp.call(f'echo FINAL elapsed time {chr_name}: {end_time -start_time}s', shell=True)
    results_df.to_csv(output, index=False, sep='\t')
else:  # predict on both strands
    start_time = time.time()
    seq_pos = chr_dict[chr_name]
    sp.call(f'echo PREDICT ON + STRAND {chr_name}', shell=True)
    if len(seq_pos) <= window_size:
        kmer_seq = [seq2kmer(seq_pos, 6)]
        pipe = TextClassificationPipeline(model=model, tokenizer=tokenizer, device=device, batch_size=batch_size)
        tokenizer_kwargs = {'padding': True}
        pred_batch = pipe(kmer_seq, **tokenizer_kwargs)
        pred_label = [pred['label'] for pred in pred_batch if pred['label'] == 'LABEL_1']
        probability = [pred['score'] for pred in pred_batch][0]
        if len(pred_label) > 0:
            pos_strand_results = pd.DataFrame([[chr_name, '+', 0, window_size, probability]], columns=df_cols)
        else:
            pos_strand_results = pd.DataFrame(columns=df_cols)
    else:
        pos_strand_results = predict(chr_dict, seq_pos, window_size, batch_size, '+', mem_batch=10)

    seq_neg = chr_dict_neg[chr_name]
    sp.call(f'echo PREDICT ON - STRAND {chr_name}', shell=True)
    if len(seq_neg) <= window_size:
        kmer_seq = [seq2kmer(seq_neg, 6)]
        pipe = TextClassificationPipeline(model=model, tokenizer=tokenizer, device=device, batch_size=batch_size)
        tokenizer_kwargs = {'padding': True}
        pred_batch = pipe(kmer_seq, **tokenizer_kwargs)
        pred_label = [pred['label'] for pred in pred_batch if pred['label'] == 'LABEL_1']
        probability = [pred['score'] for pred in pred_batch][0]
        if len(pred_label) > 0:
            neg_strand_results = pd.DataFrame([[chr_name, '-', 0, window_size, probability]], columns=df_cols)
        else:
            neg_strand_results = pd.DataFrame(columns=df_cols)
    else:
        neg_strand_results = predict(chr_dict_neg, seq_neg, window_size, batch_size, '-')

    # Correct for the actual start and ends of snoRNAs based on the first nt not the last (for - strand only)
    neg_strand_results['start'] = len(seq_neg) - neg_strand_results['start'] + 1
    neg_strand_results['end'] = len(seq_neg) - neg_strand_results['end'] + 1
    # Switch start and end because it is the opposite on the - strand
    neg_strand_results = neg_strand_results.rename(columns={'start': 'end', 'end': 'start'})
    neg_strand_results = neg_strand_results[df_cols]
    results_df = pd.concat([pos_strand_results, neg_strand_results]).sort_values(by=['start', 'end']).reset_index(drop=True)
    end_time = time.time()
    sp.call(f'echo FINAL elapsed time {chr_name}: {end_time -start_time}s', shell=True)
    results_df.to_csv(output, index=False, sep='\t')


