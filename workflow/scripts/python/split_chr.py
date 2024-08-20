#!/usr/bin/python3
from Bio import SeqIO
import os
import subprocess as sp

### Set a flag to disable chunk formation

def multi_fasta(input_fasta, output_dir):
    # If multiple sequences in the input fasta file, split in separate 
    # fasta files per chr
    entries = list(SeqIO.parse(input_fasta, "fasta"))
    if len(entries) > 1:
        os.makedirs(f"{output_dir}fasta", exist_ok=True)
        print(f"Splitting {input_fasta} into single fasta files:")
        for entry in entries:
            output_file = f"{output_dir}fasta/{entry.id}.fa"
            SeqIO.write(entry, output_file, "fasta")
            print(f"Created: {output_file}")
        fasta_dir = f"{output_dir}fasta"
        return fasta_dir


def split_and_overlap(input_fasta, max_size=5000000, fixed_length=194):
    # Split a given seq in smaller overlapping chunks (in individual fastas)
    # only if len(seq) is > max_size
    input_dir = input_fasta.rsplit('/', maxsplit=1)[0]
    seq_ = SeqIO.read(input_fasta, "fasta").seq
    chr_id = SeqIO.read(input_fasta, "fasta").id

    # Get the number of fasta files (chunks) that will be created
    last_chunk = 1
    if len(seq_) % max_size == 0:  # no last chunk (max size is a 
        last_chunk = 0             # factor of len(seq_))
    num_files = (len(seq_) // max_size) + last_chunk
    
    # Each split must overlap by a window size - 1 nt to not miss any window
    # when splitting
    overlap = fixed_length - 1

    # For seq of length between max_size and 2*max_size nt long, split in 
    # two ~ equal chunks
    if max_size < len(seq_) <= max_size * 2:
        print(f"\nConverting {chr_id}.fa into 2 smaller fasta chunks:")
        # Create first split file
        output_file = f"{input_dir}/{chr_id}_chunk_0.fa"
        with open(output_file, 'w') as f:
            f.write(f">{chr_id}_chunk_0  prev_size=0\n")
            f.write(str(seq_[0:(len(seq_)//2)]))
        # Create 2nd split file
        output_file2 = f"{input_dir}/{chr_id}_chunk_1.fa"
        # Nb of nt before that chunk (total of previous chunks)
        prev_size = len(str(seq_[0:(len(seq_)//2)]))
        with open(output_file2, 'w') as f2:
            f2.write(f">{chr_id}_chunk_1  prev_size={prev_size}\n")
            f2.write(str(seq_[(len(seq_)//2)-overlap:]))
        print(f"Created: {output_file} {output_file2}")

        # Remove initial input fasta as it has been converted to smaller chunks
        sp.call(f"rm {input_fasta}", shell=True)


    # For seq > 2*max_size nt
    elif len(seq_) > max_size * 2:
        print(f"\nConverting {chr_id}.fa into {num_files} smaller fasta chunks:")
        print(len(seq_))
        #print(chr_id, len(seq_), num_files)
        prev_size = 0
        for i in range(num_files):
            start = max(0, i * max_size - overlap)
            end = min((i + 1) * max_size, len(seq_))
            print(start, end, len(seq_[start:end]), prev_size)
            # Create a new file for each split
            output_file3 = f"{input_dir}/{chr_id}_chunk_{i}.fa"
            with open(output_file3, 'w') as f3:
                f3.write(f">{chr_id}_chunk_{i}  prev_size={prev_size}\n")
                f3.write(str(seq_[start:end]))
                print(f"Created: {output_file3}")
            prev_size += len(seq_[start:end])

        # Remove initial input fasta as it has been converted to smaller chunks
        sp.call(f"rm {input_fasta}", shell=True)
    

if __name__ == "__main__":
    input_ = "data/references/genome_fa/arabidopsis_thaliana_genome.fa"
    #input_ = "data/references/genome_fa/schizosaccharomyces_pombe_genome.fa"
    #input_ = "data/references/genome_fa/candida_albicans/Ca22chrM_C_albicans_SC5314.fa"
    #output_prefix = "2R"
    #num_files = 10
    out_ = 'out_test_sp/'
    fasta_dir = multi_fasta(input_, out_)
    if fasta_dir is not None:
        records = [i for i in os.listdir(fasta_dir) if '_chunk_' not in i]
        for fa in records:
            #print(len(SeqIO.read(f"{fasta_dir}/{fa}", "fasta").seq))
            split_and_overlap(f"{fasta_dir}/{fa}", max_size=10000000)
    #else:
        #split the single fasta file
    #split_and_overlap(input_fasta, output_prefix, num_files)

