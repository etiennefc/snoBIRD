#!/usr/bin/python3
import numpy
from Bio import SeqIO
import os
import subprocess as sp

# Load variables
chunks = snakemake.params.chunks
chunk_size = snakemake.params.chunk_size
input_ = snakemake.input.input_fasta
out_ = snakemake.output.split_chr_dir


def multi_fasta(input_fasta, output_dir):
    # If multiple sequences in the input fasta file, split in separate 
    # fasta files per chr
    entries = list(SeqIO.parse(input_fasta, "fasta"))
    if len(entries) > 1:
        os.makedirs(f"{output_dir}/fasta", exist_ok=True)
        print(f"Splitting {input_fasta} into single fasta files:")
        for entry in entries:
            output_file = f"{output_dir}/fasta/{entry.id}.fa"
            SeqIO.write(entry, output_file, "fasta")
            print(f"Created: workflow/{output_file}")
        fasta_dir = f"{output_dir}/fasta"
        return fasta_dir


def split_and_overlap(input_fasta, max_size=5000000, fixed_length=194):
    # Split a given seq in smaller overlapping chunks (in individual fastas)
    # only if len(seq) is > max_size
    input_dir = input_fasta.rsplit('/', maxsplit=1)[0]
    seq_ = SeqIO.read(input_fasta, "fasta").seq
    chr_id = SeqIO.read(input_fasta, "fasta").id

    # Get the number of fasta files (chunks) that will be created
    last_chunk = 0
    # If last chunk is >= fixed_length, create a real last chunk; otherwise, 
    # the last x nt (where x<fixed_length) are concatenated to the previous 
    # chunk to form the real last chunk. This is only for seq_ > 2 * max_size
    if (len(seq_) % max_size != 0) & (len(seq_) % max_size >= fixed_length) & (
        len(seq_) > max_size * 2):
        last_chunk = 1
    if len(seq_) < max_size:
        last_chunk = 1  # last and only chunk
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
            f.write(f">{chr_id}_chunk_0  prev_size=0 total_size={len(seq_)}\n")
            f.write(str(seq_[0:(len(seq_)//2)]))
        # Create 2nd split file
        output_file2 = f"{input_dir}/{chr_id}_chunk_1.fa"
        # Nb of nt before that chunk (total of previous chunks)
        prev_size = len(str(seq_[0:(len(seq_)//2)]))
        with open(output_file2, 'w') as f2:
            f2.write(f">{chr_id}_chunk_1  prev_size={prev_size} "+
                    f"total_size={len(seq_)}\n")
            f2.write(str(seq_[(len(seq_)//2)-overlap:]))
        print(f"Created: workflow/{output_file} workflow/{output_file2}")

        # Remove initial input fasta as it has been converted to smaller chunks
        sp.call(f"rm {input_fasta}", shell=True)


    # For seq > 2*max_size nt
    elif len(seq_) > max_size * 2:
        print(f"\nConverting {chr_id}.fa into {num_files} smaller fasta chunks:")
        for i in range(num_files):
            start = max(0, i * max_size - overlap)
            end = min((i + 1) * max_size, len(seq_))
            if i == num_files - 1:  # for remainders of seq < fixed_length, 
                end = len(seq_) # concat to previous sequence
            prev_size = max(0, len(seq_[0:start]) - 1) 
            # Create a new file for each split
            output_file3 = f"{input_dir}/{chr_id}_chunk_{i}.fa"
            with open(output_file3, 'w') as f3:
                f3.write(f">{chr_id}_chunk_{i}  prev_size={prev_size} "+
                        f"total_size={len(seq_)}\n")
                f3.write(str(seq_[start:end]))
                print(f"Created: workflow/{output_file3}")
            prev_size += len(seq_[start:end])

        # Remove initial input fasta as it has been converted to smaller chunks
        sp.call(f"rm {input_fasta}", shell=True)
    

if __name__ == "__main__":

    # Define chunk formation based on the chunks flag
    fasta_dir = multi_fasta(input_, out_)
    if chunks == None:  # no chunk divisions
        if fasta_dir == None:   # only 1 chr in that initial fasta, so copy 
                                # it at the right location
            sp.call(f"mkdir -p {out_}/fasta && cp {input_} {out_}/fasta/", 
                    shell=True)
            print(f'Copying your fasta file to workflow/{out_}/fasta/')
                # if multiple chr in the initial fasta and no chunk formation, 
                # multi_fasta has already split it into separate fasta files
    elif chunks == True:  # division into chunks
        if fasta_dir is not None:  # multiple chr in the initial fasta
            records = [i for i in os.listdir(fasta_dir) if '_chunk_' not in i]
            for fa in records:
                split_and_overlap(f"{fasta_dir}/{fa}", max_size=chunk_size)
        else:  # only 1 chr in initial fasta, but might be splittable in chunks
            sp.call(f"mkdir -p {out_}/fasta && cp {input_} {out_}/fasta/", 
                    shell=True)
            fa_ = os.listdir(f"{out_}/fasta/")[0]
            split_and_overlap(f"{out_}/fasta/{fa_}", max_size=chunk_size)


