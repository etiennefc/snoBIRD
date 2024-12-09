""" Define useful functions that are used across rules."""
import warnings
warnings.formatwarning = lambda msg, *args: f"\033[93m{msg}\033[0m\n"

def fasta_iter(input_fasta):
    ''' Return a generator of name, seq in a given fasta. 
        As fast as SeqIO.parse.'''
    grep_cmd = f"""grep -zoP ">[^\\n]*\\n>" {input_fasta}"""
    missing_seq = sp.run(grep_cmd, 
                    capture_output=True, text=True, shell=True).stdout
    if '>' in missing_seq:
        chr_name = missing_seq[0:-2].strip()
        raise ValueError(f'The entry "{chr_name}" in {input_fasta} has '+
                        'no sequence line below. Either remove that entry or '+
                        f'add its sequence under the "{chr_name}" line.')
    
    with open(input_fasta, 'r') as fasta_:
        # Just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in groupby(fasta_, lambda line: line[0] == ">"))

        for header in faiter:
            # drop the ">"
            headerStr = next(header)[1:].split()[0]
            if headerStr == '>':
                white_space_entry = next(header)[1:]
                raise ValueError(f'The entries in input fasta {input_fasta}'+
                    ' must be as follows: ">chr12", ">12" or ">X", thereby '+
                    'having NO space, tab or newline between the ">" and '+
                    f'entry name. Please fix the entry "{white_space_entry}"'+
                    ' (and maybe other subsequent sequences) to remove the '+
                    'whitespace between ">" and entry_name.')

            # join all sequence lines to one.
            seq_ = "".join(s.strip() for s in next(faiter))
            #print(headerStr, seq_)
            if len(seq_) == 0:
                raise ValueError(f'The entry >"{headerStr}" in {input_fasta}'+
                        ' has an empty sequence line below. Either remove '+
                        'that entry or add its sequence under the '+
                        f'">{headerStr}" line.')

            yield (headerStr, len(seq_))



def validate_fasta(input_fasta, valid_ext=('.fasta', '.fa')):
    ''' Validate fasta file name, name of entries and file format.'''
    # Check file extension
    if not input_fasta.endswith(valid_ext):
        raise ValueError("Invalid file extension for input fasta. "+
                f"Expected one of {valid_ext}.")



def is_sbatch_installed2():
    """ Find if SLURM's sbatch command is installed to see if the user is on a 
        local computer or a SLURM cluster."""
    try:
        # Try to run sbatch --version to check if sbatch is installed
        result = sp.run(['sbatch', '--version'], stdout=sp.PIPE, 
                    stderr=sp.PIPE, text=True)
        
        # If the command was successful, return True
        if result.returncode == 0:
            return True
        else:
            return False
    except FileNotFoundError:
        # If sbatch is not found, return False
        return False



def get_chunk_names(len_seq, chr_name, max_size=5000000, fixed_length=194):
    ''' Get the number of fasta chunks that will be created with --chunks 
        (including the untouched fastas); return a list of chr/chunk names.'''
    # Get the number of fasta files (chunks) that will be created
    last_chunk = 0
    # If last chunk is >= fixed_length, create a real last chunk; otherwise, 
    # the last x nt (where x<fixed_length) are concatenated to the previous 
    # chunk to form the real last chunk. This is only for seq_ > 2 * max_size
    if (len_seq % max_size != 0) & (len_seq % max_size >= fixed_length) & (
        len_seq > max_size):
        last_chunk = 1
    elif len_seq < max_size:
        last_chunk = 1  # last and only chunk
    num_files = (len_seq // max_size) + last_chunk
    
    # Return chunk names
    if num_files == 1:
        return [chr_name], [len_seq]
    if num_files == 2:
        return [f'{chr_name}_chunk_{i}' for i in range(num_files)], [
            int(len_seq/2)] * 2
    if num_files > 2:
        return [f'{chr_name}_chunk_{i}' for i in range(num_files)], [
            max_size] * num_files



def get_chr_names(input_fasta, chunk_value, chunk_size, dryrun=None, 
                gpu="A100"):
    ''' Main function to get the number of fasta files (one for each chr and/or
        chunks of chr) that will be created by the rule split_chr from an 
        input fasta. This serves the purpose of creating the Snakemake DAG and 
        knowing which chr_ wildcards will be created without the need of 
        checkpoints.'''

    # Check file extension
    validate_fasta(input_fasta)

    # Find all chr_names and sizes
    all_chr = {}
    for chr_, size in fasta_iter(input_fasta):
        if chr_ not in all_chr.keys():
            all_chr[chr_] = size
        else:
            raise ValueError(f"The following duplicate entry names, '>{chr_}'"+
            f", were found in input fasta {input_fasta}. Please change "+
            "sequence names to have only unique entries across the fasta file."
            )

    if len(all_chr.keys()) == 0:
        raise ValueError(f'No entries in input fasta {input_fasta} '+
        'could be found. Entry names should start with a ">" as follows: '+
        '">chr12", ">12" or ">X".')
    elif (len(all_chr.keys()) > 25) & (dryrun==True):
        num_chr = len(all_chr.keys())
        # Flag potential scaffolds at the end of the input fasta 
        # (if more chr than expected)
        warnings.warn(
        "\nUserWarning: It seems like your input fasta contains '"+
        str(num_chr)+"' entries (marked by '>'). If this is the "+
        "number of chromosomes/sequences you expect, please ignore this "+
        "message. Otherwise, you should remove unwanted sequences in your "+
        "input (usually scaffolds at the end of genome fasta files) as it "+
        "will increase SnoBIRD runtime and use GPUs unnecessarily on sequences"+
        " that are not of interest. To see which sequence/chromosome entries "+
        "are present in your input fasta, run the following command:"+
        "\n\tgrep '>' <input_fasta.fa>\n")
    
    # Check GPU type that user provided vs what is present on the HPC
    if is_sbatch_installed2() == True:  # verify if on HPC cluster
        if gpu != 'Unknown':
            gen = '['+gpu[0].upper()+gpu[0].lower()+']'
            num_ = gpu[1:]
            result = sp.run(['sinfo -o "%N %G" | grep "'+gen+num_+'"'], shell=True, 
                    capture_output=True, text=True).stdout.strip()
            if result == '':
                warnings.warn(
                '\nUserWarning: It seems that the GPU generation that you '+
                'provided "'+gpu+'" (using -G/--gpu_generation, '+
                'default: "A100") is not present in the available GPUs on '+
                'your cluster. If you are certain that it is the right GPU '+
                'generation and that it is available, please ignore this '+
                'message. Otherwise, please choose carefully your GPU '+
                'generation as it will strongly impact SnoBIRD runtime (see '+
                'README.md for more info). You can view which GPU generation '+
                'is available using the following command:'+
                '\n\tsinfo -o "%N %G" | grep "[pPaAvVhH]100"\n')   

    

    if chunk_value == None:  # no chunks will be created
        return all_chr  # return dict of chr and sizes
    if chunk_value == True:  # chunks will be created for big chromosomes
        all_chr_chunks, all_chr_chunks_sizes = [], []
        if len(all_chr.keys()) > 1:  # multiple chr in the initial fasta
            for c, chr_size_ in all_chr.items():
                chunk_names, c_size = get_chunk_names(chr_size_, c, 
                            max_size=chunk_size)
                all_chr_chunks.extend(chunk_names)
                all_chr_chunks_sizes.extend(c_size)
        else:  # only 1 chr in initial fasta, but might be splittable in chunks
            only_chr = list(all_chr.keys())[0]
            only_chr_size = list(all_chr.values())[0]
            chunk_names, c_size = get_chunk_names(only_chr_size, only_chr, 
                        max_size=chunk_size)
            all_chr_chunks.extend(chunk_names)
            all_chr_chunks_sizes.extend(c_size)
        return all_chr_chunks, all_chr_chunks_sizes










