""" Define useful functions that are used across rules."""

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
    if not input_fasta.lower().endswith(valid_ext):
        raise ValueError("Invalid file extension for input fasta. "+
                f"Expected one of {valid_ext}.")
    if input_fasta.lower().endswith('.fasta'):  # change .fasta for .fa
        sp.call(f'mv {input_fasta} {input_fasta[:-6]}.fa', shell=True)



def get_chunk_names(len_seq, chr_name, max_size=5000000):
    ''' Get the number of fasta chunks that will be created with --chunks 
        (including the untouched fastas); return a list of chr/chunk names.'''
    # Get the number of fasta files (chunks) that will be created
    last_chunk = 1
    if len_seq % max_size == 0:  # no last chunk (max size is a 
        last_chunk = 0             # factor of len(seq_))
    num_files = (len_seq // max_size) + last_chunk
    # Return chunk names
    if num_files == 1:
        return [chr_name]
    if num_files > 1:
        return [f'{chr_name}_chunk_{i}' for i in range(num_files)]



def get_chr_names(input_fasta, chunk_value, chunk_size):
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

    
    if chunk_value == None:  # no chunks will be created
        return list(all_chr.keys())
    
    if chunk_value == True:  # chunks will be created for big chromosomes
        all_chr_chunks = []
        if len(all_chr.keys()) > 1:  # multiple chr in the initial fasta
            for c, chr_size_ in all_chr.items():
                chunk_names = get_chunk_names(chr_size_, c, 
                            max_size=chunk_size)
                all_chr_chunks.extend(chunk_names)
        else:  # only 1 chr in initial fasta, but might be splittable in chunks
            only_chr = all_chr.keys()[0]
            only_chr_size = all_chr.values()[0]
            chunk_names = get_chunk_names(only_chr_size, only_chr, 
                        max_size=chunk_size)
            all_chr_chunks.extend(chunk_names)
        return all_chr_chunks











