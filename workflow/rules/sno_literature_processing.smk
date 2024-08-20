rule get_chr_size_literature:
    """ Get the chr size from genome fasta."""
    input:
        genome = get_species_genome
    output:
        chr_size = 'data/references/chr_size/{sno_fasta}_chr_size.tsv'
    wildcard_constraints:
        sno_fasta = join_list(config['sno_fasta'], config['sno_fasta'])
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input.genome} --fai-idx {wildcards.sno_fasta}_temp.fai && "
        "cut -f1,2 {wildcards.sno_fasta}_temp.fai > {output.chr_size} && "
        "rm {wildcards.sno_fasta}_temp.fai"

rule blat_sno_genome:
    """ Get the genomic location of a given snoRNA sequence 
        in a given species genome using BLAT (expressed C/D box 
        snoRNAs collected from the literature)."""
    input:
        blat_fake_dependency = rules.download_blat.output.tmp_file,
        sno_sequences = "data/sno_literature_processing/{sno_fasta}.fa",
        genome = get_species_genome
    output:
        sno_location = "data/sno_literature_processing/{sno_fasta}_location.psl"
    params:
        blat_path = "data/references/blat/blat"
    shell:
        "{params.blat_path} {input.genome} {input.sno_sequences} "
        "{output.sno_location}"

rule format_blat_output:
    """ Format BLAT output to keep only the match with the highest 
        number of matching nucleotide according to the original 
        snoRNA sequence. Update the snoRNA sequence based on the 
        actual location in the species genome (for merged snoRNAs 
        and potential mismatches mainly)."""
    input:
        blat = rules.blat_sno_genome.output.sno_location,
        genome = get_species_genome,
        chr_size = rules.get_chr_size_literature.output.chr_size
    output:
        df = 'data/sno_literature_processing/{sno_fasta}_location_formated.tsv'
    params:
        dataset_attribute = lambda wildcards: config['dataset_attributes'][wildcards.sno_fasta],
        extension = 15
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/format_blat_output.py"       

rule merge_sno_location_species:
    """ Merge all snoRNA info (from format_blat_output) from 
        all species in a single table."""
    input:
        dfs = expand(rules.format_blat_output.output.df, **config)
    output:
        df = 'data/references/sno_location_seq_all_species.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/merge_sno_location_species.py"