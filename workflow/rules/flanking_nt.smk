rule c_d_box_location:
    """ Find C, D, C' and D' motif and location (if they exist) in C/D snoRNAs."""
    input:
        cd_fa = 'data/references/all_expressed_cd_sequences.fa'
    output:
        c_d_box_location = 'data/references/all_expressed_c_d_and_prime_location.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/cd_box_location.py"

rule cd_sno_location_seq:
    """ Get sno location using bed of all RNAs downloaded in RNAcentral. 
        Also return sno sequence and 15 nt sequence upstream of C box"""
    input:
        box_location = rules.c_d_box_location.output.c_d_box_location,
        sno_sequences = 'data/references/all_expressed_cd_sequences.fa',
        sno_df = 'data/references/all_expressed_cd_sequences_location.tsv',
        chr_size = [sp for sp in expand(rules.get_chr_size_literature.output.chr_size, sno_fasta=config['sno_fasta']) + 
                    expand(rules.get_chr_size_tgirt.output.chr_size, species=["homo_sapiens", "mus_musculus", 
                    "saccharomyces_cerevisiae"]) if "schizosaccharomyces_pombe" not in sp],
        #chr_size = [get_chr_size(sp) for sp in [s for s in config['species'] 
         #               + config['species_tgirt'] if s != "schizosaccharomyces_pombe"]],
        genome = get_all_genomes('data/references/genome_fa/*.fa')
    output:
        fifteen_upstream_c_fa = "data/references/fifteen_nt_upstream_C_box.fa",
        fifteen_downstream_d_fa = "data/references/fifteen_nt_downstream_D_box.fa"
    params:
        species_name = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/cd_sno_location.py"

rule logomaker_flanking_nt_c_d:
    """ Create logo of the 15 nt upstream of the C box of C/D snoRNAs and 
        15 nt downstream of D boxes."""
    input:
        c_box_fasta = rules.cd_sno_location_seq.output.fifteen_upstream_c_fa,
        d_box_fasta = rules.cd_sno_location_seq.output.fifteen_downstream_d_fa
    output:
        c_logo = 'results/figures/logo/fifteen_nt_upstream_of_C.svg',
        d_logo = 'results/figures/logo/fifteen_nt_downstream_of_D.svg'
    conda:
        "../envs/logomaker.yaml"
    script:
        "../scripts/python/figures/logo_up_downstream_c_d_box.py"

