# For the negatives, also include snoRNA pseudogenes (control for rfam_id?)

rule random_shuffle_sno_fixed_length:
    """ Shuffle expressed C/D sequences and use that as 
        negative examples from all species"""
    input:
        sno_sequences = rules.get_sno_sequences_fixed_length.output.df
    output:
        shuffled_sno_df = 'data/references/negatives/shuffle_sno/random_shuffle_all_expressed_cd_fixed_length_{fixed_length}nt.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/random_shuffle_sno_fixed_length.py"

rule filter_rnacentral_tRNA_snRNA_pre_miRNA_fixed_length:
    """ Filter rnacentral ncRNA beds (keep non-overlapping entries 
        of tRNA, snRNA and pre-miRNAs) and return also sequence as 1 df.
        O. tauri is not on RNAcentral (thereby excluded here), S. pombe is on 
        RNAcentral but not included in our training set (thereby excluded here) and 
        T. thermophila is on RNAcentral (but relative to the previous incomplete/unassembled 
        genome version, thereby not usable/translatable to the latest assembled version of 
        its genome with which we work here, thereby excluded here). Extend their sequence 
        to obtain x nt (where x is a fixed length as with the C/D snoRNAs, 
        e.g., 211 nt (181 nt + 15 nt + 15 nt), where 181 = len(snoRNA below the 95th 
        percentile))."""
    input:
        beds = expand(rules.download_rnacentral_ncRNA.output.bed, 
                species=[sp for sp in config['species'] + config['species_tgirt'] 
                            if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe', 
                            'tetrahymena_thermophila', 'candida_albicans']]),
        genomes = get_all_genomes('data/references/genome_fa/*.fa'),
        chr_sizes = get_all_genomes('data/references/chr_size/*.tsv'),
        gallus_gallus_gff = rules.download_gallus_gallus_gff.output.gff
    output:
        df = 'data/references/rnacentral/filtered_tRNA_snRNA_pre_miRNA_fixed_length_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_rnacentral_tRNA_snRNA_pre_miRNA_fixed_length.py"

rule format_blat_haca_output_fixed_length:
    """ Format BLAT output to keep only the match with the highest 
        number of matching nucleotide according to the original 
        snoRNA sequence. Update the snoRNA sequence based on the 
        actual location in the species genome (for merged snoRNAs 
        and potential mismatches mainly) and extend to fixed_length nt of length."""
    input:
        blat = rules.blat_haca_genome.output.sno_location,
        genome = get_species_genome,
        chr_size = get_chr_size
    output:
        df = 'data/references/HACA/{species}_HACA_location_formated_fixed_length_{fixed_length}nt.tsv'
    params:
        dataset_attribute = config['dataset_attributes']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/format_blat_haca_output_fixed_length.py"

rule select_random_intergenic_intronic_exonic_regions_fixed_length:
    """ Select random intronic, exonic and intergenic regions in the genomes 
        of various species. Make sure that these random regions 
        have a length of fixed_length nt"""
    input:
        expressed_cd_all_sets = rules.tuning_train_test_split_rfam_fixed_length.output.all_positives,
        intronic_regions = rules.get_intergenic_intronic_exonic_regions.output.intronic_regions_bed,
        intergenic_regions = rules.get_intergenic_intronic_exonic_regions.output.intergenic_regions_bed,
        exonic_regions = rules.get_intergenic_intronic_exonic_regions.output.exonic_regions_bed,
        genome = get_species_genome
    output:
        random_intronic_regions = 'data/references/negatives/random_regions/selected_intronic_regions_{species}_fixed_length_{fixed_length}nt.bed',
        random_intergenic_regions = 'data/references/negatives/random_regions/selected_intergenic_regions_{species}_fixed_length_{fixed_length}nt.bed',
        random_exonic_regions = 'data/references/negatives/random_regions/selected_exonic_regions_{species}_fixed_length_{fixed_length}nt.bed'
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/select_random_intergenic_intronic_exonic_regions_fixed_length.py"



rule get_all_initial_negatives_fixed_length:
    """ From all negative examples (other ncRNA sequences (H/ACA, 
        tRNA, snRNA, pre-miRNA), shuffle of C/D sequences, random 
        sequences in introns, exons, and intergenic regions, and potentially
        snoRNA pseudogene sequences?), select the wanted proportion 
        of each of these negatives relative to the number of positive 
        examples (expressed C/D). The positives are used to filter out 
        shuffled CD sequences that are not present in the positives 
        (because we limited the number of sno per Rfam family to be max 10)"""
    input:
        ncRNA = rules.filter_rnacentral_tRNA_snRNA_pre_miRNA_fixed_length.output.df,
        haca = expand(rules.format_blat_haca_output_fixed_length.output.df, 
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe', 
                        'tetrahymena_thermophila']], fixed_length=config['fixed_length']),
        shuffle_sno = rules.random_shuffle_sno_fixed_length.output.shuffled_sno_df,
        random_intronic_regions = expand(rules.select_random_intergenic_intronic_exonic_regions_fixed_length.output.random_intronic_regions,
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe']], fixed_length=config['fixed_length']),
        random_intergenic_regions = expand(rules.select_random_intergenic_intronic_exonic_regions_fixed_length.output.random_intergenic_regions,
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe']], fixed_length=config['fixed_length']),
        random_exonic_regions = expand(rules.select_random_intergenic_intronic_exonic_regions_fixed_length.output.random_exonic_regions,
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe']], fixed_length=config['fixed_length']),
        positives = rules.tuning_train_test_split_rfam_fixed_length.output.all_positives,
        human_snoRNA_pseudogenes = rules.get_human_snoRNA_pseudogenes.output.pseudogenes,
        mouse_snoRNA_pseudogenes = rules.get_mouse_snoRNA_pseudogenes.output.pseudogenes,
        droso_snoRNA_pseudogenes = rules.get_drosophila_snoRNA_pseudogenes.output.pseudogenes,
        rfam_pseudo = rules.filter_infernal_pseudosno.output.df,
        rfam_clans = rules.download_rfam_clans.output.df
    output:
        tuning = 'data/references/negatives/initial/negatives_tuning_set_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/negatives/initial/negatives_training_set_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/negatives/initial/negatives_test_set_fixed_length_{fixed_length}nt.tsv',
        all_negatives = 'data/references/negatives/initial/all_negatives_fixed_length_{fixed_length}nt.tsv'
    params:
        short_name_dict = config['species_short_name'],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_all_initial_negatives_fixed_length.py"

rule get_all_initial_negatives_wo_pseudo_fixed_length:
    """ From all negative examples (other ncRNA sequences (H/ACA, 
        tRNA, snRNA, pre-miRNA), shuffle of C/D sequences, random 
        sequences in introns, exons, and intergenic regions, select the wanted proportion 
        of each of these negatives relative to the number of positive 
        examples (expressed C/D). The positives are used to filter out 
        shuffled CD sequences that are not present in the positives 
        (because we limited the number of sno per Rfam family to be max 100)"""
    input:
        ncRNA = rules.filter_rnacentral_tRNA_snRNA_pre_miRNA_fixed_length.output.df,
        haca = expand(rules.format_blat_haca_output_fixed_length.output.df, 
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe', 
                        'tetrahymena_thermophila']], fixed_length=config['fixed_length']),
        shuffle_sno = rules.random_shuffle_sno_fixed_length.output.shuffled_sno_df,
        random_intronic_regions = expand(rules.select_random_intergenic_intronic_exonic_regions_fixed_length.output.random_intronic_regions,
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe']], fixed_length=config['fixed_length']),
        random_intergenic_regions = expand(rules.select_random_intergenic_intronic_exonic_regions_fixed_length.output.random_intergenic_regions,
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe']], fixed_length=config['fixed_length']),
        random_exonic_regions = expand(rules.select_random_intergenic_intronic_exonic_regions_fixed_length.output.random_exonic_regions,
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe']], fixed_length=config['fixed_length']),
        positives = rules.tuning_train_test_split_rfam_sno_pseudo_fixed_length.output.all_positives
    output:
        tuning = 'data/references/negatives/initial/negatives_tuning_set_wo_pseudo_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/negatives/initial/negatives_training_wo_pseudo_set_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/negatives/initial/negatives_test_set_wo_pseudo_fixed_length_{fixed_length}nt.tsv',
        all_negatives = 'data/references/negatives/initial/all_negatives_wo_pseudo_fixed_length_{fixed_length}nt.tsv'
    params:
        short_name_dict = config['species_short_name'],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_all_initial_negatives_wo_pseudo_fixed_length.py"

rule get_all_initial_negatives_wo_pseudo_fixed_length_data_aug:
    """ From all negative examples (other ncRNA sequences (H/ACA, 
        tRNA, snRNA, pre-miRNA), shuffle of C/D sequences, random 
        sequences in introns, exons, and intergenic regions, select the wanted proportion 
        of each of these negatives relative to the number of positive 
        examples (expressed C/D). The positives are used to filter out 
        shuffled CD sequences that are not present in the positives 
        (because we limited the number of sno per Rfam family to be max 100)"""
    input:
        ncRNA = rules.filter_rnacentral_tRNA_snRNA_pre_miRNA_fixed_length.output.df,
        haca = expand(rules.format_blat_haca_output_fixed_length.output.df, 
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe', 
                        'tetrahymena_thermophila']], fixed_length=config['fixed_length']),
        shuffle_sno = rules.random_shuffle_sno_fixed_length.output.shuffled_sno_df,
        random_intronic_regions = expand(rules.select_random_intergenic_intronic_exonic_regions_fixed_length.output.random_intronic_regions,
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe']], fixed_length=config['fixed_length']),
        random_intergenic_regions = expand(rules.select_random_intergenic_intronic_exonic_regions_fixed_length.output.random_intergenic_regions,
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe']], fixed_length=config['fixed_length']),
        random_exonic_regions = expand(rules.select_random_intergenic_intronic_exonic_regions_fixed_length.output.random_exonic_regions,
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe']], fixed_length=config['fixed_length']),
        positives = rules.tuning_train_test_split_rfam_sno_pseudo_fixed_length_data_aug.output.all_positives
    output:
        tuning = 'data/references/negatives/data_augmentation/negatives_tuning_set_wo_pseudo_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/negatives/data_augmentation/negatives_training_wo_pseudo_set_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/negatives/data_augmentation/negatives_test_set_wo_pseudo_fixed_length_{fixed_length}nt.tsv',
        all_negatives = 'data/references/negatives/data_augmentation/all_negatives_wo_pseudo_fixed_length_{fixed_length}nt.tsv'
    params:
        short_name_dict = config['species_short_name'],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_all_initial_negatives_wo_pseudo_fixed_length_data_aug.py"    
