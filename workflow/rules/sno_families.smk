rule get_sno_sequences:
    """ Get all expressed C/D snoRNA sequences (from the literature and 
        TGIRT-seq) regrouped in a fasta file."""
    input:
        sno_literature = rules.merge_sno_location_species.output.df,
        sno_tgirt = expand(rules.get_expressed_snoRNAs_location.output.expressed_sno_df, 
                    species=['homo_sapiens', 'mus_musculus', 'saccharomyces_cerevisiae']),
        sno_tgirt_droso = expand(rules.get_D_melanogaster_expressed_snoRNAs_location.output.expressed_sno_df, 
                        species='drosophila_melanogaster')
    output:
        fa = 'data/references/all_expressed_cd_sequences.fa',
        df = 'data/references/all_expressed_cd_sequences_location.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_sno_sequences.py"

rule get_sno_sequences_fixed_length:
    """ Get all expressed C/D snoRNA sequences (from the literature and 
        TGIRT-seq) regrouped in a fasta file. Get a fixed extended 
        length (i.e. longest expressed snoRNA (below the 95th 
        percentile) + 15 nt up/downstream)"""
    input:
        sno_literature = rules.merge_sno_location_species.output.df,
        sno_tgirt = expand(rules.get_expressed_snoRNAs_location.output.expressed_sno_df, 
                    species=['homo_sapiens', 'mus_musculus', 'saccharomyces_cerevisiae']),
        sno_tgirt_droso = expand(rules.get_D_melanogaster_expressed_snoRNAs_location.output.expressed_sno_df, 
                        species='drosophila_melanogaster'),
        genomes = get_all_genomes('data/references/genome_fa/*.fa'),
        chr_size = get_all_genomes('data/references/chr_size/*.tsv')
    output:
        fa = 'data/references/all_expressed_cd_sequences_fixed_length_{fixed_length}nt.fa',
        df = 'data/references/all_expressed_cd_sequences_location_fixed_length_{fixed_length}nt.tsv'
    params:
        species_short_name = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_sno_sequences_fixed_length.py"

rule infernal:
    """ Use Infernal and Rfam covariance 
        models to identify the Rfam family 
        of all expressed C/D box snoRNAs."""
    input:
        sno_fasta = rules.get_sno_sequences.output.fa,
        rfam_cm = rules.download_rfam_covariance_models.output.rfam_cm
    output:
        infernal_tblout = 'data/references/infernal/sno_families.tblout',
        infernal_alignments = 'data/references/infernal/sno_families.txt'
    conda:
        "../envs/infernal.yaml"
    shell:
        "cmpress -F {input.rfam_cm} && "
        "cmscan --cut_ga --rfam --nohmmonly -o {output.infernal_alignments} "
        "--tblout {output.infernal_tblout} {input.rfam_cm} {input.sno_fasta}"

rule filter_infernal:
    """ Filter infernal output to return the Rfam family id per snoRNA."""
    input:
        infernal = rules.infernal.output.infernal_tblout,
        sno_df = rules.get_sno_sequences.output.df
    output:
        df = 'data/references/infernal/sno_families_filtered.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_infernal.py"


#rule tuning_train_test_split_rfam:
#    """ Split expressed C/D in 3 datasets (tuning (10%), training (70%) 
#        and test set (20%)). SnoRNAs of a same Rfam clan (then Rfam family) 
#        are all kept within the same set so that we limit overfitting. """
#    input:
#        sno_rfam = rules.filter_infernal.output.df,
#        sno_literature = rules.merge_sno_location_species.output.df,
#        sno_tgirt = expand(rules.get_expressed_snoRNAs_location.output.expressed_sno_df, 
#                                species=['homo_sapiens', 'mus_musculus', 'saccharomyces_cerevisiae']),
#        rfam_clans = rules.download_rfam_clans.output.df
#    output:
#        tuning = 'data/references/infernal/cd_rfam_filtered_tuning_set.tsv',
#        training = 'data/references/infernal/cd_rfam_filtered_training_set.tsv',
#        test = 'data/references/infernal/cd_rfam_filtered_test_set.tsv',
#        all_positives = 'data/references/positives/cd_rfam_filtered_all.tsv'
#    params:
#        random_seed = 42
#    conda:
#        "../envs/python_new.yaml"
#    script:
#        "../scripts/python/tuning_train_test_split_rfam.py"

rule tuning_train_test_split_rfam_fixed_length:
    """ Split expressed C/D in 3 datasets (tuning (10%), training (70%) 
        and test set (20%)). SnoRNAs of a same Rfam clan (then Rfam family) 
        are all kept within the same set so that we limit overfitting. """
    input:
        sno_rfam = rules.filter_infernal.output.df,
        sno_literature = rules.merge_sno_location_species.output.df,
        sno_tgirt = expand(rules.get_expressed_snoRNAs_location.output.expressed_sno_df, 
                                species=['homo_sapiens', 'mus_musculus', 'saccharomyces_cerevisiae']),
        sno_tgirt_droso = expand(rules.get_D_melanogaster_expressed_snoRNAs_location.output.expressed_sno_df, 
                        species='drosophila_melanogaster'),
        rfam_clans = rules.download_rfam_clans.output.df,
        extended_sno_seq = rules.get_sno_sequences_fixed_length.output.df
    output:
        tuning = 'data/references/infernal/cd_rfam_filtered_tuning_set_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/infernal/cd_rfam_filtered_training_set_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/infernal/cd_rfam_filtered_test_set_fixed_length_{fixed_length}nt.tsv',
        all_positives = 'data/references/positives/cd_rfam_filtered_all_fixed_length_{fixed_length}nt.tsv'
    params:
        random_seed = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/tuning_train_test_split_rfam_fixed_length.py"

rule get_human_snoRNA_pseudogenes:
    """ Get human snoRNA pseudogenes (i.e. not expressed) extended sequence up to fixed_length"""
    input:
        human_snoRNA_pseudogenes_dependency = expand(rules.get_expressed_snoRNAs_location.output.expressed_sno_df, species='homo_sapiens'),
        human_genome = expand(rules.download_mammal_genome.output.genome, species='homo_sapiens'),
        human_chr_size = expand(rules.get_chr_size_tgirt.output.chr_size, species='homo_sapiens')
    output:
        pseudogenes = 'data/references/negatives/snoRNA_pseudogenes/homo_sapiens_pseudogene_snoRNAs_{fixed_length}nt.tsv'
    params:
        human_snoRNA_pseudogenes = rules.get_expressed_snoRNAs_location.params.human_pseudosno
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_human_snoRNA_pseudogenes.py"

rule get_mouse_snoRNA_pseudogenes:
    """ Get mouse snoRNA pseudogenes (i.e. not expressed) extended sequence up to fixed_length"""
    input:
        mouse_snoRNA_pseudogenes_dependency = expand(rules.get_expressed_snoRNAs_location.output.expressed_sno_df, species='mus_musculus'),
        mouse_genome = expand(rules.download_mammal_genome.output.genome, species='mus_musculus'),
        mouse_chr_size = expand(rules.get_chr_size_tgirt.output.chr_size, species='mus_musculus')
    output:
        pseudogenes = 'data/references/negatives/snoRNA_pseudogenes/mus_musculus_pseudogene_snoRNAs_{fixed_length}nt.tsv'
    params:
        mouse_snoRNA_pseudogenes = rules.get_expressed_snoRNAs_location.params.mouse_pseudosno
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_mouse_snoRNA_pseudogenes.py"

rule get_drosophila_snoRNA_pseudogenes:
    """ Get drosophila snoRNA pseudogenes (i.e. not expressed) extended sequence up to fixed_length"""
    input:
        droso_snoRNA_pseudogenes = expand(rules.get_D_melanogaster_expressed_snoRNAs_location.output.droso_pseudosno, 
                                            species='drosophila_melanogaster'),
        droso_genome = expand(rules.download_genome.output.genome, species='drosophila_melanogaster'),
        droso_chr_size = expand(rules.get_chr_size_tgirt.output.chr_size, species='drosophila_melanogaster')
    output:
        pseudogenes = 'data/references/negatives/snoRNA_pseudogenes/drosophila_melanogaster_pseudogene_snoRNAs_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_drosophila_snoRNA_pseudogenes.py"

rule fasta_pseudosno:
    """ Create a fasta file of snoRNA pseudogenes (mouse, human and drosphila)."""
    input:
        pseudo_human = expand(rules.get_human_snoRNA_pseudogenes.output.pseudogenes, **config),
        pseudo_mouse = expand(rules.get_mouse_snoRNA_pseudogenes.output.pseudogenes, **config),
        pseudo_droso = expand(rules.get_drosophila_snoRNA_pseudogenes.output.pseudogenes, **config)
    output:
        fa = 'data/references/negatives/snoRNA_pseudogenes/all_pseudosno_actual_sequence.fa'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/fasta_pseudosno.py"

rule infernal_pseudosno:
    """ Use Infernal and Rfam covariance 
        models to identify the Rfam family 
        of all C/D box snoRNA pseudogenes."""
    input:
        sno_fasta = expand(rules.fasta_pseudosno.output.fa, **config),
        rfam_cm = rules.download_rfam_covariance_models.output.rfam_cm
    output:
        infernal_tblout = 'data/references/infernal/pseudosno_families.tblout',
        infernal_alignments = 'data/references/infernal/pseudosno_families.txt'
    conda:
        "../envs/infernal.yaml"
    shell:
        "cmpress -F {input.rfam_cm} && "
        "cmscan --cut_ga --rfam --nohmmonly -o {output.infernal_alignments} "
        "--tblout {output.infernal_tblout} {input.rfam_cm} {input.sno_fasta}"

rule filter_infernal_pseudosno:
    """ Filter infernal output to return the Rfam family id per snoRNA."""
    input:
        infernal = rules.infernal_pseudosno.output.infernal_tblout,
        human_snoRNA_pseudogenes = expand(rules.get_human_snoRNA_pseudogenes.output.pseudogenes, **config),
        mouse_snoRNA_pseudogenes = expand(rules.get_mouse_snoRNA_pseudogenes.output.pseudogenes, **config),
        droso_snoRNA_pseudogenes = expand(rules.get_drosophila_snoRNA_pseudogenes.output.pseudogenes, **config)
    output:
        df = 'data/references/infernal/pseudosno_families_filtered.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_infernal_pseudosno.py"

rule tuning_train_test_split_rfam_sno_pseudo_fixed_length:
    """ Split expressed C/D and pseudogenes in 3 datasets (tuning (10%), training (70%) 
        and test set (20%)). SnoRNAs of a same Rfam clan (then Rfam family) 
        are all kept within the same set so that we limit overfitting. """
    input:
        sno_rfam = rules.filter_infernal.output.df,
        sno_literature = rules.merge_sno_location_species.output.df,
        sno_tgirt = expand(rules.get_expressed_snoRNAs_location.output.expressed_sno_df, 
                                species=['homo_sapiens', 'mus_musculus', 'saccharomyces_cerevisiae']),
        sno_tgirt_droso = expand(rules.get_D_melanogaster_expressed_snoRNAs_location.output.expressed_sno_df, 
                        species='drosophila_melanogaster'),
        rfam_clans = rules.download_rfam_clans.output.df,
        extended_sno_seq = rules.get_sno_sequences_fixed_length.output.df,
        rfam_pseudo = rules.filter_infernal_pseudosno.output.df
    output:
        tuning = 'data/references/infernal/cd_rfam_filtered_tuning_set_sno_pseudo_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/infernal/cd_rfam_filtered_training_set_sno_pseudo_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/infernal/cd_rfam_filtered_test_set_sno_pseudo_fixed_length_{fixed_length}nt.tsv',
        all_positives = 'data/references/positives/cd_rfam_filtered_all_sno_pseudo_fixed_length_{fixed_length}nt.tsv'
    params:
        random_seed = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/tuning_train_test_split_rfam_sno_pseudo_fixed_length.py"

rule tuning_train_test_split_rfam_sno_pseudo_fixed_length_data_aug:
    """ Split expressed C/D and pseudogenes in 3 datasets (tuning (10%), training (70%) 
        and test set (20%)). SnoRNAs of a same Rfam clan (then Rfam family) 
        are all kept within the same set so that we limit overfitting. After split, do 
        data augmentation by taking 5 windows before and 5 after the actual sno window."""
    input:
        tuning = rules.tuning_train_test_split_rfam_sno_pseudo_fixed_length.output.tuning,
        training = rules.tuning_train_test_split_rfam_sno_pseudo_fixed_length.output.training,
        test = rules.tuning_train_test_split_rfam_sno_pseudo_fixed_length.output.test,
        genomes = get_all_genomes('data/references/genome_fa/*.fa')
    output:
        tuning = 'data/references/data_augmentation/cd_rfam_filtered_tuning_set_sno_pseudo_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/data_augmentation/cd_rfam_filtered_training_set_sno_pseudo_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/data_augmentation/cd_rfam_filtered_test_set_sno_pseudo_fixed_length_{fixed_length}nt.tsv',
        all_positives = 'data/references/data_augmentation/cd_rfam_filtered_all_sno_pseudo_fixed_length_{fixed_length}nt.tsv'
    params:
        random_seed = 42,
        species_dict = config['species_short_name'],
        chr_size_dir = 'data/references/chr_size/',
        data_aug_num = 5  # how many nt on each side do we want to increase sno length
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/tuning_train_test_split_rfam_sno_pseudo_fixed_length_data_aug.py"

rule tuning_train_test_split_rfam_sno_pseudo_fixed_length_data_aug_equal_ratio:
    """ Create sets for the training of the second model (expressed sno vs pseudo).
        Split expressed C/D and pseudogenes in 3 datasets (tuning (10%), training (70%) 
        and test set (20%)). SnoRNAs of a same Rfam clan (then Rfam family) 
        are all kept within the same set so that we limit overfitting. After split, do 
        data augmentation by taking 15 windows before and 15 after the actual pseudosno 
        windows and 5 windows for expressed C/D (so to keep the ~ same number of examples per class)."""
    input:
        tuning = rules.tuning_train_test_split_rfam_sno_pseudo_fixed_length.output.tuning,
        training = rules.tuning_train_test_split_rfam_sno_pseudo_fixed_length.output.training,
        test = rules.tuning_train_test_split_rfam_sno_pseudo_fixed_length.output.test,
        genomes = get_all_genomes('data/references/genome_fa/*.fa')
    output:
        tuning = 'data/references/data_augmentation_equal_ratio/cd_rfam_filtered_tuning_set_sno_pseudo_equal_ratio_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/data_augmentation_equal_ratio/cd_rfam_filtered_training_set_sno_pseudo_equal_ratio_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/data_augmentation_equal_ratio/cd_rfam_filtered_test_set_sno_pseudo_equal_ratio_fixed_length_{fixed_length}nt.tsv',
        all_positives = 'data/references/data_augmentation_equal_ratio/cd_rfam_filtered_all_sno_pseudo_equal_ratio_fixed_length_{fixed_length}nt.tsv'
    params:
        random_seed = 42,
        species_dict = config['species_short_name'],
        chr_size_dir = 'data/references/chr_size/',
        data_aug_num = 5,  # how many nt on each side do we want to increase sno length
        data_aug_num_pseudo = 15
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/tuning_train_test_split_rfam_sno_pseudo_fixed_length_data_aug_equal_ratio.py"