#rule separate_chrom:
#rule window_step_analysis:
#    """ Find what is the optimal step size based on the number of consecutive positive windows predicted as C/D snoRNAs based on the correctly predicted snoRNAs (expressed and pseudogenes)."""
#    input:
#        snoBIRD = 'data/references/trained_transformer_2_classes_4e-5_4e-7_16_30/transformer_2_classes_LR_schedule_trained_fold_9.pt',
#        genome = 'data/references/genome_fa/',
#        preds = 'data/references/trained_transformer_2_classes_4e-5_4e-7_16_30/transformer_2_classes_LR_schedule_test_predictions_190nt_fold_9.tsv',
#        sno_pos = 'data/references/cd_rfam_filtered_all_fixed_length_190nt.tsv',
#        pseudo_pos = 'data/references/human_mouse_pseudogene_snoRNAs_190nt.tsv'
#    output:
#        window_preds = 'results/predictions/transformer/window_step_analysis/prediction_50nt_up_downstream_positives.tsv'
#    params:
#        random_state = 42,
#        pretrained_model = "zhihan1996/DNA_bert_6",
#        fixed_length = 190,
#        step_size = 1,
#        sp_name_dict = config['species_short_name']
#    conda:
#        "../envs/python_new2.yaml"
#    script:
#        "../scripts/python/window_step_analysis.py"
#
#rule window_prediction_dist:
#    """ Create a lineplot that shows the proportion of predicted windows 
#        as positives 50 nt before and after the actaul snoRNA window 
#        (for accurately predicted CD snoRNAs in the test set)."""
#    input:
#        df = rules.window_step_analysis.output.window_preds
#    output:
#        lineplot = 'results/figures/lineplot/transformer/190nt/window_prediction_dist.svg'
#    conda:
#        "../envs/python_new.yaml"
#    script:
#        "../scripts/python/figures/window_prediction_dist.py"


rule genome_windows_separate_chrom:
    """ Predict windows of 194 nt on the genome."""
    input:
        snoBIRD = 'data/references/trained_transformer_2_classes_4e-5_4e-7_16_30/transformer_2_classes_LR_schedule_trained_fold_9.pt',
        genome = 'data/references/genome_fa/drosophila_melanogaster/{chr}.fa'  # rule separate_chrom
    output:
        windows = 'results/predictions/snoBIRD/drosophila_melanogaster/pos_windows_{chr}.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        fixed_length = 190,
        step_size = 1,
	strand = 'both',
	python_script = "scripts/python/genome_windows_separate_chrom.py"
    shell:
        "bash scripts/bash/genome_windows_separate_chrom.sh "
	"{input.snoBIRD} {input.genome} {output.windows} "
	"{params.pretrained_model} "
	"{params.fixed_length} {params.step_size} "
	"{params.strand} {params.python_script}"

rule genome_windows_separate_chrom_1X:
    """ Predict windows of 194 nt on the genome."""
    input:
        snoBIRD = 'data/references/trained_transformer_2_classes_2e-5_2e-6_16_20_1X/transformer_2_classes_LR_schedule_trained_fold_5.pt',
        genome = 'data/references/genome_fa/saccharomyces_cerevisiae/{chr}.fa'  
    output:
        windows = 'results/predictions/snoBIRD/saccharomyces_cerevisiae/pos_windows_1X_{chr}.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        fixed_length = 190,
        step_size = 1,
        strand = 'both',
        python_script = "scripts/python/genome_windows_separate_chrom.py"
    shell:
        "bash scripts/bash/genome_windows_separate_chrom.sh "
        "{input.snoBIRD} {input.genome} {output.windows} "
        "{params.pretrained_model} "
        "{params.fixed_length} {params.step_size} "
        "{params.strand} {params.python_script}"

rule genome_windows_separate_chrom_2X:
    """ Predict windows of 194 nt on the genome."""
    input:
        snoBIRD = 'data/references/trained_transformer_2_classes_2e-5_2e-6_16_20_2X/transformer_2_classes_LR_schedule_trained_fold_10.pt',
        genome = 'data/references/genome_fa/saccharomyces_cerevisiae/{chr}.fa'
    output:
        windows = 'results/predictions/snoBIRD/saccharomyces_cerevisiae/pos_windows_2X_{chr}.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        fixed_length = 190,
        step_size = 1,
        strand = 'both',
        python_script = "scripts/python/genome_windows_separate_chrom.py"
    shell:
        "bash scripts/bash/genome_windows_separate_chrom.sh "
        "{input.snoBIRD} {input.genome} {output.windows} "
        "{params.pretrained_model} "
        "{params.fixed_length} {params.step_size} "
        "{params.strand} {params.python_script}"

rule genome_windows_separate_chrom_data_aug:
    """ Predict windows of 194 nt on the genome."""
    input:
        snoBIRD = 'data/references/trained_transformer_2_classes_2e-5_2e-6_16_4_data_aug/transformer_2_classes_LR_schedule_trained_fold_6.pt',
        genome = 'data/references/genome_fa/saccharomyces_cerevisiae/{chr}.fa'  
    output:
        windows = 'results/predictions/snoBIRD/saccharomyces_cerevisiae/pos_windows_data_aug_fold6_{chr}.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        fixed_length = 190,
        step_size = 1,
        strand = 'both',
        python_script = "scripts/python/genome_windows_separate_chrom.py"
    shell:
        "bash scripts/bash/genome_windows_separate_chrom.sh "
        "{input.snoBIRD} {input.genome} {output.windows} "
        "{params.pretrained_model} "
        "{params.fixed_length} {params.step_size} "
        "{params.strand} {params.python_script}"

rule genome_windows_separate_chrom_data_aug2:
    """ Predict windows of 194 nt on the genome."""
    input:
        snoBIRD = 'data/references/trained_transformer_2_classes_2e-5_2e-6_16_4_data_aug/transformer_2_classes_LR_schedule_trained_fold_1.pt',
        genome = 'data/references/genome_fa/saccharomyces_cerevisiae/{chr}.fa'
    output:
        windows = 'results/predictions/snoBIRD/saccharomyces_cerevisiae/pos_windows_data_aug_fold1_{chr}.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        fixed_length = 190,
        step_size = 1,
        strand = 'both',
        python_script = "scripts/python/genome_windows_separate_chrom.py"
    shell:
        "bash scripts/bash/genome_windows_separate_chrom.sh "
        "{input.snoBIRD} {input.genome} {output.windows} "
        "{params.pretrained_model} "
        "{params.fixed_length} {params.step_size} "
        "{params.strand} {params.python_script}"

rule genome_windows_separate_chrom_equal:
    """ Predict windows of 194 nt on the genome."""
    input:
        snoBIRD = 'data/references/trained_transformer_2_classes_2e-5_2e-6_16_20_equal/transformer_2_classes_LR_schedule_trained_fold_3.pt',
        genome = 'data/references/genome_fa/saccharomyces_cerevisiae/{chr}.fa'
    output:
        windows = 'results/predictions/snoBIRD/saccharomyces_cerevisiae/pos_windows_equal_{chr}.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        fixed_length = 190,
        step_size = 1,
        strand = 'both',
        python_script = "scripts/python/genome_windows_separate_chrom.py"
    shell:
        "bash scripts/bash/genome_windows_separate_chrom.sh "
        "{input.snoBIRD} {input.genome} {output.windows} "
        "{params.pretrained_model} "
        "{params.fixed_length} {params.step_size} "
        "{params.strand} {params.python_script}"


#rule genome_windows_separate_chrom_onnx:
#    """ Create windows of 190 the Candida genome."""
#    input:
#        snoBIRD = 'data/references/trained_transformer_2_classes_4e-5_4e-7_16_30/transformer_2_classes_LR_schedule_trained_fold_9.pt',
#        genome = 'data/references/genome_fa/candida_albicans/{chr}.fa'  # rule separate_chrom
#    output:
#        windows = 'results/onnx/windows_onnx_{chr}.tsv'
#    params:
#        random_state = 42,
#        pretrained_model = "zhihan1996/DNA_bert_6",
#        fixed_length = 190,
#        step_size = 1,
#        strand = 'both'
#    conda:
#        "../envs/python_new3.yaml"
#    script:
#        "../scripts/python/genome_windows_separate_chrom_onnx.py"
