rule split_chr:
    """ Separate a multi-fasta file into separate fasta files per chr and 
        divide into smaller chunks if necessary."""
    input:
        input_fasta = config.get("input_fasta")
    output:
        split_chr_dir = directory("data/references/genome_fa/")
    params:
        chunks = config.get("chunks"),
        chunk_size = config.get("chunk_size"),
        python_script = 'scripts/python/split_chr.py'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/split_chr.py"

rule genome_prediction:
    """ Predict with SnoBIRD's first model C/D box snoRNA genes in the input 
        fasta(s)."""
    input:
        genome_dir = rules.split_chr.output.split_chr_dir,
        snoBIRD = rules.download_models.output.model1,
        pretrained_model = rules.download_DNA_BERT.output.dnabert,
        tokenizer = rules.download_DNA_BERT.output.tokenizer
    output:
        windows = "results/predictions/first_model/positive_windows_{chr_}.tsv"
    params:
        step_size = config.get('step_size'),
        fixed_length = config.get('fixed_length'),
        strand = config.get('strand'),
        python_script = 'scripts/python/genome_prediction.py'
    #conda:
    #    "../envs/python_new.yaml"
    shell:
        "bash scripts/bash/genome_prediction.sh "
        "{input.snoBIRD} {input.genome_dir}/fasta/{wildcards.chr_}.fa "
        "{input.pretrained_model} {input.tokenizer} "
        "{output.windows} "
        "{params.fixed_length} {params.step_size} "
        "{params.strand} {params.python_script}"


rule merge_filter_windows:
    """ From the positive windows predicted by the first model of SnoBIRD, 
        concat these windows in one file (across chr and/or chr chunks), filter 
        the windows by score and consecutive blocks of nt and merge afterwards 
        into one block per prediction.
        Step_size=5 is the best option for S. pombe."""
    input:
        #preds = expand(rules.genome_prediction.output.windows, chr_=CHR_),
        input_fasta = config.get("input_fasta")
    output:
        filtered_preds = 'results/predictions/snoBIRD/schizosaccharomyces_pombe/filtered_preds_step5.bed',
        center_preds = 'results/predictions/snoBIRD/schizosaccharomyces_pombe/filtered_center_preds_step5.bed',
        density_block_length = 'results/figures/density/snoBIRD/pred_block_length_s_pombe.svg',
        bed_overlap_sno = 'results/predictions/snoBIRD/schizosaccharomyces_pombe/overlap_snoBIRD_annotated_CD.bed'
    params:
        fixed_length = 194
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/pombe_filter_windows.py"


'''
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
'''