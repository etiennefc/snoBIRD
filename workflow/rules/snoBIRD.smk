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
        windows = "results/intermediate/predictions/first_model/positive_windows_{chr_}.tsv"
    params:
        step_size = config.get('step_size'),
        fixed_length = config.get('fixed_length'),
        strand = config.get('strand'),
        python_script = 'scripts/python/genome_prediction.py',
        gpu = config.get('gpu_generation'),
        chunks = config.get("chunks"),
        chunk_size = config.get("chunk_size"),
        chr_dict = config.get('CHR_dict'),
        batch_size = config.get("batch_size"),
        num_labels = config.get('num_labels')
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
        into one block per prediction (and one center window surrounding the 
        snoRNA)."""
    input:
        predictions = expand(rules.genome_prediction.output.windows, 
                            chr_=config.get('CHR_')),
        input_fasta_dir = rules.split_chr.output.split_chr_dir,
        input_fasta = config.get("input_fasta")
    output:
        filtered_preds = 'results/intermediate/predictions/first_model/filtered_positive_windows.bed',
        center_preds = 'results/intermediate/predictions/first_model/filtered_center_positive_windows.bed'
    params:
        fixed_length = config.get('fixed_length'),
        step_size = config.get("step_size"),
        chunk_size = config.get("chunk_size"),
        strand = config.get("strand"),
        prob_threshold = config.get("min_probability_threshold_first_model"),
        min_consecutive_windows_threshold = config.get("min_consecutive_windows_threshold")
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/merge_filter_windows.py"

rule shap_snoBIRD:
    """ Compute SHAP values for each predicted C/D snoRNA. This is needed to 
        find the C and D boxes and define the snoRNA limits (start/end) in the 
        positive window that encompass the snoRNA. The SHAP values give an 
        indication of which part of the sequence is important for the C/D 
        snoRNA prediction (at the nucleotide resolution)."""
    input:
        snoBIRD = rules.download_models.output.model1,
        preds = rules.merge_filter_windows.output.center_preds,
        pretrained_model = rules.download_DNA_BERT.output.dnabert,
        tokenizer = rules.download_DNA_BERT.output.tokenizer
    output:
        shap_df = "results/intermediate/predictions/first_model/SHAP/shap_values_all_predictions.tsv"
    params:
        fixed_length = config.get('fixed_length'),
        python_script = 'scripts/python/shap_snoBIRD.py',
        gpu = config.get('gpu_generation'),
        batch_size = config.get("batch_size"),
        num_labels = config.get('num_labels')
    #conda:
    #    "../envs/python_new.yaml"
    shell:
        "bash scripts/bash/shap_snoBIRD.sh "
        "{input.snoBIRD} {input.preds} "
        "{input.pretrained_model} {input.tokenizer} "
        "{params.fixed_length} {params.python_script} "
        "{output.shap_df}"

rule find_sno_limits_shap_minimal:
    """ Run this rule if the user wants to run ONLY the first SnoBIRD model 
        (minimal case in which the pipeline stops when the C/D snoRNA 
        predictions are complete (no info on expressed vs pseudogene)). 
        Based on the SHAP values, find the C and/or the D box to delimit the 
        snoRNA start and end, as well as the C' and D' boxes. Return also the 
        box score (i.e. sum of mutations across the C/D/C'/D' boxes). """
    input:
        shap_df = rules.shap_snoBIRD.output.shap_df,
        preds = rules.merge_filter_windows.output.center_preds
    output:
        minimal_output = 'results/final/snoBIRD_complete_predictions.{output_type}'
    params:
        output_type = config.get("output_type"),
        fixed_length = config.get("fixed_length"),
        min_box_dist = config.get("min_box_dist"),
        flanking_nt = config.get("flanking_nt")
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/find_sno_limits_shap_minimal.py"

rule find_sno_limits_shap:
    """ Based on the SHAP values, find the C and/or the D box to delimit the 
        snoRNA start and end, as well as the C' and D' boxes. This creates the 
        dataframe needed for the second step of SnoBIRD predictions"""
    input:
        shap_df = rules.shap_snoBIRD.output.shap_df,
        preds = rules.merge_filter_windows.output.center_preds
    output:
        df = 'results/intermediate/predictions/first_model/SHAP/all_cd_predicted_sno_limits.tsv'
    params:
        fixed_length = config.get("fixed_length"),
        min_box_dist = config.get("min_box_dist"),
        flanking_nt = config.get("flanking_nt")
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/find_sno_limits_shap.py"

rule sno_pseudo_prediction:
    """ Predict with SnoBIRD's second model if the C/D box snoRNA genes 
        identified by the first model are expressed or pseudogenes."""
    input:
        snoBIRD = rules.download_models.output.model2,
        pretrained_model = rules.download_DNA_BERT.output.dnabert,
        tokenizer = rules.download_DNA_BERT.output.tokenizer,
        preds = rules.merge_filter_windows.output.center_preds
    output:
        windows = "results/intermediate/predictions/second_model/sno_pseudo_predictions.tsv"
    params:
        fixed_length = config.get('fixed_length'),
        python_script = 'scripts/python/sno_pseudo_prediction.py',
        gpu = config.get('gpu_generation'),
        batch_size = config.get("batch_size"),
        num_labels = config.get('num_labels')
    #conda:
    #    "../envs/python_new.yaml"
    shell:
        "bash scripts/bash/sno_pseudo_prediction.sh "
        "{input.pretrained_model} {input.tokenizer} "
        "{input.preds} {input.snoBIRD} "
        "{output.windows} "
        "{params.fixed_length} {params.python_script}"

rule filter_sno_pseudo_predictions_with_features:
    """ Compute the snoRNA normalized structure stability as well as its 
        terminal stem combined score. Use these metrics with the box score and 
        the second model's predictions to filter even more which C/D snoRNAs 
        are predicted as expressed or snoRNA pseudogenes. """
    input:
        sno_limits = rules.find_sno_limits_shap.output.df,
        sno_pseudo_preds = rules.sno_pseudo_prediction.output.windows
    output:
        final_output = 'results/final/snoBIRD_complete_predictions.{output_type}'
    params:
        output_type = config.get("output_type"),
        prob_second_model = config.get('min_probability_threshold_second_model'),
        box_score_thresh = config.get("box_score_threshold"),
        score_c_thresh = config.get("score_c_threshold"),
        score_d_thresh = config.get("score_d_threshold"),
        terminal_stem_score_thresh = config.get("terminal_stem_score_threshold"),
        normalized_sno_stability_thresh = config.get("normalized_sno_stability_threshold")
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_sno_pseudo_predictions_with_features.py"
