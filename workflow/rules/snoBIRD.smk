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
        python_script = 'scripts/python/split_chr.py',
        profile = config.get("profile"),
        cluster_env = rules.create_env.output.env
    shell:
        "if [ {params.profile} = local ]; then "
        "conda run -p snoBIRD_env/ "
        "python3 {params.python_script} {input.input_fasta} "
        "{output.split_chr_dir} {params.chunks} {params.chunk_size}; else "
        "bash scripts/bash/split_chr.sh {params.python_script} "
        "{input.input_fasta} {output.split_chr_dir} {params.chunks} "
        "{params.chunk_size} {params.cluster_env}; "
        "fi"

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
        real_chunk_chr_size = lambda wildcards: config['CHR_sizes'][wildcards.chr_],
        chunks = config.get("chunks"),
        chunk_size = config.get("chunk_size"),  # max chunk size (user-defined)
        batch_size = config.get("batch_size"),
        num_labels = config.get('num_labels'),
        cluster_env = rules.create_env.output.env,
        profile = config.get("profile"),
        chr_chunk_name = '{chr_}'
    shell:
        "if [ {params.profile} = local ]; then "
        "conda run -p snoBIRD_env/ "
        "python3 {params.python_script} {input.snoBIRD} "
        "{input.genome_dir}/fasta/{wildcards.chr_}.fa "
        "{input.pretrained_model} {input.tokenizer} {output.windows} "
        "{params.fixed_length} {params.step_size} {params.strand} "
        "{params.batch_size} {params.num_labels} {params.profile}; else "
        "bash scripts/bash/genome_prediction.sh "
        "{input.snoBIRD} {input.genome_dir}/fasta/{wildcards.chr_}.fa "
        "{input.pretrained_model} {input.tokenizer} "
        "{output.windows} "
        "{params.fixed_length} {params.step_size} "
        "{params.strand} {params.python_script} {params.batch_size} "
        "{params.num_labels} {params.cluster_env} {params.profile}; "
        "fi"

rule merge_filter_windows:
    """ From the positive windows predicted by the first model of SnoBIRD, 
        concat these windows in one file (across chr and/or chr chunks), filter 
        the windows by score and consecutive blocks of nt and merge afterwards 
        into one block per prediction (and one center window surrounding the 
        snoRNA). Predict, if step_size>1, if the center window is also 
        predicted as a C/D to filter out even more predictions. For the 
        predictions input, the files paths are joined into 1 string so that 
        this variable can then be expanded adequately later in the script."""
    input:
        predictions = expand(rules.genome_prediction.output.windows, 
                        chr_=config.get('CHR_')),
        input_fasta_dir = rules.split_chr.output.split_chr_dir,
        input_fasta = config.get("input_fasta"),
        pretrained_model = rules.download_DNA_BERT.output.dnabert,
        tokenizer = rules.download_DNA_BERT.output.tokenizer,
        snoBIRD = rules.download_models.output.model1
    output:
        filtered_preds = 'results/intermediate/predictions/first_model/filtered_positive_windows.bed',
        center_preds = 'results/intermediate/predictions/first_model/filtered_center_positive_windows.bed'
    params:
        input_preds = '__PRED_SEP__'.join(
                        expand(rules.genome_prediction.output.windows, 
                        chr_=config.get('CHR_'))),
        fixed_length = config.get('fixed_length'),
        step_size = config.get("step_size"),
        chunk_size = config.get("chunk_size"),
        gpu = config.get('gpu_generation'),
        batch_size = config.get("batch_size"),
        num_labels = config.get('num_labels'),
        prob_threshold = config.get("min_probability_threshold_first_model"),
        min_consecutive_windows_threshold = config.get("min_consecutive_windows_threshold"),
        python_script = "scripts/python/merge_filter_windows.py",
        cluster_env = rules.create_env.output.env,
        profile = config.get("profile"),
        final_output = config.get('output_name')  # used to create a unique job name (
                    # in case multiple snoBIRD instances are run on the same cluster)
    shell:
        "if [ {params.profile} = local ]; then "
        "conda run -p snoBIRD_env/ "
        "python3 {params.python_script} {params.input_preds} "
        "{input.input_fasta_dir} {input.input_fasta} {input.pretrained_model} "
        "{input.tokenizer} {input.snoBIRD} {output.filtered_preds} "
        "{output.center_preds} {params.fixed_length} {params.step_size} "
        "{params.chunk_size} {params.batch_size} {params.num_labels} "
        "{params.prob_threshold} {params.min_consecutive_windows_threshold} "
        "{params.profile}; else "
        "bash scripts/bash/merge_filter_windows.sh "
        "{params.input_preds} {input.input_fasta_dir} {input.input_fasta} "
        "{input.pretrained_model} {input.tokenizer} {input.snoBIRD} "
        "{output.filtered_preds} {output.center_preds} {params.fixed_length} "
        "{params.step_size} {params.chunk_size} {params.batch_size} "
        "{params.num_labels} {params.prob_threshold} "
        "{params.min_consecutive_windows_threshold} {params.python_script} "
        "{params.cluster_env} {params.profile} {params.gpu} "
        "{params.final_output}; "
        "fi"

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
        num_labels = config.get('num_labels'),
        cluster_env = rules.create_env.output.env,
        profile = config.get("profile"),
        final_output = config.get('output_name')  # used to create a unique job name (
                    # in case multiple snoBIRD instances are run on the same cluster)
    shell:
        "if [ {params.profile} = local ]; then "
        "conda run -p snoBIRD_env/ "
        "python3 {params.python_script} {input.snoBIRD} {input.preds} "
        "{input.pretrained_model} {input.tokenizer} {output.shap_df} "
        "{params.fixed_length} {params.batch_size} {params.num_labels} "
        "{params.profile}; else "
        "bash scripts/bash/shap_snoBIRD.sh "
        "{input.snoBIRD} {input.preds} "
        "{input.pretrained_model} {input.tokenizer} {output.shap_df} "
        "{params.fixed_length} {params.python_script} "
        "{params.batch_size} {params.num_labels} {params.cluster_env} "
        "{params.profile}; "
        "fi"

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
        minimal_output = 'results/final/{output_name}.{output_type}'
    params:
        output_type = config.get("output_type"),
        fixed_length = config.get("fixed_length"),
        python_script = "scripts/python/find_sno_limits_shap_minimal.py",
        min_box_dist = config.get("min_box_dist"),
        flanking_nt = config.get("flanking_nt"),
        cluster_env = rules.create_env.output.env,
        profile = config.get("profile")
    shell:
        "if [ {params.profile} = local ]; then "
        "conda run -p snoBIRD_env/ "
        "python3 {params.python_script} {input.shap_df} {input.preds} "
        "{output.minimal_output} {params.output_type} {params.fixed_length} "
        "{params.min_box_dist} {params.flanking_nt}; else "
        "bash scripts/bash/find_sno_limits_shap_minimal.sh "
        "{input.shap_df} {input.preds} "
        "{output.minimal_output} {params.output_type} {params.fixed_length} "
        "{params.python_script} {params.min_box_dist} {params.flanking_nt} "
        "{params.cluster_env}; "
        "fi"

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
        flanking_nt = config.get("flanking_nt"),
        python_script = "scripts/python/find_sno_limits_shap.py",
        cluster_env = rules.create_env.output.env,
        profile = config.get("profile")
    shell:
        "if [ {params.profile} = local ]; then "
        "conda run -p snoBIRD_env/ "
        "python3 {params.python_script} {input.shap_df} {input.preds} "
        "{output.df} {params.fixed_length} "
        "{params.min_box_dist} {params.flanking_nt}; else "
        "bash scripts/bash/find_sno_limits_shap.sh "
        "{input.shap_df} {input.preds} "
        "{output.df} {params.fixed_length} "
        "{params.python_script} {params.min_box_dist} {params.flanking_nt} "
        "{params.cluster_env}; "
        "fi"

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
        num_labels = config.get('num_labels'),
        cluster_env = rules.create_env.output.env,
        profile = config.get("profile")
    shell:
        "if [ {params.profile} = local ]; then "
        "conda run -p snoBIRD_env/ "
        "python3 {params.python_script} {input.snoBIRD} "
        "{input.pretrained_model} {input.tokenizer} {input.preds} "
        "{output.windows} {params.fixed_length} "
        "{params.batch_size} {params.num_labels} {params.profile}; else "
        "bash scripts/bash/sno_pseudo_prediction.sh "
        "{input.snoBIRD} {input.pretrained_model} {input.tokenizer} "
        "{input.preds} {output.windows} "
        "{params.fixed_length} {params.python_script} "
        "{params.batch_size} {params.num_labels} "
        "{params.cluster_env} {params.profile}; "
        "fi"

rule filter_sno_pseudo_predictions_with_features:
    """ Compute the snoRNA normalized structure stability as well as its 
        terminal stem combined score. Use these metrics with the box score and 
        the second model's predictions to filter even more which C/D snoRNAs 
        are predicted as expressed or snoRNA pseudogenes. """
    input:
        sno_limits = rules.find_sno_limits_shap.output.df,
        sno_pseudo_preds = rules.sno_pseudo_prediction.output.windows
    output:
        final_output = 'results/final/{output_name}.{output_type}'
    params:
        output_type = config.get("output_type"),
        fixed_length = config.get("fixed_length"),
        python_script = "scripts/python/filter_sno_pseudo_predictions_with_features.py",
        prob_second_model = config.get('min_probability_threshold_second_model'),
        box_score_thresh = config.get("box_score_threshold"),
        score_c_thresh = config.get("score_c_threshold"),
        score_d_thresh = config.get("score_d_threshold"),
        terminal_stem_score_thresh = config.get("terminal_stem_score_threshold"),
        normalized_sno_stability_thresh = config.get("normalized_sno_stability_threshold"),
        cluster_env = rules.create_env.output.env,
        profile = config.get("profile")
    shell:
        "if [ {params.profile} = local ]; then "
        "conda run -p snoBIRD_env/ "
        "python3 {params.python_script} {input.sno_limits} "
        "{input.sno_pseudo_preds} {output.final_output} "
        "{params.output_type} {params.fixed_length} "
        "{params.prob_second_model} {params.box_score_thresh} "
        "{params.score_c_thresh} {params.score_d_thresh} "
        "{params.terminal_stem_score_thresh} "
        "{params.normalized_sno_stability_thresh}; else "
        "bash scripts/bash/filter_sno_pseudo_predictions_with_features.sh "
        "{input.sno_limits} {input.sno_pseudo_preds} "
        "{output.final_output} {params.output_type} {params.fixed_length} "
        "{params.python_script} {params.prob_second_model} "
        "{params.box_score_thresh} {params.score_c_thresh} "
        "{params.score_d_thresh} {params.terminal_stem_score_thresh} "
        "{params.normalized_sno_stability_thresh} {params.cluster_env}; "
        "fi"
        
