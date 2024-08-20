#rule hypertuning_hmm:
#    """ Define the best hyperparameters of a Hidden Markov model."""
#    input:
#        X_tuning = rules.onehot_encode_normalize_initial_fixed_length.output.normalized_tuning,
#        y_tuning = rules.onehot_encode_normalize_initial_fixed_length.output.target_tuning
#    output:
#        best_hyperparams = 'results/predictions/hmm/hmm_best_params.tsv'
#    params:
#        #hyperparams_search_space = config['hyperparameter_space_hmm'],
#        random_state = 42
#    conda:
#        "../envs/python_new.yaml"
#    script:
#        "../scripts/python/hypertuning_hmm.py"

rule training_hmm:
    """ Train the HMM with the best hyperparams."""
    input:
        #positives = 'data/references/positives/cd_rfam_filtered_all_fixed_length_211nt.tsv',
        #pseudo = 'data/references/negatives/snoRNA_pseudogenes/homo_sapiens_pseudogene_snoRNAs_211nt.tsv',
        #negatives = expand(rules.get_all_initial_negatives.output, **config),
        X_train = expand(rules.get_three_sets_added_features_fixed_length.output.training, **config),
        y_train = expand(rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_training, **config),
        X_test = expand(rules.get_three_sets_added_features_fixed_length.output.test, **config),
        y_test = expand(rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_test, **config)
    output:
        model = expand('results/predictions/hmm/hmm_trained_fold_{fold_num}.pt', 
                                        fold_num=[str(i) for i in range(1,11)])
        #training_metrics_per_fold = 'results/predictions/gru/{fixed_length}nt/gru_training_metrics_{fixed_length}nt_last_epoch_per_fold_w_avg.tsv',
        #learning_curves = expand('results/figures/lineplot/gru/{fixed_length}nt/gru_training_f1_score_{fixed_length}nt_fold_{fold_num}.svg',
        #                                fold_num=[str(i) for i in range(1,11)], allow_missing=True)                                        
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_hmm.py"
