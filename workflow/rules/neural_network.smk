rule hypertuning_gru:
    """ Define the best hyperparameters of a gated recurrent units (GRU), 
        a type of recurrent neural network. By default, we use the 
        Tree-structured Parzen estimator (TPE), which is a bayesian optimization 
        method that is in general more efficent and effective than grid search.
        We use a 3-fold cross-validation to evaluate the f1-score to find the 
        hyperparams combination which maximizes the f1-score (which accounts for 
        all 3 classes equally)."""
    input:
        X_tuning = rules.onehot_encode_normalize_initial_fixed_length.output.normalized_tuning,
        y_tuning = rules.onehot_encode_normalize_initial_fixed_length.output.target_tuning
    output:
        best_hyperparams = 'results/predictions/gru/{fixed_length}nt/gru_best_params_{fixed_length}nt.tsv'
    params:
        hyperparams_search_space = config['hyperparameter_space_GRU'],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/hypertuning_gru.py"

rule hypertuning_gru_added_features:
    """ Define the best hyperparameters of a gated recurrent units (GRU), 
        a type of recurrent neural network. By default, we use the 
        Tree-structured Parzen estimator (TPE), which is a bayesian optimization 
        method that is in general more efficent and effective than grid search.
        We use a 3-fold cross-validation to evaluate the f1-score to find the 
        hyperparams combination which maximizes the f1-score (which accounts for 
        all 3 classes equally). The input features are the sequence and 4
        intrinsic features (box_score, structure_mfe, terminal_stem_mfe and length)."""
    input:
        X_tuning = rules.onehot_encode_normalize_added_features_fixed_length.output.normalized_tuning,
        y_tuning = rules.onehot_encode_normalize_added_features_fixed_length.output.target_tuning
    output:
        best_hyperparams = 'results/predictions/gru/{fixed_length}nt/added_features/gru_best_params_{fixed_length}nt.tsv'
    params:
        hyperparams_search_space = config['hyperparameter_space_GRU'],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/hypertuning_gru.py"

rule hypertuning_gru_added_features_wo_seq:
    """ Define the best hyperparameters of a gated recurrent units (GRU), 
        a type of recurrent neural network. By default, we use the 
        Tree-structured Parzen estimator (TPE), which is a bayesian optimization 
        method that is in general more efficent and effective than grid search.
        We use a 3-fold cross-validation to evaluate the f1-score to find the 
        hyperparams combination which maximizes the f1-score (which accounts for 
        all 3 classes equally). The input features are ONLY the 4
        intrinsic features (box_score, structure_mfe, terminal_stem_mfe and length)."""
    input:
        X_tuning = rules.onehot_encode_normalize_added_features_fixed_length.output.normalized_tuning,
        y_tuning = rules.onehot_encode_normalize_added_features_fixed_length.output.target_tuning
    output:
        best_hyperparams = 'results/predictions/gru/{fixed_length}nt/added_features/gru_best_params_{fixed_length}nt_wo_seq.tsv'
    params:
        hyperparams_search_space = config['hyperparameter_space_GRU'],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/hypertuning_gru_wo_seq.py"

rule hypertuning_gru_added_features_simplified:
    """ Define the best hyperparameters of a gated recurrent units (GRU), 
        a type of recurrent neural network. By default, we use the 
        Tree-structured Parzen estimator (TPE), which is a bayesian optimization 
        method that is in general more efficent and effective than grid search.
        We use a 3-fold cross-validation to evaluate the f1-score to find the 
        hyperparams combination which maximizes the f1-score (which accounts for 
        all 3 classes equally). The input features are the sequence and 4
        intrinsic features (box_score, structure_mfe, terminal_stem_mfe and length). 
        We use here a simplified hyperparams search space (smaller number of layers, 
        nodes, and higher learning rate and dropout rate)."""
    input:
        X_tuning = rules.onehot_encode_normalize_added_features_fixed_length.output.normalized_tuning,
        y_tuning = rules.onehot_encode_normalize_added_features_fixed_length.output.target_tuning
    output:
        best_hyperparams = 'results/predictions/gru/{fixed_length}nt/added_features/gru_best_params_{fixed_length}nt_simplified.tsv'
    params:
        hyperparams_search_space = config['hyperparameter_space_GRU_simplified'],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/hypertuning_gru.py"

rule hypertuning_gru_added_features_simplified2:
    """ Define the best hyperparameters of a gated recurrent units (GRU), 
        a type of recurrent neural network. By default, we use the 
        Tree-structured Parzen estimator (TPE), which is a bayesian optimization 
        method that is in general more efficent and effective than grid search.
        We use a 3-fold cross-validation to evaluate the f1-score to find the 
        hyperparams combination which maximizes the f1-score (which accounts for 
        all 3 classes equally). The input features are the sequence and 4
        intrinsic features (box_score, structure_mfe, terminal_stem_mfe and length). 
        We use here a simplified hyperparams search space (smaller number of layers, 
        nodes, and higher learning rate and dropout rate)."""
    input:
        X_tuning = rules.onehot_encode_normalize_added_features_fixed_length.output.normalized_tuning,
        y_tuning = rules.onehot_encode_normalize_added_features_fixed_length.output.target_tuning
    output:
        best_hyperparams = 'results/predictions/gru/{fixed_length}nt/added_features/gru_best_params_{fixed_length}nt_simplified2.tsv'
    params:
        hyperparams_search_space = config['hyperparameter_space_GRU_simplified2'],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/hypertuning_gru.py"

rule hypertuning_gru_added_features_half_normalized_simplified:
    """ Define the best hyperparameters of a gated recurrent units (GRU), 
        a type of recurrent neural network. By default, we use the 
        Tree-structured Parzen estimator (TPE), which is a bayesian optimization 
        method that is in general more efficent and effective than grid search.
        We use a 3-fold cross-validation to evaluate the f1-score to find the 
        hyperparams combination which maximizes the f1-score (which accounts for 
        all 3 classes equally). The input features are the sequence and 4
        intrinsic features (box_score, structure_mfe, terminal_stem_mfe and length). 
        We use here a simplified hyperparams search space (smaller number of layers, 
        nodes, and higher learning rate and dropout rate). The sequence features 
        are NOT NORMALIZED (i.e. only 0 or 1, no standardization)."""
    input:
        X_tuning = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_tuning,
        y_tuning = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_tuning
    output:
        best_hyperparams = 'results/predictions/gru/{fixed_length}nt/added_features_half_normalized/gru_best_params_{fixed_length}nt_simplified.tsv'
    params:
        hyperparams_search_space = config['hyperparameter_space_GRU_simplified'],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/hypertuning_gru.py"

rule hypertuning_gru_added_features_half_normalized_simplified_ATCG:
    """ Define the best hyperparameters of a gated recurrent units (GRU), 
        a type of recurrent neural network. By default, we use the 
        Tree-structured Parzen estimator (TPE), which is a bayesian optimization 
        method that is in general more efficent and effective than grid search.
        We use a 3-fold cross-validation to evaluate the f1-score to find the 
        hyperparams combination which maximizes the f1-score (which accounts for 
        all 3 classes equally). The input features are the sequence (ONLY ATCG, no N)) and 4
        intrinsic features (box_score, structure_mfe, terminal_stem_mfe and length). 
        We use here a simplified hyperparams search space (smaller number of layers, 
        nodes, and higher learning rate and dropout rate). The sequence features 
        are NOT NORMALIZED (i.e. only 0 or 1, no standardization)."""
    input:
        X_tuning = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_tuning,
        y_tuning = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_tuning
    output:
        best_hyperparams = 'results/predictions/gru/{fixed_length}nt/added_features_half_normalized/gru_best_params_{fixed_length}nt_simplified_ATCG.tsv'
    params:
        hyperparams_search_space = config['hyperparameter_space_GRU_simplified'],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/hypertuning_gru_ATCG.py"

rule training_gru:
    """ Train the GRU with the best hyperparams. Use a 10-fold CV 
        approach to evaluate the average performance of the model. 
        Save the model trained on each fold and its performance 
        metrics. Save also the average metrics across the 10 folds."""
    input:
        X_train = rules.onehot_encode_normalize_initial_fixed_length.output.normalized_training,
        y_train = rules.onehot_encode_normalize_initial_fixed_length.output.target_training,
        best_hyperparams = rules.hypertuning_gru.output.best_hyperparams
    output:
        trained_model = expand('results/predictions/gru/{fixed_length}nt/gru_trained_{fixed_length}nt_fold_{fold_num}.pt', 
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        training_metrics_per_fold = 'results/predictions/gru/{fixed_length}nt/gru_training_metrics_{fixed_length}nt_last_epoch_per_fold_w_avg.tsv',
        learning_curves = expand('results/figures/lineplot/gru/{fixed_length}nt/gru_training_f1_score_{fixed_length}nt_fold_{fold_num}.svg',
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True)                                        
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_gru.py"

rule training_gru_added_features:
    """ Train the GRU with the best hyperparams. Use a 10-fold CV 
        approach to evaluate the average performance of the model. 
        Save the model trained on each fold and its performance 
        metrics. Save also the average metrics across the 10 folds.
        The input features are the sequence and 4 intrinsic features 
        (box_score, structure_mfe, terminal_stem_mfe and length)."""
    input:
        X_train = rules.onehot_encode_normalize_added_features_fixed_length.output.normalized_training,
        y_train = rules.onehot_encode_normalize_added_features_fixed_length.output.target_training,
        best_hyperparams = rules.hypertuning_gru_added_features.output.best_hyperparams
    output:
        trained_model = expand('results/predictions/gru/{fixed_length}nt/added_features/gru_trained_{fixed_length}nt_fold_{fold_num}.pt', 
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        training_metrics_per_fold = 'results/predictions/gru/{fixed_length}nt/added_features/gru_training_metrics_{fixed_length}nt_last_epoch_per_fold_w_avg.tsv',
        learning_curves = expand('results/figures/lineplot/gru/{fixed_length}nt/added_features/gru_training_f1_score_{fixed_length}nt_fold_{fold_num}.svg',
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True)                                        
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_gru.py"

rule training_gru_added_features_simplified:
    """ Train the GRU with the best hyperparams. Use a 10-fold CV 
        approach to evaluate the average performance of the model. 
        Save the model trained on each fold and its performance 
        metrics. Save also the average metrics across the 10 folds.
        The input features are the sequence and 4 intrinsic features 
        (box_score, structure_mfe, terminal_stem_mfe and length)."""
    input:
        X_train = rules.onehot_encode_normalize_added_features_fixed_length.output.normalized_training,
        y_train = rules.onehot_encode_normalize_added_features_fixed_length.output.target_training,
        best_hyperparams = rules.hypertuning_gru_added_features_simplified.output.best_hyperparams
    output:
        trained_model = expand('results/predictions/gru/{fixed_length}nt/added_features/gru_simplified_trained_{fixed_length}nt_fold_{fold_num}.pt', 
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        training_metrics_per_fold = 'results/predictions/gru/{fixed_length}nt/added_features/gru_simplified_training_metrics_{fixed_length}nt_last_epoch_per_fold_w_avg.tsv',
        learning_curves = expand('results/figures/lineplot/gru/{fixed_length}nt/added_features/gru_simplified_training_f1_score_{fixed_length}nt_fold_{fold_num}.svg',
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True)
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_gru.py"

rule training_gru_added_features_simplified2:
    """ Train the GRU with the best hyperparams. Use a 10-fold CV 
        approach to evaluate the average performance of the model. 
        Save the model trained on each fold and its performance 
        metrics. Save also the average metrics across the 10 folds.
        The input features are the sequence and 4 intrinsic features 
        (box_score, structure_mfe, terminal_stem_mfe and length). 
        Increased dropout and learning rate."""
    input:
        X_train = rules.onehot_encode_normalize_added_features_fixed_length.output.normalized_training,
        y_train = rules.onehot_encode_normalize_added_features_fixed_length.output.target_training,
        best_hyperparams = rules.hypertuning_gru_added_features_simplified2.output.best_hyperparams
    output:
        trained_model = expand('results/predictions/gru/{fixed_length}nt/added_features/gru_simplified2_trained_{fixed_length}nt_fold_{fold_num}.pt', 
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        training_metrics_per_fold = 'results/predictions/gru/{fixed_length}nt/added_features/gru_simplified2_training_metrics_{fixed_length}nt_last_epoch_per_fold_w_avg.tsv',
        learning_curves = expand('results/figures/lineplot/gru/{fixed_length}nt/added_features/gru_simplified2_training_f1_score_{fixed_length}nt_fold_{fold_num}.svg',
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        all_fold_epochs_df = 'results/predictions/gru/{fixed_length}nt/added_features/gru_simplified2_training_f1_score_{fixed_length}nt_per_epoch_per_fold.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_gru.py"

rule training_gru_added_features_wo_seq:
    """ Train the GRU with the best hyperparams. Use a 10-fold CV 
        approach to evaluate the average performance of the model. 
        Save the model trained on each fold and its performance 
        metrics. Save also the average metrics across the 10 folds.
        The input features are ONLY the 4 intrinsic features 
        (box_score, structure_mfe, terminal_stem_mfe and length)."""
    input:
        X_train = rules.onehot_encode_normalize_added_features_fixed_length.output.normalized_training,
        y_train = rules.onehot_encode_normalize_added_features_fixed_length.output.target_training,
        best_hyperparams = rules.hypertuning_gru_added_features_wo_seq.output.best_hyperparams
    output:
        trained_model = expand('results/predictions/gru/{fixed_length}nt/added_features/gru_trained_{fixed_length}nt_fold_{fold_num}_wo_seq.pt', 
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        training_metrics_per_fold = 'results/predictions/gru/{fixed_length}nt/added_features/gru_training_metrics_{fixed_length}nt_last_epoch_per_fold_w_avg_wo_seq.tsv',
        learning_curves = expand('results/figures/lineplot/gru/{fixed_length}nt/added_features/gru_training_f1_score_{fixed_length}nt_fold_{fold_num}_wo_seq.svg',
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True)                                        
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_gru_wo_seq.py"

rule training_gru_added_features_half_normalized_simplified:
    """ Train the GRU with the best hyperparams. Use a 10-fold CV 
        approach to evaluate the average performance of the model. 
        Save the model trained on each fold and its performance 
        metrics. Save also the average metrics across the 10 folds.
        The input features are the sequence (one-hot encoded NOT 
        normalized) and 4 normalized intrinsic features 
        (box_score, structure_mfe, terminal_stem_mfe and length)."""
    input:
        X_train = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_training,
        y_train = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_training,
        best_hyperparams = rules.hypertuning_gru_added_features_half_normalized_simplified.output.best_hyperparams
	    #best_hyperparams = 'gru_best_params_211nt_MODIFIED.tsv'
    output:
        trained_model = expand('results/predictions/gru/{fixed_length}nt/added_features_half_normalized/gru_simplified_trained_{fixed_length}nt_fold_{fold_num}.pt', 
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        training_metrics_per_fold = 'results/predictions/gru/{fixed_length}nt/added_features_half_normalized/gru_simplified_training_metrics_{fixed_length}nt_last_epoch_per_fold_w_avg.tsv',
        learning_curves = expand('results/figures/lineplot/gru/{fixed_length}nt/added_features_half_normalized/gru_simplified_training_f1_score_{fixed_length}nt_fold_{fold_num}.svg',
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        all_fold_epochs_df = 'results/predictions/gru/{fixed_length}nt/added_features_half_normalized/gru_simplified_training_f1_score_{fixed_length}nt_per_epoch_per_fold.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_gru.py"

rule training_gru_added_features_half_normalized_simplified_ATCG:
    """ Train the GRU with the best hyperparams. Use a 10-fold CV 
        approach to evaluate the average performance of the model. 
        Save the model trained on each fold and its performance 
        metrics. Save also the average metrics across the 10 folds.
        The input features are the sequence (one-hot encoded NOT 
        normalized (ONLY ATCG, no N)) and 4 normalized intrinsic features 
        (box_score, structure_mfe, terminal_stem_mfe and length)."""
    input:
        X_train = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_training,
        y_train = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_training,
        best_hyperparams = rules.hypertuning_gru_added_features_half_normalized_simplified_ATCG.output.best_hyperparams
	    #best_hyperparams = 'gru_best_params_211nt_MODIFIED.tsv'
    output:
        trained_model = expand('results/predictions/gru/{fixed_length}nt/added_features_half_normalized/gru_simplified_trained_{fixed_length}nt_fold_{fold_num}_ATCG.pt', 
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        training_metrics_per_fold = 'results/predictions/gru/{fixed_length}nt/added_features_half_normalized/gru_simplified_training_metrics_{fixed_length}nt_last_epoch_per_fold_w_avg_ATCG.tsv',
        learning_curves = expand('results/figures/lineplot/gru/{fixed_length}nt/added_features_half_normalized/gru_simplified_training_f1_score_{fixed_length}nt_fold_{fold_num}_ATCG.svg',
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        all_fold_epochs_df = 'results/predictions/gru/{fixed_length}nt/added_features_half_normalized/gru_simplified_training_f1_score_{fixed_length}nt_per_epoch_per_fold_ATCG.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_gru_ATCG.py"

rule test_gru:
    """ Test the performance of the trained GRU on the actual test 
        set. Don't forget to use the best hyperparams (and specify
        the model architecture before loading the weights/parameters
        learned during training)."""
    input:
        X_test = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_test,
        y_test = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_test,
        best_hyperparams = rules.hypertuning_gru_added_features_half_normalized_simplified.output.best_hyperparams,
	    #best_hyperparams = 'gru_best_params_211nt_MODIFIED.tsv',
        training_metrics = rules.training_gru_added_features_half_normalized_simplified.output.training_metrics_per_fold,
        all_models = expand(rules.training_gru_added_features_half_normalized_simplified.output.trained_model, 
                            fold_num=[str(i) for i in range(1,11)], allow_missing=True)
    output:
        df_metrics_on_test = 'results/predictions/gru/{fixed_length}nt/added_features_half_normalized/gru_test_metrics_simplifed2_{fixed_length}nt.tsv',
        test_predictions = 'results/predictions/gru/{fixed_length}nt/added_features_half_normalized/gru_test_predictions_simplifed2_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/test_gru.py"

