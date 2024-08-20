rule hypertuning_lstm:
    """ Define the best hyperparameters of a long short-term memory (LSTM), 
        a type of recurrent neural network. By default, we use the 
        Tree-structured Parzen estimator (TPE), which is a bayesian optimization 
        method that is in general more efficent and effective than grid search.
        We use a 3-fold cross-validation to evaluate the f1-score to find the 
        hyperparams combination which maximizes the f1-score (which accounts for 
        all 3 classes equally). The input features are the one-hot encoded sequence 
        and normalized intrinsic features (box_score, structure_mfe, 
        terminal_stem_mfe and length)."""
    input:
        X_tuning = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_tuning,
        y_tuning = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_tuning
    output:
        best_hyperparams = 'results/predictions/lstm/{fixed_length}nt/lstm_best_params_{fixed_length}nt.tsv'
    params:
        hyperparams_search_space = config['hyperparameter_space_LSTM_complex'],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/hypertuning_lstm.py"

rule lr_boundaries_lstm:
    """ Define the learning rate boundaries."""
    input:
        X_train = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_training,
        y_train = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_training,
        best_hyperparams = rules.hypertuning_lstm.output.best_hyperparams
    output:
        lr_boundaries = 'results/predictions/lstm/{fixed_length}nt/lstm_learning_rate_boundaries_{fixed_length}nt.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/lr_boundaries_lstm.py"

rule training_lstm:
    """ Train the LSTM with the best hyperparams. Use a 10-fold CV 
        approach to evaluate the average performance of the model. 
        Save the model trained on each fold and its performance 
        metrics. Save also the average metrics across the 10 folds."""
    input:
        X_train = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_training,
        y_train = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_training,
        best_hyperparams = rules.hypertuning_lstm.output.best_hyperparams  # 2 layers (135, 191 nodes), dropout=0.508, lr=0.00279, Adam
        #best_hyperparams = "lstm_best_params_211nt_MODIFIED.tsv"
    output:
        trained_model = expand('results/predictions/lstm/{fixed_length}nt/lstm_trained_{fixed_length}nt_fold_{fold_num}.pt', 
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        training_metrics_per_fold = 'results/predictions/lstm/{fixed_length}nt/lstm_training_metrics_{fixed_length}nt_last_epoch_per_fold_w_avg.tsv',
        learning_curves = expand('results/figures/lineplot/lstm/{fixed_length}nt/lstm_training_f1_score_{fixed_length}nt_fold_{fold_num}.svg',
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        all_fold_epochs_df = 'results/predictions/lstm/{fixed_length}nt/lstm_training_f1_score_{fixed_length}nt_per_epoch_per_fold.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_lstm.py"

rule training_lstm_cyclic_lr:
    """ Train the LSTM with the best hyperparams. Use a 10-fold CV 
        approach to evaluate the average performance of the model. 
        Save the model trained on each fold and its performance 
        metrics. Save also the average metrics across the 10 folds. 
        Use a cyclic learning rate approach to overcome plateaus in 
        performance"""
    input:
        X_train = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_training,
        y_train = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_training,
        best_hyperparams = rules.hypertuning_lstm.output.best_hyperparams
    output:
        trained_model = expand('results/predictions/lstm/{fixed_length}nt/lstm_cyclic_lr_trained_{fixed_length}nt_fold_{fold_num}.pt', 
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        training_metrics_per_fold = 'results/predictions/lstm/{fixed_length}nt/lstm_cyclic_lr_training_metrics_{fixed_length}nt_last_epoch_per_fold_w_avg.tsv',
        learning_curves = expand('results/figures/lineplot/lstm/{fixed_length}nt/lstm_cyclic_lr_training_f1_score_{fixed_length}nt_fold_{fold_num}.svg',
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        all_fold_epochs_df = 'results/predictions/lstm/{fixed_length}nt/lstm_cyclic_lr_training_f1_score_{fixed_length}nt_per_epoch_per_fold.tsv'
    params:
        random_state = 42,
        low_lr = 0.0005,
        high_lr = 0.02
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_lstm_cyclic_lr.py"

rule test_lstm:
    """ Test the performance of the trained LSTM on the actual test 
        set. Don't forget to use the best hyperparams (and specify
        the model architecture before loading the weights/parameters
        learned during training)."""
    input:
        X_test = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_test,
        y_test = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_test,
        best_hyperparams = rules.hypertuning_lstm.output.best_hyperparams,
        #best_hyperparams = "lstm_best_params_211nt_MODIFIED.tsv",
        training_metrics = rules.training_lstm.output.training_metrics_per_fold,
        all_models = expand(rules.training_lstm.output.trained_model, 
                            fold_num=[str(i) for i in range(1,11)], allow_missing=True)
    output:
        df_metrics_on_test = 'results/predictions/lstm/{fixed_length}nt/lstm_test_metrics_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/test_lstm.py"

