# Try with a simple SVM, log_reg and random forest to see if they're capable of high accuracy or not

rule hyperparameter_tuning_cv_simple_models:
    """ Tune the hyperparameters of each model (Logistic regression, Support
        vector classifier, Random Forest, Gradient boosting classifier 
        and K-nearest neighbors) before even training them, using GridSearchCV 
        with stratified k-fold. The CV will take place as a stratifed k-fold 
        (3 fold) using GridSearchCV and will return the best hyperparameters 
        for each tested model. Tune to optimize based on the f1_score metrics 
        (unweighted average between the 3 classes, because we want it to be as 
        good on all classes, not more on the 'other' class!"""
    input:
        X_cv = rules.onehot_encode_normalize_initial_fixed_length.output.normalized_tuning,
        y_cv = rules.onehot_encode_normalize_initial_fixed_length.output.target_tuning
    output:
        best_hyperparameters = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_best_params_{fixed_length}nt.tsv'
    params:
        hyperparameter_space = lambda wildcards: config['hyperparameter_space'][wildcards.simple_models],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/hyperparameter_tuning_cv_simple_models.py"

rule train_simple_models:
    """ Train (fit) each model on the training set using the best
        hyperparameters found by hyperparameter_tuning_cv_simple_models.
        Pickle these fitted models (into .sav files) so that they can be reused
        after without the need to retrain them all over again."""
    input:
        X_train = rules.onehot_encode_normalize_initial_fixed_length.output.normalized_training,
        y_train = rules.onehot_encode_normalize_initial_fixed_length.output.target_training,
        best_hyperparameters = rules.hyperparameter_tuning_cv_simple_models.output.best_hyperparameters
    output:
        pickled_trained_model = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_trained_{fixed_length}nt.sav',
        training_accuracy = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_training_accuracy_{fixed_length}nt.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/train_simple_models.py"

rule test_simple_models:
    """ Test model performance on unseen test data and return their accuracies. Also return 
        the actual predictions, as well as the performance on the snoRNA pseudogenes only."""
    input:
        X_test = rules.onehot_encode_normalize_initial_fixed_length.output.normalized_test,
        y_test = rules.onehot_encode_normalize_initial_fixed_length.output.target_test,
        pickled_trained_model = rules.train_simple_models.output.pickled_trained_model
    output:
        test_accuracy = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_test_accuracy_{fixed_length}nt.tsv',
        y_preds = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_test_predictions_{fixed_length}nt.tsv',
        pseudosno_performance = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_performance_on_snoRNA_pseudogene_{fixed_length}nt.tsv',
        pseudosno_preds = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_pseudosno_predictions_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/test_simple_models.py"

rule hyperparameter_tuning_cv_simple_models_4_features:
    """ Tune the hyperparameters of each model (Logistic regression, Support
        vector classifier, Random Forest, Gradient boosting classifier 
        and K-nearest neighbors) before even training them, using GridSearchCV 
        with stratified k-fold. The CV will take place as a stratifed k-fold 
        (3 fold) using GridSearchCV and will return the best hyperparameters 
        for each tested model. Tune to optimize based on the f1_score metrics 
        (unweighted average between the 3 classes, because we want it to be as 
        good on all classes, not more on the 'other' class. Use only the 4 intrinsic features."""
    input:
        X_cv = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_tuning,
        y_cv = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_tuning
    output:
        best_hyperparameters = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_best_params_4_features_{fixed_length}nt.tsv'
    params:
        hyperparameter_space = lambda wildcards: config['hyperparameter_space'][wildcards.simple_models],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/hyperparameter_tuning_cv_simple_models_4_features.py"

rule train_simple_models_4_features:
    """ Train (fit) each model on the training set using the best
        hyperparameters found by hyperparameter_tuning_cv_simple_models.
        Pickle these fitted models (into .sav files) so that they can be reused
        after without the need to retrain them all over again. Train on 4 intrinsic features only"""
    input:
        X_train = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_training,
        y_train = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_training,
        best_hyperparameters = rules.hyperparameter_tuning_cv_simple_models_4_features.output.best_hyperparameters
    output:
        pickled_trained_model = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_trained_4_features_{fixed_length}nt.sav',
        training_accuracy = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_training_accuracy_4_features_{fixed_length}nt.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/train_simple_models_4_features.py"

rule test_simple_models_4_features:
    """ Test model performance on unseen test data and return their accuracies. Also return 
        the actual predictions, as well as the performance on the snoRNA pseudogenes only. 
        On 4 intrinsic features only."""
    input:
        X_test = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.encoded_test,
        y_test = rules.onehot_encode_normalize_added_features_half_normalized_fixed_length.output.target_test,
        pickled_trained_model = rules.train_simple_models_4_features.output.pickled_trained_model
    output:
        test_accuracy = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_test_accuracy_4_features_{fixed_length}nt.tsv',
        y_preds = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_test_predictions_4_features_{fixed_length}nt.tsv',
        pseudosno_performance = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_performance_on_snoRNA_pseudogene_4_features_{fixed_length}nt.tsv',
        pseudosno_preds = 'results/predictions/{simple_models}/{fixed_length}nt/{simple_models}_pseudosno_predictions_4_features_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/test_simple_models_4_features.py"

#### return the performance on the snoRNA pseudogenes
#### return confusion matrix as with existing_cd_predictors