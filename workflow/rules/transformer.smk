rule training_transformer:
    """ Train the Transformer with the best hyperparams. Must be connected to internet to load the pretrained model for the first time"""
    input:
        X_train = 'data/references/positives_and_negatives/added_features/added_features_training_set_fixed_length_{fixed_length}nt.tsv',
        y_train = 'data/references/positives_and_negatives/added_features/added_features_training_target_{fixed_length}nt.tsv'
    output:
        model = 'results/predictions/transformer/{fixed_length}/transformer_trained_fold_{fold_num}.pt',
        fold_loss = 'results/predictions/transformer/{fixed_length}/transformer_trained_fold_{fold_num}_loss_per_epoch.tsv',
        fold_f1_score = 'results/predictions/transformer/{fixed_length}/transformer_trained_fold_{fold_num}_f1_score_per_epoch.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/transformer_training_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_training.sh "
	"{params.pretrained_model} {wildcards.fold_num} "
	"{params.random_state} "
	"{input.X_train} {input.y_train} "
	"{output.model} {output.fold_loss} {output.fold_f1_score} &> {log}"


rule hypertuning_transformer_2_classes:
    """ Get the best hyperparameters of the transformer with 2 classes using Grid Search. With 20 epochs, we get lr=4e-5 and batch_size=16"""
    input:
        X_tuning = rules.get_three_sets_initial_fixed_length.output.tuning,
        y_tuning = rules.get_three_sets_initial_fixed_length.output.tuning_target
    output:
        best_hyperparams = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_best_hyperparams.tsv'
    params:
        hyperparams_space = config['hyperparameter_space_transformer_2_classes'],
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/hypertuning_transformer_2_classes_{fixed_length}nt.log"
    shell:
        "bash scripts/bash/hypertuning_transformer_2_classes.sh "
        "{params.random_state} {params.pretrained_model} "
        "{input.X_tuning} {input.y_tuning} "
        "{output.best_hyperparams} &> {log}"


rule training_transformer_2_classes:
    """ Train the Transformer with the best hyperparams. ONLY 2 classes to predict: sno/pseudosno (1) or other (0). Must be connected to internet to load the pretrained model for the first time. ORGIINAL hyperparams: lr=2e-5, batch_size=100"""
    input:
        X_train = rules.get_three_sets_initial_fixed_length.output.training,
        y_train = rules.get_three_sets_initial_fixed_length.output.training_target,
        best_hyperparams = rules.hypertuning_transformer_2_classes.output.best_hyperparams
    output:
        model = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_trained_fold_{fold_num}.pt',
        fold_loss = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_trained_fold_{fold_num}_loss_per_epoch.tsv',
        fold_f1_score = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_trained_fold_{fold_num}_f1_score_per_epoch.tsv'
    params:
        python_script = "scripts/python/training_2_classes_transformer.py",
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/transformer_2_classes_training_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_2_classes_training.sh "
        "{params.pretrained_model} {wildcards.fold_num} "
        "{params.random_state} "
        "{input.X_train} {input.y_train} {input.best_hyperparams} "
        "{output.model} {output.fold_loss} {output.fold_f1_score} &> {log}"

rule training_transformer_2_classes_LR_schedule:
    """ Train the Transformer with the best hyperparams. ONLY 2 classes to predict: sno/pseudosno (1) or other (0). Must be connected to internet to load the pretrained model for the first time. The learning rate is firt linearly warmed-up from 0 to 4e-5 in the first epoch, and then linearly decayed to 0 at the last epoch of training (as was done in Ji et al., 2021, Bioinformatics)"""
    input:
        #X_train = rules.get_three_sets_initial_fixed_length.output.training,
        #y_train = rules.get_three_sets_initial_fixed_length.output.training_target,
        #X_train = 'data/references/positives_and_negatives/initial/initial_training_set_fixed_length_194nt.tsv',
        #y_train = 'data/references/positives_and_negatives/initial/initial_training_target_fixed_length_194nt.tsv',
        X_train = 'data/references/positives_and_negatives/data_augmentation/training_set_fixed_length_194nt.tsv',
        y_train = 'data/references/positives_and_negatives/data_augmentation/training_target_fixed_length_194nt.tsv',
        #best_hyperparams = rules.hypertuning_transformer_2_classes.output.best_hyperparams
        best_hyperparams = 'results/predictions/transformer/194/transformer_2_classes_best_hyperparams.tsv'
    output:
        model = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_LR_schedule_trained_fold_{fold_num}.pt',
        fold_loss = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_LR_schedule_trained_fold_{fold_num}_loss_per_epoch.tsv',
        fold_f1_score = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_LR_schedule_trained_fold_{fold_num}_f1_score_per_epoch.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/transformer_2_classes_training_LR_schedule_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_2_classes_training_LR_schedule.sh "
        "{params.pretrained_model} {wildcards.fold_num} "
        "{params.random_state} "
        "{input.X_train} {input.y_train} {input.best_hyperparams} "
        "{output.model} {output.fold_loss} {output.fold_f1_score} &> {log}"

rule test_before_training_transformer_2_classes:
    """ Get the Transformer predictions before training (baseline!) with the best hyperparams. ONLY 2 classes to predict: sno/pseudosno (1) or other (0). Must be connected to internet to load the pretrained model for the first time."""
    input:
        X_train = rules.get_three_sets_initial_fixed_length.output.training,        
        y_train = rules.get_three_sets_initial_fixed_length.output.training_target,
        best_hyperparams = rules.hypertuning_transformer_2_classes.output.best_hyperparams
    output:
        fold_loss = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_Before_trained_fold_{fold_num}_loss_per_epoch.tsv',
        fold_f1_score = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_Before_trained_fold_{fold_num}_f1_score_per_epoch.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/test_before_training_transformer_2_classes_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/test_before_training_transformer_2_classes.sh "
        "{params.pretrained_model} {wildcards.fold_num} "
        "{params.random_state} "
        "{input.X_train} {input.y_train} {input.best_hyperparams} "
        "{output.fold_loss} {output.fold_f1_score} &> {log}"

rule training_transformer_w_features:
    """ Train the Transformer with sequence and numerical features. Must be connected to internet to load the pretrained model for the first time"""
    input:
        X_train = 'data/references/positives_and_negatives/added_features/added_features_training_set_fixed_length_{fixed_length}nt.tsv',
        y_train = 'data/references/positives_and_negatives/added_features/added_features_training_target_{fixed_length}nt.tsv'
    output:
        model = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_fold_{fold_num}.pt',
        fold_loss = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_fold_{fold_num}_loss_per_epoch.tsv',
        fold_f1_score = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_fold_{fold_num}_f1_score_per_epoch.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/transformer_w_features/transformer_training_w_features_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_training_w_features.sh "
        "{params.pretrained_model} {wildcards.fold_num} "
        "{params.random_state} "
        "{input.X_train} {input.y_train} "
        "{output.model} {output.fold_loss} {output.fold_f1_score} &> {log}"


rule training_transformer_w_features_2_classes:
    """ Train the Transformer with sequence and numerical features. Must be connected to internet to load the pretrained model for the first time"""
    input:
        X_train = 'data/references/positives_and_negatives/added_features/added_features_training_set_fixed_length_{fixed_length}nt.tsv',
        y_train = 'data/references/positives_and_negatives/added_features/added_features_training_target_{fixed_length}nt.tsv'
    output:
        model = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_2_classes_fold_{fold_num}.pt',
        fold_loss = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_2_classes_fold_{fold_num}_loss_per_epoch.tsv',
        fold_f1_score = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_2_classes_fold_{fold_num}_f1_score_per_epoch.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/transformer_w_features_2_classes/transformer_training_w_features_2_classes_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_training_w_features_2_classes.sh "
        "{params.pretrained_model} {wildcards.fold_num} "
        "{params.random_state} "
        "{input.X_train} {input.y_train} "
        "{output.model} {output.fold_loss} {output.fold_f1_score} &> {log}"


rule test_transformer_2_classes:
    """ Test the performance of the trained transform on the actual test
        set. Don't forget to use the best hyperparams (and specify
        the model architecture before loading the weights/parameters
        learned during training). This transformer predicts only 2 classes (other vs sno (sno|pseudosno))."""
    input:
        X_test = rules.get_three_sets_initial_fixed_length.output.test,
        y_test = rules.get_three_sets_initial_fixed_length.output.test_target,
        best_hyperparams = rules.hypertuning_transformer_2_classes.output.best_hyperparams,
        model = rules.training_transformer_2_classes.output.model 
    output:
        df_metrics_on_test = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_test_metrics_{fixed_length}nt_fold_{fold_num}.tsv',
        test_predictions = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_test_predictions_{fixed_length}nt_fold_{fold_num}.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        python_script = "scripts/python/test_transformer_2_classes.py"
    log:
        "logs/test_transformer_2_classes_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_test_2_classes.sh "
        "{params.pretrained_model} "
        "{input.X_test} {input.y_test} {input.best_hyperparams} "
        "{input.model} {output.df_metrics_on_test} {output.test_predictions} "
        "{params.python_script} &> {log}"

rule test_transformer_2_classes_LR_schedule:
    """ Test the performance of the trained transform on the actual test
        set. Don't forget to use the best hyperparams (and specify
        the model architecture before loading the weights/parameters
        learned during training). This transformer predicts only 2 classes (other vs sno (sno|pseudosno))."""
    input:
        #X_test = rules.get_three_sets_initial_fixed_length.output.test,
        #y_test = rules.get_three_sets_initial_fixed_length.output.test_target,
        #X_test = 'data/references/positives_and_negatives/initial/initial_test_set_fixed_length_194nt.tsv',
        #y_test = 'data/references/positives_and_negatives/initial/initial_test_target_fixed_length_194nt.tsv',
        X_test = 'data/references/positives_and_negatives/data_augmentation/test_set_fixed_length_194nt.tsv',
        y_test = 'data/references/positives_and_negatives/data_augmentation/test_target_fixed_length_194nt.tsv',
        #best_hyperparams = rules.hypertuning_transformer_2_classes.output.best_hyperparams,
        best_hyperparams = 'results/predictions/transformer/194/transformer_2_classes_best_hyperparams.tsv',
        model = rules.training_transformer_2_classes_LR_schedule.output.model
    output:
        df_metrics_on_test = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_LR_schedule_test_metrics_{fixed_length}nt_fold_{fold_num}.tsv',
        test_predictions = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_LR_schedule_test_predictions_{fixed_length}nt_fold_{fold_num}.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        python_script = "scripts/python/test_transformer_2_classes.py"
    log:
        "logs/test_transformer_2_classes_LR_schedule_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_test_2_classes.sh "
        "{params.pretrained_model} "
        "{input.X_test} {input.y_test} {input.best_hyperparams} "
        "{input.model} {output.df_metrics_on_test} {output.test_predictions} "
        "{params.python_script} &> {log}"

rule learning_curve_avg_f1_score_training_transformer:
    """ Create average learning curve (of avg f1-score across 3 classes) 
        across 10 folds on training set."""
    input:
        f1_score_tsv = glob.glob('/home/etienne/Narval/scratch/cd_predictor/workflow/results/predictions/transformer/211/transformer_trained_fold_*f1_score_per_epoch.tsv')
    output:
        learning_curve = 'results/figures/lineplot/transformer/211nt/transformer_training_f1_score_avg_across_fold.svg'
    params:
        num_epoch = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/learning_curve_avg_f1_score_training_transformer.py"

rule learning_curve_avg_f1_score_training_transformer_w_features:
    """ Create average learning curve (of avg f1-score across 3 classes) 
        across 10 folds on training set for transformer trained also with 
        4 numerical features."""
    input:
        f1_score_tsv = glob.glob('/home/etienne/Narval/scratch/cd_predictor/workflow/results/predictions/transformer/211/transformer_trained_w_features_*f1_score_per_epoch.tsv')
    output:
        learning_curve = 'results/figures/lineplot/transformer/211nt/transformer_training_f1_score_avg_across_fold_w_features.svg'
    params:
        num_epoch = 25
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/learning_curve_avg_f1_score_training_transformer.py"

rule learning_curve_avg_f1_score_training_transformer_2_classes:
    """ Create average learning curve (of avg f1-score across 2 classes (other vs sno (sno|pseudosno))) 
        across 10 folds on training set for transformer trained w sequence only."""
    input:
        f1_before_train = glob.glob('/home/etienne/Narval/scratch/cd_predictor/workflow/results/predictions/transformer/190/transformer_2_classes_Before_t*f1_score_per_epoch.tsv'),
        f1_score_tsv = glob.glob('/home/etienne/Narval/scratch/cd_predictor/workflow/results/predictions/transformer/190/transformer_2_classes_*LR*f1_score_per_epoch.tsv')
    output:
        learning_curve = 'results/figures/lineplot/transformer/190nt/transformer_2_classes_training_f1_score_avg_across_fold.svg'
    params:
        num_epoch = 30
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/learning_curve_avg_f1_score_training_transformer.py"

rule learning_curve_avg_f1_score_training_transformer_w_features_2_classes:
    """ Create average learning curve (of avg f1-score across 2 classes (other vs sno (sno|pseudosno)))
        across 10 folds on training set for transformer trained w sequence and 4 intrinsic features."""
    input:
        f1_score_tsv = glob.glob('/home/etienne/Narval/scratch/cd_predictor/workflow/results/predictions/transformer/211/transformer_trained_w_features_2_classes_fold_*f1_score_per_epoch.tsv')
    output:
        learning_curve = 'results/figures/lineplot/transformer/211nt/transformer_w_features_2_classes_training_f1_score_avg_across_fold.svg'
    params:
        num_epoch = 200
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/learning_curve_avg_f1_score_training_transformer.py"

