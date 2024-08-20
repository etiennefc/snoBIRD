#rule onehot_encode_normalize_initial:
#    """ One-hot encode the sequence of examples in the tuning, 
#        training and test sets. Then normalize/standardize 
#        SEPARATELY the datasets (ensuring no data leakage 
#        (because of the mean/stdev that would be computed on 
#        the whole datasets otherwise). Transform also the target 
#        column into either 1 (expressed_CD_snoRNA) or 0 (other)."""
#    input:
#        tuning = rules.get_three_sets_initial.output.tuning,
#        training = rules.get_three_sets_initial.output.training,
#        test = rules.get_three_sets_initial.output.test
#    output:
#        encoded_tuning = 'data/references/positives_and_negatives/initial/initial_one_hot_encoded_tuning_set.tsv',
#        encoded_training = 'data/references/positives_and_negatives/initial/initial_one_hot_encoded_training_set.tsv',
#        encoded_test = 'data/references/positives_and_negatives/initial/initial_one_hot_encoded_test_set.tsv',
#        normalized_tuning = 'data/references/positives_and_negatives/initial/initial_encoded_scaled_tuning_set.tsv',
#        normalized_training = 'data/references/positives_and_negatives/initial/initial_encoded_scaled_training_set.tsv',
#        normalized_test = 'data/references/positives_and_negatives/initial/initial_encoded_scaled_test_set.tsv',
#        target_tuning = 'data/references/positives_and_negatives/initial/initial_tuning_target.tsv',
#        target_training = 'data/references/positives_and_negatives/initial/initial_training_target.tsv',
#        target_test = 'data/references/positives_and_negatives/initial/initial_test_target.tsv'
#    conda:
#        "../envs/python_new.yaml"
#    script:
#        "../scripts/python/onehot_encode_normalize_initial.py"

rule onehot_encode_normalize_initial_fixed_length:
    """ One-hot encode the sequence of examples in the tuning, 
        training and test sets. Then normalize/standardize 
        SEPARATELY the datasets (ensuring no data leakage 
        (because of the mean/stdev that would be computed on 
        the whole datasets otherwise). Transform also the target 
        column into either 2 (expressed_CD_snoRNA), 1 (CD_snoRNA_pseudogene) or 0 (other)."""
    input:
        tuning = rules.get_three_sets_initial_fixed_length.output.tuning,
        training = rules.get_three_sets_initial_fixed_length.output.training,
        test = rules.get_three_sets_initial_fixed_length.output.test
    output:
        encoded_tuning = 'data/references/positives_and_negatives/initial/initial_one_hot_encoded_tuning_set_{fixed_length}nt.tsv',
        encoded_training = 'data/references/positives_and_negatives/initial/initial_one_hot_encoded_training_set_{fixed_length}nt.tsv',
        encoded_test = 'data/references/positives_and_negatives/initial/initial_one_hot_encoded_test_set_{fixed_length}nt.tsv',
        normalized_tuning = 'data/references/positives_and_negatives/initial/initial_encoded_scaled_tuning_set_{fixed_length}nt.tsv',
        normalized_training = 'data/references/positives_and_negatives/initial/initial_encoded_scaled_training_set_{fixed_length}nt.tsv',
        normalized_test = 'data/references/positives_and_negatives/initial/initial_encoded_scaled_test_set_{fixed_length}nt.tsv',
        target_tuning = 'data/references/positives_and_negatives/initial/initial_tuning_target_{fixed_length}nt.tsv',
        target_training = 'data/references/positives_and_negatives/initial/initial_training_target_{fixed_length}nt.tsv',
        target_test = 'data/references/positives_and_negatives/initial/initial_test_target_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/onehot_encode_normalize_initial_fixed_length.py"

rule onehot_encode_normalize_added_features_fixed_length:
    """ One-hot encode the sequence of examples in the tuning, 
        training and test sets. Then normalize/standardize 
        SEPARATELY the datasets (ensuring no data leakage 
        (because of the mean/stdev that would be computed on 
        the whole datasets otherwise), including the added features. 
        Transform also the target column into either 2 
        (expressed_CD_snoRNA), 1 (CD_snoRNA_pseudogene) or 0 (other)."""
    input:
        tuning = rules.get_three_sets_added_features_fixed_length.output.tuning,
        training = rules.get_three_sets_added_features_fixed_length.output.training,
        test = rules.get_three_sets_added_features_fixed_length.output.test
    output:
        encoded_tuning = 'data/references/positives_and_negatives/added_features/added_features_one_hot_encoded_tuning_set_{fixed_length}nt.tsv',
        encoded_training = 'data/references/positives_and_negatives/added_features/added_features_one_hot_encoded_training_set_{fixed_length}nt.tsv',
        encoded_test = 'data/references/positives_and_negatives/added_features/added_features_one_hot_encoded_test_set_{fixed_length}nt.tsv',
        normalized_tuning = 'data/references/positives_and_negatives/added_features/added_features_encoded_scaled_tuning_set_{fixed_length}nt.tsv',
        normalized_training = 'data/references/positives_and_negatives/added_features/added_features_encoded_scaled_training_set_{fixed_length}nt.tsv',
        normalized_test = 'data/references/positives_and_negatives/added_features/added_features_encoded_scaled_test_set_{fixed_length}nt.tsv',
        target_tuning = 'data/references/positives_and_negatives/added_features/added_features_tuning_target_{fixed_length}nt.tsv',
        target_training = 'data/references/positives_and_negatives/added_features/added_features_training_target_{fixed_length}nt.tsv',
        target_test = 'data/references/positives_and_negatives/added_features/added_features_test_target_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/onehot_encode_normalize_initial_fixed_length.py"

rule onehot_encode_normalize_added_features_half_normalized_fixed_length:
    """ One-hot encode the sequence of examples in the tuning, 
        training and test sets. Then normalize/standardize 
        SEPARATELY the datasets (ensuring no data leakage 
        (because of the mean/stdev that would be computed on 
        the whole datasets otherwise) just for the the added 
        features (NOT the one-hot encoded sequence features). 
        Transform also the target column into either 2 
        (expressed_CD_snoRNA), 1 (CD_snoRNA_pseudogene) or 0 (other)."""
    input:
        tuning = rules.get_three_sets_added_features_fixed_length.output.tuning,
        training = rules.get_three_sets_added_features_fixed_length.output.training,
        test = rules.get_three_sets_added_features_fixed_length.output.test
    output:
        encoded_tuning = 'data/references/positives_and_negatives/added_features_half_normalized/added_features_one_hot_encoded_tuning_set_{fixed_length}nt.tsv',
        encoded_training = 'data/references/positives_and_negatives/added_features_half_normalized/added_features_one_hot_encoded_training_set_{fixed_length}nt.tsv',
        encoded_test = 'data/references/positives_and_negatives/added_features_half_normalized/added_features_one_hot_encoded_test_set_{fixed_length}nt.tsv',
        target_tuning = 'data/references/positives_and_negatives/added_features_half_normalized/added_features_tuning_target_{fixed_length}nt.tsv',
        target_training = 'data/references/positives_and_negatives/added_features_half_normalized/added_features_training_target_{fixed_length}nt.tsv',
        target_test = 'data/references/positives_and_negatives/added_features_half_normalized/added_features_test_target_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/onehot_encode_half_normalized_fixed_length.py"