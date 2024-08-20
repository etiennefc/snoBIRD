#rule get_three_sets_initial:
#    """ Concat the positives and negative examples into the 
#        tuning (10%), training (70%) and test (20%) sets. This is 
#        the initial split based on sequence only without any 
#        snoRNA pseudogenes in the negative examples. It is a 
#        stratified split (same proportion of pos/neg examples 
#        across the 3 sets, i.e. ~ 20:1 (negatives:positives))."""
#    input:
#        positives = rules.tuning_train_test_split_rfam.output,
#        negatives = rules.get_all_initial_negatives.output
#    output:
#        tuning = 'data/references/positives_and_negatives/initial/initial_tuning_set.tsv',
#        training = 'data/references/positives_and_negatives/initial/initial_training_set.tsv',
#        test = 'data/references/positives_and_negatives/initial/initial_test_set.tsv'
#    params:
#        random_state = 42,
#        short_name_dict = config['species_short_name']
#    conda:
#        "../envs/python_new.yaml"
#    script:
#        "../scripts/python/get_three_sets_initial.py"

rule get_three_sets_initial_fixed_length:
    """ Concat the positives and negative examples into the 
        tuning (10%), training (70%) and test (20%) sets. This is 
        the initial split based on sequence only. It is a 
        stratified split (same proportion of pos/pseudosno/neg examples 
        across the 3 sets, i.e. ~ 80:1 (negatives:positives)) 
        (25/25/25 intronic/exonic/intergenic, the rest is 1 for 
        shuffled sno, HACA, tRNA, snRNA, premiRNA)."""
    input:
        positives = rules.tuning_train_test_split_rfam_sno_pseudo_fixed_length.output,
        negatives = rules.get_all_initial_negatives_wo_pseudo_fixed_length.output
    output:
        tuning = 'data/references/positives_and_negatives/initial/initial_tuning_set_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/positives_and_negatives/initial/initial_training_set_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/positives_and_negatives/initial/initial_test_set_fixed_length_{fixed_length}nt.tsv',
        tuning_target = 'data/references/positives_and_negatives/initial/initial_tuning_target_fixed_length_{fixed_length}nt.tsv',
        training_target = 'data/references/positives_and_negatives/initial/initial_training_target_fixed_length_{fixed_length}nt.tsv',
        test_target = 'data/references/positives_and_negatives/initial/initial_test_target_fixed_length_{fixed_length}nt.tsv'
    params:
        random_state = 42,
        short_name_dict = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_three_sets_initial_fixed_length.py"

rule get_three_sets_initial_fixed_length_data_aug:
    """ Concat the positives and negative examples into the 
        tuning (10%), training (70%) and test (20%) sets. This is 
        the initial split based on sequence only. It is a 
        stratified split (same proportion of pos/pseudosno/neg examples 
        across the 3 sets, i.e. ~ 80:1 (negatives:positives)) 
        (25/25/25 intronic/exonic/intergenic, the rest is 1 for 
        shuffled sno, HACA, tRNA, snRNA, premiRNA)."""
    input:
        positives = rules.tuning_train_test_split_rfam_sno_pseudo_fixed_length_data_aug.output,
        negatives = rules.get_all_initial_negatives_wo_pseudo_fixed_length_data_aug.output
    output:
        tuning = 'data/references/positives_and_negatives/data_augmentation/tuning_set_1_ratio_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/positives_and_negatives/data_augmentation/training_set_1_ratio_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/positives_and_negatives/data_augmentation/test_set_1_ratio_fixed_length_{fixed_length}nt.tsv',
        tuning_target = 'data/references/positives_and_negatives/data_augmentation/tuning_target_1_ratio_fixed_length_{fixed_length}nt.tsv',
        training_target = 'data/references/positives_and_negatives/data_augmentation/training_target_1_ratio_fixed_length_{fixed_length}nt.tsv',
        test_target = 'data/references/positives_and_negatives/data_augmentation/test_target_1_ratio_fixed_length_{fixed_length}nt.tsv'
    params:
        random_state = 42,
        short_name_dict = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_three_sets_initial_fixed_length_data_aug.py"

rule get_three_sets_initial_fixed_length_data_aug_equal_ratio:
    """ Concat the positives and negative examples into the 
        tuning (10%), training (70%) and test (20%) sets. This is 
        the initial split.."""
    input:
        positives = rules.tuning_train_test_split_rfam_sno_pseudo_fixed_length_data_aug_equal_ratio.output
    output:
        tuning = 'data/references/positives_and_negatives/data_augmentation_equal_ratio/tuning_set_equal_ratio_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/positives_and_negatives/data_augmentation_equal_ratio/training_set_equal_ratio_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/positives_and_negatives/data_augmentation_equal_ratio/test_set_equal_ratio_fixed_length_{fixed_length}nt.tsv',
        tuning_target = 'data/references/positives_and_negatives/data_augmentation_equal_ratio/tuning_target_equal_ratio_fixed_length_{fixed_length}nt.tsv',
        training_target = 'data/references/positives_and_negatives/data_augmentation_equal_ratio/training_target_equal_ratio_fixed_length_{fixed_length}nt.tsv',
        test_target = 'data/references/positives_and_negatives/data_augmentation_equal_ratio/test_target_equal_ratio_fixed_length_{fixed_length}nt.tsv'
    params:
        random_state = 42,
        short_name_dict = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_three_sets_initial_fixed_length_data_aug_equal_ratio.py"

rule get_three_sets_initial_fixed_length_data_aug_equal_ratio_combined_tune_train:
    """ Concat the positives and negative examples into the 
        tuning (10%), training (70%) and test (20%) sets. This is 
        the initial split.."""
    input:
        positives = rules.tuning_train_test_split_rfam_sno_pseudo_fixed_length_data_aug_equal_ratio.output
    output:
        training = 'data/references/positives_and_negatives/data_augmentation_equal_ratio_combined_tune_train/training_set_equal_ratio_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/positives_and_negatives/data_augmentation_equal_ratio_combined_tune_train/test_set_equal_ratio_fixed_length_{fixed_length}nt.tsv',
        training_target = 'data/references/positives_and_negatives/data_augmentation_equal_ratio_combined_tune_train/training_target_equal_ratio_fixed_length_{fixed_length}nt.tsv',
        test_target = 'data/references/positives_and_negatives/data_augmentation_equal_ratio_combined_tune_train/test_target_equal_ratio_fixed_length_{fixed_length}nt.tsv'
    params:
        random_state = 42,
        short_name_dict = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_three_sets_initial_fixed_length_data_aug_equal_ratio_combined_tune_train.py"

rule get_three_sets_added_features_fixed_length:
    """ Concat the positives and negative examples into the 
        tuning (10%), training (70%) and test (20%) sets. This is 
        the latest split based on sequence, box_score, structure and 
        terminal_stem stabilities. It is a stratified split 
        (same proportion of pos/pseudosno/neg examples 
        across the 3 sets, i.e. ~ 20:1 (negatives:positives)). 
        We actually do the split based on the sequence features only, 
        then merge the intrinsinc features to the resulting dfs."""
    input:
        positives = rules.tuning_train_test_split_rfam_fixed_length.output,
        negatives = rules.get_all_initial_negatives_fixed_length.output,
        box_score = rules.box_score.output,
        structure_stability = rules.structure_stability.output,
        terminal_stem_stability = rules.terminal_stem_stability.output,
        length = rules.predicted_length.output
    output:
        tuning = 'data/references/positives_and_negatives/added_features/added_features_tuning_set_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/positives_and_negatives/added_features/added_features_training_set_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/positives_and_negatives/added_features/added_features_test_set_fixed_length_{fixed_length}nt.tsv'
    params:
        random_state = 42,
        short_name_dict = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_three_sets_added_features_fixed_length.py"    
