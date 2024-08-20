rule length_selection_positives_density:
    """ Create a density plot showing the length 
        distribution of positive examples. Also 
        show different thresholds based on 
        percentile/quartiles to choose a cutoff. 
        The chosen percentile is 95% (< 183 nt)"""
    input:
        positives = rules.get_sno_sequences.output.df
    output:
        density = 'results/figures/density/length_selection_positives.svg'
    params:
        percent_colors = config["colors"]['percent_colors']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/length_selection_positives_density.py"

rule roc_curve_cd_predictors:
    """ Create a ROC curve of the existing CD predictors.*******NOT COMPLETED****"""
    input:
        snoreport = rules.filter_snoreport_predictions.output.predictions_tsv,
        snoscan = rules.filter_snoscan_predictions.output.predictions_tsv,
        infernal_rfam = rules.filter_rfam_infernal_predictions.output.predictions_tsv
    output:
        roc = 'results/figures/roc/existing_cd_predictors_{fixed_length}.svg'
    params:
        colors = config['colors']['predictors']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/roc_curve_cd_predictors.py"
        

#rule pr_curve_cd_predictor:

rule metrics_lineplot_predictors:
    """ Compute the precision, recall, accuracy and f1-score 
        on the test set (relative to other vs expressed_CD_snoRNA) 
        and also the precision/recall on the snoRNA_pseudogene class 
        for the existing cd predictors and simple models (knn, gbm, 
        logreg, svc and rf) and display it as a dot plot."""
    input:
        #snoreport = rules.filter_snoreport_predictions.output.predictions_tsv,
        #snoscan = rules.filter_snoscan_predictions.output.predictions_tsv,
        #infernal_rfam = rules.filter_rfam_infernal_predictions.output.predictions_tsv,
        #snoBIRD = 'results/predictions/snoBIRD/transformer/194/3e-5_3e-6_32_4_data_aug_1_ratio/transformer_2_classes_LR_schedule_test_predictions_194nt_fold_8.tsv'
    output:
        dotplot_all = 'results/figures/lineplot/metrics_existing_cd_predictors_{fixed_length}_w_pseudogenes.svg',
        dotplot_expressed_cd = 'results/figures/lineplot/metrics_existing_cd_predictors_{fixed_length}_expressed_CD_only.svg',
        dotplot_pseudo = 'results/figures/lineplot/metrics_existing_cd_predictors_{fixed_length}_pseudogenes_only.svg'
    params:
        predictors_colors = config['colors']['predictors'],
        simple_models_colors = config['colors']['simple_models']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/metrics_lineplot_predictors.py"

rule s_pombe_prediction_time_barplot:
    """ Make a bar plot to show the prediction time of the 
        predictors vs snoBIRD on the entirety of the S. pombe genome."""
    output:
        barplot = 'results/figures/barplot/s_pombe_prediction_time.svg'
    params:
        predictors_colors = config['colors']['predictors']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/s_pombe_prediction_time_barplot.py"

rule FP_FN_initial_analyses_pie:
    """ Create pie charts showing the proportion of species
        or type of negatives predicted as false positives/negatives 
        (FP/FN) by existing CD_predictors and by the chosen trained GRU NN. """
    input:
        snoreport = rules.filter_snoreport_predictions.output.predictions_tsv,
        snoscan = rules.filter_snoscan_predictions.output.predictions_tsv,
        infernal_rfam = rules.filter_rfam_infernal_predictions.output.predictions_tsv,
        gru_nn = rules.test_gru.output.test_predictions,
        features = rules.get_three_sets_added_features_fixed_length.output
    output:
        pie_species = 'results/figures/pie/{error}_per_species_existing_cd_predictors_{fixed_length}.svg',
        pie_neg_type = 'results/figures/pie/{error}_per_negative_type_existing_cd_predictors_{fixed_length}.svg'
    params:
        species_colors = config['colors']['species'],
        biotype_colors = config['colors']['biotypes']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/FP_FN_initial_analyses_pie.py"

rule donut_positives_negatives:
    """ Create donut charts showing the proportion of species
        for the positive examples and the proportion of 
        negative types for the negatives."""
    input:
        #negatives = rules.get_all_initial_negatives_fixed_length.output,
        #positives = rules.tuning_train_test_split_rfam_fixed_length.output.all_positives
    output:
        pie_species = 'results/figures/pie/positives_per_species_{fixed_length}.svg',
        pie_neg_type = 'results/figures/pie/negatives_per_biotype_{fixed_length}.svg'
    params:
        species_colors = config['colors']['species'],
        biotype_colors = config['colors']['biotypes'],
        species_name = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/donut_positives_negatives.py"

rule bar_confusion_value_per_species_test:
    """ Create a stacked bar chart showing the proportion of test set examples 
        predicted as FP, FN, TP, TN per species, with the total number of 
        examples in the set for the given species above the bars."""
    input:
        snoreport = rules.filter_snoreport_predictions.output.predictions_tsv,
        snoscan = rules.filter_snoscan_predictions.output.predictions_tsv,
        infernal_rfam = rules.filter_rfam_infernal_predictions.output.predictions_tsv
    output:
        bar_all = 'results/figures/barplot/confusion_values_per_species_{cd_predictors}_{fixed_length}.svg',
        bar_FN_FP = 'results/figures/barplot/FN_FP_per_species_{cd_predictors}_{fixed_length}.svg'
    params:
        species_colors = config['colors']['species'],
        conf_value_colors = config['colors']['confusion_value']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/bar_confusion_value_per_species_test.py" 


rule violin_sno_length_per_species:
    """ Create a violin plot to show the length distribution of expressed C/D 
        and snoRNA pseudogenes per species with a threshold of fixed_length."""
    input:
        expressed_cd = rules.get_sno_sequences_fixed_length.output.df,
        human_pseudo = rules.get_human_snoRNA_pseudogenes.output.pseudogenes,
        mouse_pseudo = rules.get_mouse_snoRNA_pseudogenes.output.pseudogenes,
        drosophila_pseudo = rules.get_drosophila_snoRNA_pseudogenes.output.pseudogenes 
    output:
        violin = 'results/figures/violinplot/sno_pseudo_length_per_species_{fixed_length}.svg'
    params:
        target_colors = config['colors']['target'],
        short_sp_names = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/violin_sno_length_per_species.py" 


rule density_confusion_value_per_predictor:
    """ Create a density plot of the 4 intrinsic features 
        per confusion value (TP, TN, FP and FN) for each predictor 
        (existing tools and GRU_nn)."""
    input:
        snoreport = rules.filter_snoreport_predictions.output.predictions_tsv,
        snoscan = rules.filter_snoscan_predictions.output.predictions_tsv,
        infernal_rfam = rules.filter_rfam_infernal_predictions.output.predictions_tsv,
        gru_nn = rules.test_gru.output.test_predictions,
        features = rules.get_three_sets_added_features_fixed_length.output.test
    output:
        density = 'results/figures/density/{intrinsic_feature}_confusion_value_all_models_{fixed_length}.svg'
    params:
        conf_value_colors = config['colors']['confusion_value']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/density_confusion_value_per_predictor.py" 


rule upset_confusion_value_models:
    """ Create a upset plot for each confusion value to see if the existing predictors 
        and GRU NN display the same prediction errors on the test set."""
    input:
        snoreport = rules.filter_snoreport_predictions.output.predictions_tsv,
        snoscan = rules.filter_snoscan_predictions.output.predictions_tsv,
        infernal_rfam = rules.filter_rfam_infernal_predictions.output.predictions_tsv,
        gru_nn = rules.test_gru.output.test_predictions
    output:
        upset_TP = 'results/figures/upset/upset_TP_all_models_{fixed_length}nt_species.svg',
        upset_FP = 'results/figures/upset/upset_FP_all_models_{fixed_length}nt_species.svg',
        upset_FN = 'results/figures/upset/upset_FN_all_models_{fixed_length}nt_species.svg',
        upset_TN = 'results/figures/upset/upset_TN_all_models_{fixed_length}nt_species.svg',
        upset_TN_biotype = 'results/figures/upset/upset_TN_all_models_{fixed_length}nt_biotype.svg',
        upset_FP_biotype = 'results/figures/upset/upset_FP_all_models_{fixed_length}nt_biotype.svg'
    params:
        species_colors = config['colors']['species'],
        biotype_colors = config['colors']['biotypes']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/upset_confusion_value_models.py"



rule density_box_score:
    """ Create a density plot of the box score distribution of positives, 
        pseudosno and negatives. Create also another density plot of 
        negatives only with a hue of gene_biotype."""
    input:
        positives = rules.box_score.output.positives,
        negatives = rules.box_score.output.negatives,
        biotype_df = rules.get_all_initial_negatives_fixed_length.output
    output:
        density_all = 'results/figures/density/box_score_positives_negatives_{fixed_length}nt.svg',
        density_negatives = 'results/figures/density/box_score_negatives_gene_biotype_{fixed_length}nt.svg'
    params:
        biotype_colors = config['colors']['biotypes'],
        target_colors = config['colors']['target']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/density_box_score.py"

rule density_stem_length:
    """ Predict the C and D box positions in expressed C/D 
        and find the number of nt up/downstream until the 
        end of the annotated snoRNA. Return the distribution
        of length for both positions. This will serve to define 
        the terminal stem length threshold when computing the 
        terminal stem stability of a given window."""
    input:
        positives_df = rules.get_sno_sequences_fixed_length.output.df
    output:
        density = 'results/figures/density/terminal_stem_length_distributions_positives_{fixed_length}nt.svg'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/density_stem_length.py"

rule density_structure_stability_length:
    """ Create a density plot of the structure stability distribution 
        of positives, pseudosno and negatives. Create also another 
        density plot of negatives only with a hue of gene_biotype. 
        Do also this for normalized MFE (by length) and for the
        snoRNA length only """
    input:
        positives = rules.structure_stability.output.positives_tsv,
        negatives = rules.structure_stability.output.negatives_tsv,
        positives_fa = rules.structure_stability.output.positives_fa,
        negatives_fa = rules.structure_stability.output.negatives_fa,
        biotype_df = rules.get_all_initial_negatives_fixed_length.output
    output:
        density_all = 'results/figures/density/structure_stability_positives_negatives_{fixed_length}nt.svg',
        density_negatives = 'results/figures/density/structure_stability_negatives_gene_biotype_{fixed_length}nt.svg',
        density_all_normalized = 'results/figures/density/normalized_structure_stability_positives_negatives_{fixed_length}nt.svg',
        density_negatives_normalized = 'results/figures/density/normalized_structure_stability_negatives_gene_biotype_{fixed_length}nt.svg',
        density_all_length = 'results/figures/density/length_positives_negatives_{fixed_length}nt.svg',
        density_negatives_length = 'results/figures/density/length_negatives_gene_biotype_{fixed_length}nt.svg'
    params:
        biotype_colors = config['colors']['biotypes'],
        target_colors = config['colors']['target']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/density_structure_stability_length.py"

rule density_terminal_stem_stability:
    """ Create a density plot of the terminal stem stability distribution 
        of positives, pseudosno and negatives. Create also another 
        density plot of negatives only with a hue of gene_biotype."""
    input:
        positives = rules.terminal_stem_stability.output.positives_tsv,
        negatives = rules.terminal_stem_stability.output.negatives_tsv,
        biotype_df = rules.get_all_initial_negatives_fixed_length.output
    output:
        density_all = 'results/figures/density/terminal_stem_stability_positives_negatives_{fixed_length}nt.svg',
        density_negatives = 'results/figures/density/terminal_stem_stability_negatives_gene_biotype_{fixed_length}nt.svg'
    params:
        biotype_colors = config['colors']['biotypes'],
        target_colors = config['colors']['target']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/density_terminal_stem_stability.py"

rule learning_curve_avg_f1_score_training_added_features:
    """ Create average learning curve (of avg f1-score across 3 classes) across 10 folds on training set. 
        Input features = sequence + 4 added features."""
    input:
        f1_score_tsv = rules.training_gru_added_features_simplified2.output.all_fold_epochs_df
    output:
        learning_curve = 'results/figures/lineplot/gru/{fixed_length}nt/added_features/gru_training_f1_score_simplified2_{fixed_length}nt_avg_across_fold.svg'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/learning_curve_avg_f1_score_training_added_features.py"

rule learning_curve_avg_f1_score_training_added_features_half_normalized:
    """ Create average learning curve (of avg f1-score across 3 classes) across 10 folds on training set.
        Input features = sequence + 4 added features."""
    input:
        f1_score_tsv = rules.training_gru_added_features_half_normalized_simplified.output.all_fold_epochs_df
    output:
        learning_curve = 'results/figures/lineplot/gru/{fixed_length}nt/added_features_half_normalized/gru_training_f1_score_simplified_{fixed_length}nt_avg_across_fold.svg'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/learning_curve_avg_f1_score_training_added_features.py"

rule bar_cd_expression_species:
    """ Create a stacked bar chart per species with TGIRT-Seq dataset to show 
        the proportion of expressed CD vs pseudogenes."""
    input:
        fake_dep = expand(rules.get_expressed_snoRNAs_location.output, 
                    species=['homo_sapiens', 'saccharomyces_cerevisiae', 'mus_musculus']),
        fake_dep2 = expand(rules.get_D_melanogaster_expressed_snoRNAs_location.output, species='drosophila_melanogaster')
    output:
        bar = 'results/figures/barplot/cd_expression_per_species.svg'
    params:
        colors = config['colors']['target'],
        tgirt_dir = 'data/references/tgirt_seq_output/',
        species_short_name = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/bar_cd_expression_species.py"













