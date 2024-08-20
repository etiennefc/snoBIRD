rule box_score:
    """ Compute a box score for all examples (score of 0, 
        means all the C/D/C'/D' boxes are perfect, and 
        with each mutation wih regards to the consensus, 
        the score increases)."""
    input:
        positives_fa = rules.tuning_train_test_split_rfam_fixed_length.output.all_positives,
        negatives_fa = rules.get_all_initial_negatives_fixed_length.output
    output:
        positives = 'data/references/box_score/all_positives_box_score_{fixed_length}nt.tsv',
        negatives = 'data/references/box_score/all_negatives_and_pseudosno_box_score_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/box_score.py"

rule structure_stability:
    """ Compute the structure stability of the snoRNA. Its boundaries are 
        defined by 5nt before/after the C/D boxes respectively. The 5 nt 
        is the median distance between a sno boundary and either the C or 
        D box based on all expressed C/D snoRNA distributions (see 
        density_stem_length output)."""
    input:
        positives_df = rules.tuning_train_test_split_rfam_fixed_length.output.all_positives,
        negatives_df = rules.get_all_initial_negatives_fixed_length.output,
        box_score_positives = rules.box_score.output.positives,
        box_score_negatives = rules.box_score.output.negatives
    output:
        positives_fa = 'data/references/structure/all_positives_structure_stability_{fixed_length}nt.fa',
        negatives_fa = 'data/references/structure/all_negatives_and_pseudosno_structure_stability_{fixed_length}nt.fa',
        positives_tsv = 'data/references/structure/all_positives_structure_stability_{fixed_length}nt.tsv',
        negatives_tsv = 'data/references/structure/all_negatives_and_pseudosno_structure_stability_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/structure_stability.py"

rule predicted_length:
    """ Get the length of the snoRNA based on the boxes found 
        in the structure_stability rule."""
    input:
        positives_fa = rules.structure_stability.output.positives_fa,
        negatives_fa = rules.structure_stability.output.negatives_fa
    output:
        positives_df = 'data/references/predicted_length/all_positives_length_{fixed_length}nt.tsv',
        negatives_df = 'data/references/predicted_length/all_negatives_length_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/predicted_length.py"


rule terminal_stem_stability:
    """ Get the terminal stem stability of the snoRNA 
        (15 nt outside of the snoRNA boundaries and 5 
        internal nt, if possible)."""
    input:
        positives_windows = rules.tuning_train_test_split_rfam_fixed_length.output.all_positives,
        negatives_windows = rules.get_all_initial_negatives_fixed_length.output,
        positives_len_fa = rules.structure_stability.output.positives_fa,
        negatives_len_fa = rules.structure_stability.output.negatives_fa
    output:
        positives_fa = 'data/references/terminal_stem/all_positives_terminal_stem_stability_{fixed_length}nt.fa',
        negatives_fa = 'data/references/terminal_stem/all_negatives_and_pseudosno_terminal_stem_stability_{fixed_length}nt.fa',
        positives_tsv = 'data/references/terminal_stem/all_positives_terminal_stem_stability_{fixed_length}nt.tsv',
        negatives_tsv = 'data/references/terminal_stem/all_negatives_and_pseudosno_terminal_stem_stability_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/terminal_stem_stability.py"
