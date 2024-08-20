rule compare_snoreport_w_snoBIRD:
    """ Test the performance of snoreport2 on C/D in one chr in C. albicans 
        to compare its speed with snoBIRD's. **WARNING**: Snoreport2 must 
        be installed manually on your computer before 
        using this rule and you might need to change the
        path in the command line below to access snoreport2."""
    input:
        test_fa = "data/references/genome_fa/candida_albicans/Ca22chrM_C_albicans_SC5314.fa"
    output:
        predictions = 'results/predictions/snoBIRD_comparison/snoreport2_predicted_cd.txt',
        time_predictions = 'results/predictions/snoBIRD_comparison/snoreport2_time_predicted_cd.txt'
    shell:
        """/usr/bin/time -f "%E" -o {output.time_predictions} ~/snoReport_2/snoreport_2 -i {input.test_fa} """
        """-CD """
        """--positives """
        """-o {output.predictions}"""

rule test_snoreport:
    """ Test the performance of snoreport2 on the test 
        set of fixed C/D length (by first converting the 
        test set to fasta). **WARNING**: Snoreport2 must 
        be installed manually on your computer before 
        using this rule and you might need to change the
        path in the command line below to access snoreport2."""
    input:
        test_set = rules.get_three_sets_initial_fixed_length_data_aug.output.test
    output:
        predictions = 'results/predictions/snoreport2/fixed_length_{fixed_length}nt/predicted_cd.fa'
    shell:
        """awk 'NR>1 {{print ">"$1" "$2" "$3"\\n"$4}}' {input.test_set} > """
        """temp_snoreport_{wildcards.fixed_length}.fa && """
        """~/snoReport_2/snoreport_2 -i temp_snoreport_{wildcards.fixed_length}.fa """
        """-CD """
        """--positives """
        """-o {output.predictions} && """
        """rm temp_snoreport_{wildcards.fixed_length}.fa"""

rule filter_snoreport_predictions:
    """ Filter snoreport2 predictions to remove duplicates 
        and return a clear target prediction for each example."""
    input:
        predictions_fa = rules.test_snoreport.output.predictions,
        test_set = rules.get_three_sets_initial_fixed_length_data_aug.output.test
    output:
        predictions_tsv = 'results/predictions/snoreport2/fixed_length_{fixed_length}nt/test_predictions.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_snoreport_predictions.py"

rule test_snoreport_pseudosno:
    """ Test the performance of snoreport2 on snoRNA pseudogenes 
        (i.e. not expressed human snoRNAs) (by first converting the 
        table to fasta). **WARNING**: Snoreport2 must 
        be installed manually on your computer before 
        using this rule and you might need to change the
        path in the command line below to access snoreport2."""
    input:
        pseudosno = rules.get_human_snoRNA_pseudogenes.output.pseudogenes
    output:
        predictions = 'results/predictions/snoreport2/fixed_length_{fixed_length}nt/predicted_snoRNA_pseudogenes.fa'
    shell:
        """awk 'NR>1 {{print ">"$1" "$3" "$4"\\n"$11}}' {input.pseudosno} > """
        """temp_snoreport_pseudosno_{wildcards.fixed_length}.fa && """
        """~/snoReport_2/snoreport_2 -i temp_snoreport_pseudosno_{wildcards.fixed_length}.fa """
        """-CD """
        """--positives """
        """-o {output.predictions} && """
        """rm temp_snoreport_pseudosno_{wildcards.fixed_length}.fa"""

rule compare_snoscan_w_snoBIRD:
    """ Predict with snoscan the presence of 
        C/D in one chr in C. albicans to compare its speed with snoBIRD's."""
    input:
        target_rDNA = "data/references/genome_fa/candida_albicans/rDNA.fa",
        test_fa = "data/references/genome_fa/candida_albicans/Ca22chrM_C_albicans_SC5314.fa"
    output:
        predictions = 'results/predictions/snoBIRD_comparison/snoscan_predicted_cd.txt',
        time_predictions = 'results/predictions/snoBIRD_comparison/snoscan_time_predicted_cd.txt'
    conda: 
        "../envs/snoscan.yaml"
    script:
        "../scripts/python/compare_snoscan_w_snoBIRD.py"

rule test_snoscan:
    """ Predict with snoscan the presence of 
        expressed C/D in the test set."""
    input:
        target_rDNA = rules.download_rDNA.output.rDNA_fa,
        test_set = rules.get_three_sets_initial_fixed_length_data_aug.output.test
    output:
        predictions = 'results/predictions/snoscan/fixed_length_{fixed_length}nt/predicted_cd.txt'
    conda: 
        "../envs/snoscan.yaml"
    script:
        "../scripts/python/test_sno_scan.py"

rule filter_snoscan_predictions:
    """ Filter snoscan predictions to remove duplicates 
        and return a clear target prediction for each example."""
    input:
        predictions_fa = rules.test_snoscan.output.predictions,
        test_set = rules.get_three_sets_initial_fixed_length_data_aug.output.test
    output:
        predictions_tsv = 'results/predictions/snoscan/fixed_length_{fixed_length}nt/test_predictions.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_snoscan_predictions.py"

rule test_snoscan_pseudosno:
    """ Predict with snoscan the presence of 
        expressed C/D in the snoRNA pseudogenes."""
    input:
        target_rDNA = rules.download_rDNA.output.rDNA_fa,
        pseudosno = rules.get_human_snoRNA_pseudogenes.output.pseudogenes
    output:
        predictions = 'results/predictions/snoscan/fixed_length_{fixed_length}nt/predicted_snoRNA_pseudogenes.txt'
    conda: 
        "../envs/snoscan.yaml"
    script:
        "../scripts/python/test_sno_scan_pseudosno.py"

rule compare_rfam_infernal_w_snoBIRD:
    """ Predict with infernal the presence of 
        C/D in one chr in C. albicans to compare its speed with snoBIRD's."""
    input:
        rfam_cm = rules.download_rfam_covariance_models.output.rfam_cm,
        test_fa = "data/references/genome_fa/candida_albicans/Ca22chrM_C_albicans_SC5314.fa"
    output:
        predictions_tblout = 'results/predictions/snoBIRD_comparison/infernal_rfam_predicted_cd.tblout',
        predictions_alignments = 'results/predictions/snoBIRD_comparison/infernal_rfam_predicted_cd.txt',
        time_predictions = 'results/predictions/snoBIRD_comparison/infernal_rfam_time_predicted_cd.txt'
    conda: 
        "../envs/infernal.yaml"
    shell:
        """/usr/bin/time -f "%E" -o {output.time_predictions} """
        """cmpress -F {input.rfam_cm} && """
        """/usr/bin/time -f "%E" -a -o {output.time_predictions} cmscan --cut_ga --rfam --nohmmonly -o {output.predictions_alignments} """
        """--tblout {output.predictions_tblout} {input.rfam_cm} {input.test_fa} """

rule test_rfam_infernal:
    """ Use Infernal and Rfam covariance 
        models to predict if testset examples 
        are part of a C/D Rfam family."""
    input:
        test_set = rules.get_three_sets_initial_fixed_length_data_aug.output.test,
        rfam_cm = rules.download_rfam_covariance_models.output.rfam_cm
    output:
        infernal_tblout = 'results/predictions/infernal_rfam/fixed_length_{fixed_length}nt/predicted_cd.tblout',
        infernal_alignments = 'results/predictions/infernal_rfam/fixed_length_{fixed_length}nt/predicted_cd.txt'
    conda:
        "../envs/infernal.yaml"
    shell:
        """awk 'NR>1 {{print ">"$1" "$2" "$3"\\n"$4}}' {input.test_set} > infernal_temp.fa && """
        """cmpress -F {input.rfam_cm} && """
        """cmscan --cut_ga --rfam --nohmmonly -o {output.infernal_alignments} """
        """--tblout {output.infernal_tblout} {input.rfam_cm} infernal_temp.fa && """
        """rm infernal_temp.fa"""

rule filter_rfam_infernal_predictions:
    """ Filter infernal/rfam predictions to remove duplicates 
        and return a clear target prediction for each example."""
    input:
        predictions_table = rules.test_rfam_infernal.output.infernal_tblout,
        test_set = rules.get_three_sets_initial_fixed_length_data_aug.output.test
    output:
        predictions_tsv = 'results/predictions/infernal_rfam/fixed_length_{fixed_length}nt/test_predictions.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_rfam_infernal_predictions.py"

rule test_rfam_infernal_pseudosno:
    """ Use Infernal and Rfam covariance 
        models to predict if snoRNA pseudogenes 
        are predicted as part of a C/D Rfam family."""
    input:
        pseudosno = rules.get_human_snoRNA_pseudogenes.output.pseudogenes,
        rfam_cm = rules.download_rfam_covariance_models.output.rfam_cm
    output:
        infernal_tblout = 'results/predictions/infernal_rfam/fixed_length_{fixed_length}nt/predicted_snoRNA_pseudogenes.tblout',
        infernal_alignments = 'results/predictions/infernal_rfam/fixed_length_{fixed_length}nt/predicted_snoRNA_pseudogenes.txt'
    conda:
        "../envs/infernal.yaml"
    shell:
        """awk 'NR>1 {{print ">"$1" "$3" "$4"\\n"$11}}' {input.pseudosno} > infernal_temp_pseudosno.fa && """
        """cmpress -F {input.rfam_cm} && """
        """cmscan --cut_ga --rfam --nohmmonly -o {output.infernal_alignments} """
        """--tblout {output.infernal_tblout} {input.rfam_cm} infernal_temp_pseudosno.fa && """
        """rm infernal_temp_pseudosno.fa"""

rule filter_cd_predictors_pseudosno:
    """ Filter the results of snoreport2, snoscan and 
        infernal on snoRNA pseudogenes."""
    input:
        snoreport = rules.test_snoreport_pseudosno.output.predictions,
        snoscan = rules.test_snoscan_pseudosno.output.predictions,
        infernal_rfam = rules.test_rfam_infernal_pseudosno.output.infernal_tblout,
        all_pseudosno = rules.get_human_snoRNA_pseudogenes.output.pseudogenes
    output:
        pseudosno_preds = 'results/predictions/snoRNA_pseudogenes/fixed_length_{fixed_length}nt/snoreport_snoscan_infernal_snoRNA_pseudogenes_predictions.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_cd_predictors_pseudosno.py"

rule confusion_matrix_cd_predictors:
    """ Generate a confusion matrix for each existing cd predictors based on 
        the expressed_CD vs other (not considering the snoRNA pseudogene 
        class here)."""
    input:
        snoreport = rules.filter_snoreport_predictions.output.predictions_tsv,
        snoscan = rules.filter_snoscan_predictions.output.predictions_tsv,
        infernal_rfam = rules.filter_rfam_infernal_predictions.output.predictions_tsv
    output:
        matrix_snoreport = 'results/predictions/snoreport2/fixed_length_{fixed_length}nt/confusion_matrix.tsv',
        matrix_snoscan = 'results/predictions/snoscan/fixed_length_{fixed_length}nt/confusion_matrix.tsv',
        matrix_infernal_rfam = 'results/predictions/infernal_rfam/fixed_length_{fixed_length}nt/confusion_matrix.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/confusion_matrix_cd_predictors.py"

rule test_rfam_infernal_cerevisiae:
    """ Use Infernal and Rfam covariance 
        models to predict on the genome of S. cerevisiae to see if some sequences 
        are part of a C/D Rfam family."""
    input:
        test_genome = 'data/references/genome_fa/saccharomyces_cerevisiae_genome.fa',
        rfam_cm = rules.download_rfam_covariance_models.output.rfam_cm
    output:
        infernal_tblout = 'results/predictions/infernal_rfam/S_cerevisiae/predicted_cd.tblout',
        infernal_alignments = 'results/predictions/infernal_rfam/S_cerevisiae/predicted_cd.txt'
    conda:
        "../envs/infernal.yaml"
    shell:
        """cmpress -F {input.rfam_cm} && """
        """cmscan --cut_ga --rfam --nohmmonly -o {output.infernal_alignments} """
        """--tblout {output.infernal_tblout} {input.rfam_cm} {input.test_genome}"""

rule overlap_rfam_infernal_cerevisiae_CD:
    """ Filter rfam_infernal predictions on S. cerevisiae genome to see the overlap between 
        the predictions and the existing expressed C/D. Return the precision and recall."""
    input:
        infernal_rfam = rules.test_rfam_infernal_cerevisiae.output.infernal_tblout,
        cd_sno = 'data/references/positives/cd_rfam_filtered_all_fixed_length_190nt.tsv',
        haca_yeast = 'data/references/sno_type_df/saccharomyces_cerevisiae_snotype_umass.tsv'
    output:
        df = 'results/predictions/infernal_rfam/S_cerevisiae/precision_recall_cerevisiae.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/overlap_rfam_infernal_cerevisiae_CD.py"

rule test_snoscan_cerevisiae:
    """ Predict with snoscan the presence of 
        C/D in S. cerevisiae genome."""
    input:
        target_rDNA = rules.download_rDNA.output.rDNA_fa,
        test_genome = 'data/references/genome_fa/saccharomyces_cerevisiae_genome.fa'
    output:
        predictions = 'results/predictions/snoscan/S_cerevisiae/predicted_cd.txt'
    conda: 
        "../envs/snoscan.yaml"
    shell:
        "snoscan {input.target_rDNA} {input.test_genome} -o {output.predictions}"

rule overlap_snoscan_cerevisiae_CD:
    """ Filter snoscan predictions on S. cerevisiae genome to see the overlap between 
        the predictions and the existing expressed C/D. Return the precision and recall."""
    input:
        snoscan = rules.test_snoscan_cerevisiae.output.predictions,
        cd_sno = 'data/references/positives/cd_rfam_filtered_all_fixed_length_190nt.tsv',
    output:
        df = 'results/predictions/snoscan/S_cerevisiae/precision_recall_cerevisiae.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/overlap_snoscan_cerevisiae_CD.py"

rule test_snoreport_cerevisiae:
    """ Test the performance of snoreport2 on S cerevisiae genome. **WARNING**: Snoreport2 must 
        be installed manually on your computer before 
        using this rule and you might need to change the
        path in the command line below to access snoreport2."""
    input:
        pos_strand_genome = 'data/references/genome_fa/saccharomyces_cerevisiae_genome.fa',
        neg_strand_genome = 'data/references/genome_fa/negative_strand/saccharomyces_cerevisiae_genome_negative_strand.fa'
    output:
        predictions_strand_pos = 'results/predictions/snoreport2/S_cerevisiae/predicted_cd.fa',
        predictions_strand_neg = 'results/predictions/snoreport2/S_cerevisiae/predicted_cd_negative_strand.fa'
    shell:
        """~/snoReport_2/snoreport_2 -i {input.pos_strand_genome} """
        """-CD """
        """--positives """
        """-o {output.predictions_strand_pos} && """
        """~/snoReport_2/snoreport_2 -i {input.neg_strand_genome} """
        """-CD """
        """--positives """
        """-o {output.predictions_strand_neg} """

rule overlap_snoreport_cerevisiae_CD:
    """ Filter snoreport predictions on S. cerevisiae genome to see the overlap between 
        the predictions and the existing expressed C/D. Return the precision and recall."""
    input:
        snoreport_pos = rules.test_snoreport_cerevisiae.output.predictions_strand_pos,
        snoreport_neg = rules.test_snoreport_cerevisiae.output.predictions_strand_neg,
        cd_sno = 'data/references/positives/cd_rfam_filtered_all_fixed_length_190nt.tsv',
        chr_size = 'data/references/chr_size/saccharomyces_cerevisiae_chr_size.tsv'
    output:
        df = 'results/predictions/snoreport2/S_cerevisiae/precision_recall_cerevisiae.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/overlap_snoreport_cerevisiae_CD.py"

