rule download_pombe_gtf:
    """ Download the reference gtf file of S. pombe
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/schizosaccharomyces_pombe.gtf'
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-55/gtf/schizosaccharomyces_pombe/*5.gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_pombe_genome:
    """ Download the reference genome (fasta file) of S. pombe, 
        from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/schizosaccharomyces_pombe_genome.fa'
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-55/fasta/schizosaccharomyces_pombe/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule get_s_pombe_CD_snoRNAs:
    """ From Pombase, get the gene_id of snoRNAs in 
        S. pombe. From the gtf, extract the coordinates 
        of these snoRNAs (we don't filter directly from 
        the gtf because it's simpler than dealing with 
        the many snoRNA genes duplicates/overlaps). From RNA central, 
        find if these snoRNAs are C/D box snoRNAs."""
    input:
        s_pombe_genes = rules.download_all_genes_pombe.output.df,
        gtf = rules.download_pombe_gtf.output.gtf,
        genome = rules.download_pombe_genome.output.genome,
        sno_type_pombe = rules.download_CD_pombe_review.output.df
    output:
        bed = 'data/references/s_pombe/filtered_CD_snoRNAs.bed',
        fa = 'data/references/s_pombe/filtered_CD_snoRNAs.fa'
    params:
        link = 'https://rnacentral.org/link/pombase:'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_s_pombe_CD_snoRNAs.py"

rule test_rfam_infernal_pombe:
    """ Use Infernal and Rfam covariance 
        models to predict on the genome of S. pombe to see if some sequences 
        are part of a C/D Rfam family."""
    input:
        test_genome = 'data/references/genome_fa/schizosaccharomyces_pombe_genome.fa',
        rfam_cm = rules.download_rfam_covariance_models.output.rfam_cm
    output:
        infernal_tblout = 'results/predictions/infernal_rfam/S_pombe/predicted_cd.tblout',
        infernal_alignments = 'results/predictions/infernal_rfam/S_pombe/predicted_cd.txt'
    conda:
        "../envs/infernal.yaml"
    shell:
        """cmpress -F {input.rfam_cm} && """
        """cmscan --cut_ga --rfam --nohmmonly -o {output.infernal_alignments} """
        """--tblout {output.infernal_tblout} {input.rfam_cm} {input.test_genome}"""

rule test_snoscan_pombe:
    """ Predict with snoscan the presence of 
        C/D in S. pombe genome."""
    input:
        target_rDNA = 'data/references/rDNA/s_pombe_rDNA.fa',
        test_genome = 'data/references/genome_fa/schizosaccharomyces_pombe_genome.fa'
    output:
        predictions = 'results/predictions/snoscan/S_pombe/predicted_cd.txt'
    conda: 
        "../envs/snoscan.yaml"
    shell:
        "snoscan {input.target_rDNA} {input.test_genome} -o {output.predictions}"

rule test_snoreport_pombe:
    """ Test the performance of snoreport2 on S pombe genome. **WARNING**: Snoreport2 must 
        be installed manually on your computer before 
        using this rule and you might need to change the
        path in the command line below to access snoreport2."""
    input:
        pos_strand_genome = 'data/references/genome_fa/schizosaccharomyces_pombe_genome.fa',
        neg_strand_genome = 'data/references/genome_fa/negative_strand/schizosaccharomyces_pombe_genome_negative_strand.fa'
    output:
        predictions_strand_pos = 'results/predictions/snoreport2/S_pombe/predicted_cd.fa',
        predictions_strand_neg = 'results/predictions/snoreport2/S_pombe/predicted_cd_negative_strand.fa'
    shell:
        """~/snoReport_2/snoreport_2 -i {input.pos_strand_genome} """
        """-CD """
        """--positives """
        """-o {output.predictions_strand_pos} && """
        """~/snoReport_2/snoreport_2 -i {input.neg_strand_genome} """
        """-CD """
        """--positives """
        """-o {output.predictions_strand_neg} """


rule get_pred_bed_other_tool_pombe:
    """ Create bed out of predictions of all other tools on the S. pombe genome. """
    input:
        snoreport_pos = rules.test_snoreport_pombe.output.predictions_strand_pos,
        snoreport_neg = rules.test_snoreport_pombe.output.predictions_strand_neg,
        snoscan = rules.test_snoscan_pombe.output.predictions,
        infernal = rules.test_rfam_infernal_pombe.output.infernal_tblout,
        sno_rfam = rules.download_sno_type_info.output.sno_type_rfam,
        genome_fa = expand(rules.download_yeast_genome.output.genome, species='schizosaccharomyces_pombe')
    output:
        bed_infernal = "results/predictions/infernal_rfam/S_pombe/predicted_cd.bed",
        bed_snoscan = "results/predictions/snoscan/S_pombe/predicted_cd.bed",
        bed_snoscan_filtered_target = "results/predictions/snoscan/S_pombe/predicted_cd_1_target.bed",
        bed_snoreport = "results/predictions/snoreport2/S_pombe/predicted_cd_all_strand.bed"
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_pred_bed_other_tool_pombe.py"


rule intersect_pombe_preds:
    """ For all tools, get the number of predictions that overlap with already annotated C/D snoRNAs. """
    input:
        snoreport = rules.get_pred_bed_other_tool_pombe.output.bed_snoreport,
        snoscan = rules.get_pred_bed_other_tool_pombe.output.bed_snoscan_filtered_target,
        infernal = rules.get_pred_bed_other_tool_pombe.output.bed_infernal,
        snoBIRD = 'results/predictions/snoBIRD/S_pombe/filtered_preds_step5.bed',
        cd_pombe = 'data/references/positives/pombe_CD_snoRNAs.bed'
    output:
        bed_infernal = "results/predictions/infernal_rfam/S_pombe/overlap_w_annotated_cd.bed",
        bed_snoscan = "results/predictions/snoscan/S_pombe/overlap_w_annotated_cd.bed",
        bed_snoreport = "results/predictions/snoreport2/S_pombe/overlap_w_annotated_cd.bed",
        bed_snoBIRD = "results/predictions/snoBIRD/S_pombe/overlap_w_annotated_cd.bed"
    params:
        fixed_length = 194
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/intersect_pombe_preds.py"

rule igv_batch_screenshot_pombe:
    """ From a bed, create a batch script which will be fed to IGV 
        to automate screenshot creation at each locus in the bed file."""
    input:
        snoreport_bed = rules.get_pred_bed_other_tool_pombe.output.bed_snoreport,
        snoscan_bed = rules.get_pred_bed_other_tool_pombe.output.bed_snoscan_filtered_target,
        infernal_bed = rules.get_pred_bed_other_tool_pombe.output.bed_infernal,
        snoBIRD_bed = 'results/predictions/snoBIRD/S_pombe/filtered_preds_step5.bed'
    output:
        batch_script_snoscan = 'results/igv_scripts/snoscan_s_pombe.batch',
        batch_script_snoreport = 'results/igv_scripts/snoreport_s_pombe.batch',
        batch_script_infernal = 'results/igv_scripts/infernal_rfam_s_pombe.batch',
        batch_script_snoBIRD = 'results/igv_scripts/snoBIRD_s_pombe.batch'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/igv_batch_screenshot_pombe.py"


rule upset_tools_pombe:
    """ From the snoBIRD prediction on S. pombe genome and all of the other tools, 
        find which predictions overlap, which are expressed or not and create an 
        upset plot out of it. """
    input:
        snoreport = rules.get_pred_bed_other_tool_pombe.output.bed_snoreport,
        snoscan = rules.get_pred_bed_other_tool_pombe.output.bed_snoscan_filtered_target,
        infernal_rfam = rules.get_pred_bed_other_tool_pombe.output.bed_infernal,
        snoBIRD = "results/predictions/snoBIRD/S_pombe/filtered_preds_step5.bed"  # from pombe_filter_windows on lab32GB
    output:
        upset = 'results/figures/upset/tool_comparison_s_pombe.svg',
        df = 'results/predictions/s_pombe_all_tools_upset_data.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/upset_tools_pombe.py"



rule bar_cerevisiae_pred:
    """ Create a stacked bar chart of the total number of predictions (out of which the 
        number of annotated vs not annotated C/D)."""
    #input:
        #cerevisiae_pred_all_tools = # To link here
    output:
        bar = 'results/figures/barplot/cerevisiae_predictions.svg'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/bar_cerevisiae_pred.py"

