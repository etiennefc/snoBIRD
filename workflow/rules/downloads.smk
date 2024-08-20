rule download_blat:
    """ Download BLAT to find coordinates 
        of a given sequence in a genome 
        (in fasta). Must enter sudo password at 
        a given point in the download."""
    output:
        tmp_file = "data/references/blat/blat_test.txt"
    params:
        blat = config['download']['blat']
    shell:
        "mkdir -p ./data/references/blat/ && "
        "rsync -aP {params.blat} ./data/references/blat/ &> {output.tmp_file}"

rule download_rfam_covariance_models:
    """ Download the Rfam library of covariance models that 
        will be used by Infernal."""
    output:
        rfam_cm = 'data/references/RFam.cm'
    params:
        link = config['download']['rfam_cm']
    shell:
        "wget {params.link} && "
        "gunzip Rfam.cm.gz && mv Rfam.cm {output.rfam_cm}"

rule download_genome:
    """ Download fasta of genome per species from Ensembl"""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["leishmania_major", 
              "dictyostelium_discoideum", "giardia_lamblia", "arabidopsis_thaliana", 
              "oryza_sativa", "aspergillus_fumigatus", "neurospora_crassa", "candida_albicans"], remove=True)
    params:
        link = "ftp://ftp.ensembl.org/pub/release-108/fasta/{species}/dna/*dna.toplevel.fa.gz",
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_mammal_genome:
    """ Download the reference genome (fasta file) of human and mouse
        from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species_tgirt'], ["mus_musculus",    
              "homo_sapiens"])
    params:
        link = "ftp://ftp.ensembl.org/pub/release-108/fasta/{species}/dna/*dna.primary_assembly.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_yeast_genome:
    """ Download the reference genome (fasta file) of S. pombe, 
        S. cerevisiae, N. crassa, C. albicans and A. fumigatus from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["saccharomyces_cerevisiae",
              "aspergillus_fumigatus", "neurospora_crassa", "candida_albicans"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-55/fasta/{species}/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_tetrahymena_genome:
    """ Download fasta of Tetrahymena thermophila genome per species from Zenodo/TGD"""
    output:
        genome = 'data/references/genome_fa/tetrahymena_thermophila_genome.fa'
    params:
        link = config['download']['t_thermophila_genome']
    shell:
        "wget -O {output.genome} {params.link}"

rule download_other_protist_genome:
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["leishmania_major", 
              "dictyostelium_discoideum", "giardia_lamblia"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/protists/release-55/fasta/{species}/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_plant_genome:
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["arabidopsis_thaliana", "oryza_sativa"])
    params:
        link = "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-55/fasta/{species}/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_o_tauri_genome:
    """ Download fasta of Ostreococcus tauri genome per species from Zenodo 
        (originally by concatenating the 20 chromosomes sequences from NCBI, 
        id ranging from NC_014426.2 to NC_014445.2)."""
    output:
        genome = 'data/references/genome_fa/ostreococcus_tauri_genome.fa'
    params:
        link = config['download']['o_tauri_genome']
    shell:
        "wget -O {output.genome} {params.link}"

rule download_mouse_gtf:
    """ Download the annotation of the mouse genome (gtf)
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/mus_musculus.gtf'
    params:
        link = config['download']['mouse_gtf']
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_yeast_gtf:
    """ Download the reference genome (fasta file) of S. cerevisiae, 
        S. pombe, A. fumigatus, N. crassa and C. albicans
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ["saccharomyces_cerevisiae", 
                "aspergillus_fumigatus", "neurospora_crassa", "candida_albicans"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-55/gtf/{species}/*5.gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_human_gtf:
    """ Download gtf of human genome from Zenodo. Remove trailing tabs."""
    output:
        gtf = 'data/references/gtf/homo_sapiens.gtf'
    params:
        link = config['download']['human_gtf']
    shell:
        "wget -O {output.gtf} {params.link} && "
        "sed -i 's/[\t]$//g' {output.gtf}"

rule download_tetrahymena_gtf:
    """ Download gtf of Tetrahymena thermophila genome from Zenodo. Remove trailing tabs."""
    output:
        gtf = 'data/references/gtf/tetrahymena_thermophila.gtf'
    params:
        link = config['download']['tetrahymena_gtf']
    shell:
        "wget -O {output.gtf} {params.link} && "
        "sed -i 's/[\t]$//g' {output.gtf}"

rule download_protist_gtf:
    """ Download the annotation (gtf) of different protists 
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ["leishmania_major", 
              "dictyostelium_discoideum", "giardia_lamblia"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/protists/release-55/gtf/{species}/*[^chr].gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"
    
rule download_plant_gtf:
    """ Download the annotation (gtf) of different plants
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ["oryza_sativa", "arabidopsis_thaliana"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/plants/release-55/gtf/{species}/*[^chr].gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_animal_gtf:
    """ Download the annotation (gtf) of different animals
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ['macaca_mulatta', 'ornithorhynchus_anatinus', 
                    'gallus_gallus', 'caenorhabditis_elegans', 'drosophila_melanogaster'])
    params:
        link = "ftp://ftp.ensembl.org/pub/release-108/gtf/{species}/*[^chrabinitio].gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_rnacentral_ncRNA:
    """ Download rnacentral bed file of all ncRNAs per species 
        (will be useful for tRNAs, snRNAs and pre-miRNA). 
        **Ostreococcus tauri is not present in rnacentral """
    output:
        bed = 'data/references/rnacentral/{species}.bed'
    wildcard_constraints:
        species=join_list(config['species']+config['species_tgirt'], 
                ["ostreococcus_tauri", "schizosaccharomyces_pombe"], remove=True)
    params:
        link = config['download']['rnacentral'] 
    shell:
        "wget -O temp_rnacentral_{wildcards.species}.gz {params.link}{wildcards.species}.*.bed.gz && "
        "gunzip temp_rnacentral_{wildcards.species}.gz && "
        "mv temp_rnacentral_{wildcards.species} {output.bed}"

rule download_gallus_gallus_gff:
    """ Download from rnacentral the latest version of Gallus gallus gff (bGalGal1/GRCg7b) to 
        convert the location of the old version (which is used in RNAcentral)"""
    output:
        gff = "data/references/rnacentral/gallus_gallus_bGalGal1_GRCg7b.gff3"  
    params:
        link = config['download']['gallus_gallus_gff']
    shell:
        "wget -O gff_gallus.gz {params.link} && "
        "gunzip gff_gallus.gz && mv gff_gallus {output.gff}"

rule download_eukaryote_HACA_snoRNAs:
    """ Download all H/ACA present in RNAcentral (v21) for all eukaryotes. It was done 
        manually by searching for snoRNAs and then selecting H/ACA snoRNAs (89 489 sequences). 
        ***Need to put that fasta on Zenodo for posterity***"""
    output:
        HACA_fa = 'data/references/HACA/HACA_eukaryotes.fa'  # to put on Zenodo
    params:
        link = config['download']['haca_rnacentral']
    shell:
        "wget -O {output.HACA_fa} {params.link}"

rule download_eukaryote_tRNAs:
    """ Download all tRNAs present in GtRNAdb (Release 20) for all eukaryotes"""
    output:
        tRNA_df = 'data/references/tRNA/tRNA_eukaryotes.tsv'  # to put on Zenodo
    params:
        link = 'http://gtrnadb.ucsc.edu/download/GtRNAdb/search/gtrnadb-search187174.out'  
    shell:
        "wget -O {output.tRNA_df} {params.link}"

rule download_rfam_clans:
    """ Download the clans in RFam, i.e. ~superfamilies (Release 14.9)"""
    output:
        df = 'data/references/rfam/fam_clans_families.tsv' 
    params:
        link = 'https://ftp.ebi.ac.uk/pub/databases/Rfam/14.9/database_files/clan_membership.txt.gz'  
    shell:
        "wget -O {output.df}.gz {params.link} && "
        "gunzip {output.df}.gz"

rule download_rDNA:
    """ Download ribosomal DNA fasta file per species 
        (it's needed to run snoscan). """ # TO PUT ON ZENODO!!!
    output:
        rDNA_fa = 'data/references/rDNA/rDNA_species.fa'
    params:
        link = config['download']['rDNA'] 
    shell:
        "wget -O {output.rDNA_fa} {params.link}"

rule download_all_genes_pombe:
    """ Download from Pombase all genes and their description for S. pombe. 
        This is to retrieve all potential genes known as snoRNAs."""
    output:
        df = 'data/references/s_pombe/all_genes_s_pombe.tsv'
    params:
        link = config['download']['all_genes_pombe'] 
    shell:
        "wget -O {output.df} {params.link}"

rule download_CD_pombe_review:
    """ From Fafard-Couture et al. 2024 RNA Biology, get the 
        snoRNA type for snoRNAs in S. pombe."""
    output:
        df = 'data/references/s_pombe/cd_s_pombe_review.tsv'
    params:
        link = config['download']['cd_s_pombe_review']
    shell:
        "wget -O {output.df} {params.link}"

rule download_sno_type_info:
    """ Download snoRNA type per Rfam family. The snotype per Rfam 
        family was obtained manually by saving (right-click) the 
        Rfam webpage displaying all C/D (and then H/ACA) snoRNA 
        families and parsing the rfam family id only using the following command: 
        grep ">RF" ~/Downloads/HACA_RFAM_families_2024.html | sed -E 's/.*">//g; s/<.*//g'"""
    output:
        sno_type_rfam = 'data/references/snoRNA_type_rfam_families.tsv'
    params:
        link = config['download']['sno_type_rfam']
    shell:
        "wget -O {output.sno_type_rfam} {params.link}"
