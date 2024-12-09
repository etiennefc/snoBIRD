rule download_models:
    """ Download first and second SnoBIRD models from Zenodo."""
    output:
        model1 = "data/references/models/snoBIRD_first_model.pt",
        model2 = "data/references/models/snoBIRD_second_model.pt"
    params:
        link1 = config['download']['first_model'],
        link2 = config['download']['second_model']
    shell:
        "echo 'Downloading SnoBIRD models:' && "
        "wget {params.link1} -O {output.model1} && "
        "wget {params.link2} -O {output.model2}"

rule download_DNA_BERT:
    """ Download the BERT tokenizer and the pretrained DNA_BERT parameters 
        from Zenodo (originally from 'zhihan1996/DNA_bert_6' in 
        huggingface)."""
    output:
        tokenizer = directory("data/references/DNA_BERT_6_tokenizer/"),
        dnabert = directory("data/references/DNA_BERT_6_pretrained_model/")
    params:
        link1 = config['download']['DNA_bert_6_tokenizer'],
        link2 = config['download']['DNA_bert_6_pretrained_model']
    shell:
        "echo 'Downloading DNA_BERT_6 parameters:' && mkdir -p data/references"
        " && wget {params.link1} -O temp_tokenizer.tar.gz && "
        "wget {params.link2} -O temp_dnabert.tar.gz && "
        "tar -xzf temp_tokenizer.tar.gz && "
        "mv DNA_BERT_6_tokenizer/ {output.tokenizer} && "
        "tar -xzf temp_dnabert.tar.gz && "
        "mv DNA_BERT_6_pretrained_model/ {output.dnabert} && "
        "rm temp_dnabert.tar.gz temp_tokenizer.tar.gz"

rule create_env:
    """ Create the environment in which all rules will be run. A conda env is 
        created if SnoBIRD is run locally. If SnoBIRD is run on a cluster, a 
        virtualenv (with pip) is instead created, as the equivalent conda env 
        of this virtualenv hinders GPU usage. We create a tar archive so that 
        we can copy it on cluster nodes directly more efficently and then 
        activate that environment on the node directly (which is overall more 
        efficient)."""
    output:
        env = "envs/snoBIRD_env.tar.gz"
    params:
        cluster = is_sbatch_installed2(),
        virtualenv = config["virtualenv"],
        bedtool_link = config["download"]["bedtool"],
        condaenv = config["condaenv"]
    shell:
        "if [ {params.cluster} = True ]; then "
        "module load python && "
        "echo 'Creating snoBIRD_env (it can take a while, please be patient)...' &&"
        "virtualenv --download snoBIRD_env && "
        "source snoBIRD_env/bin/activate && "
        "echo 'Downloading packages in snoBIRD_env...' && "
        "pip install --quiet -r {params.virtualenv} && "
        "wget {params.bedtool_link} && "
        "mv bedtools.static snoBIRD_env/bin/bedtools && "
        "chmod a+x snoBIRD_env/bin/bedtools && "
        "deactivate && "
        "echo 'Converting the virtualenv in a tar archive (please be patient)...' && "
        "tar -czf {output.env} snoBIRD_env && "
        "echo 'Environment creation is complete!'; else "
        "echo 'Creating snoBIRD_env and downloading packages in it...' && "
        "mamba env create -q -f {params.condaenv} -p snoBIRD_env && "
        "touch {output.env} && echo 'Environment creation is complete!'; "
        "fi"
