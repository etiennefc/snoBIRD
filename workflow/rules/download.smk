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
    """ Donwload the BERT tokenizer and the pretrained DNA_BERT parameters 
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