rule download_models:
    """ Download first and second SnoBIRD models."""
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