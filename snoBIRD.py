#!/usr/bin/python3
import sys
import subprocess as sp
import argparse
import os

def find_download(dir_="workflow/data/references/models/", quiet=False):
    "Find if SnoBIRD models have already been downloaded."
    model_download = "will be been downloaded!"
    model_download_bool = False
    if os.path.exists(dir_):
        potential_downloads = os.listdir(dir_)
        if ('snoBIRD_first_model.pt' in potential_downloads) & (
            'snoBIRD_second_model.pt' in potential_downloads):
            model_download = ("have already been downloaded. You should now "+
                        "run SnoBIRD without the -d/--download_model option.")
            model_download_bool = True
    if quiet == False:
        print(f'SnoBIRD models {model_download}')
    return model_download_bool


def main(no_arg=False):
    # Set up argument parser
    msg = ("python3 snoBIRD.py [-options] -i "+
            "</home/your_username/full_path/to/input_fasta.fa>")
    parser = argparse.ArgumentParser(add_help=False, usage=msg,
                    description=(
    "SnoBIRD identifies C/D box snoRNA genes in a given genomic sequence."))
    
    # Set up required options
    required_group = parser.add_argument_group('Required options')
    required_group.add_argument('--input_fasta', '-i', type=str, 
        help="COMPLETE (absolute) path to input fasta file containing the "+
        "sequence(s) for SnoBIRD to predict on. It can contain one or "+
        "multiple entries marked by '>' in the same file")
    
    # Set up optional options
    optional_group = parser.add_argument_group('Available options')
    optional_group.add_argument('--help', '-h', action='store_true', 
                                help="Show this help message")
    optional_group.add_argument('--version', '-v', action='store_true', 
                                help="Show SnoBIRD version")
    optional_group.add_argument('--download_model', '-d', action='store_true', 
        help="Download the 2 models that SnoBIRD uses "+
        "(should be used only once, before running snoBIRD)")
    optional_group.add_argument('--dryrun', '-n', action='store_true', 
        help="Run a snakemake dry run to just print a summary of the DAG of "+
        "jobs and verify job dependencies")
    optional_group.add_argument('--step_size', '-s', type=int, 
        help="Step size between predicted windows (default: 5); where s=1 "+
        "will predict on ALL possible windows contained in a given sequence", 
                                default=5)
    # --no-chunks is an implicit available flag 
    # (--chunks returns True; --no-chunks returns None)
    optional_group.add_argument('--chunks', 
                                action=argparse.BooleanOptionalAction, 
                                default=True,
        help="Divide (--chunks) or not (--no-chunks) large fasta files into"+
            " smaller chunks/fastas to predict more efficiently in parallel; "+
            "if --no-chunks, SnoBIRD is less efficient "+
            "(default: --chunks)")
    optional_group.add_argument('--chunk_size', '-cs', type=int, 
        help="Maximal chunk size in megabases (Mb) when --chunks is chosen "+
            "(default: 5000000)", default=5000000)
    optional_group.add_argument('--strand', '-S', type=str, 
        choices=['both', '+', '-'], 
        help="Strand on which to predict (default: both)", default="both")
    #parser.add_argument('--run_rule', type=str, help="Rule to run")
    #parser.add_argument('--configfile', type=str, help="Configuration file")
    #parser.add_argument('--cores', type=int, help="Number of cores")
    
    args = parser.parse_args()

    # Define help and version args effects
    if args.help:
        parser.print_help()
        exit()
    if args.version:
        print('SnoBIRD v0.1')
        exit()

    # No options or input were given
    if no_arg == True:
        parser.print_help()
        exit()

    ## Build Snakemake command
    snakemake_cmd = ("cd workflow && snakemake --use-conda --cores 1 "+
                    "--rerun-triggers mtime --conda-frontend mamba ")

    # Define the required args effects
    config_l = "--config "
    if args.input_fasta:
        config_l += f"input_fasta={args.input_fasta} "
        if '/' not in args.input_fasta:
            raise argparse.ArgumentTypeError(f'You must provide the complete'+ 
                   f' (absolute) path to {args.input_fasta}; for example, '+
                    '</home/your_username/full_path/to/input_fasta.fa>')
        if find_download(quiet=True) == False:
            if not args.download_model:
                print('SnoBIRD models need to be downloaded first with:\n'+
                        'python3 snoBIRD.py --download_model')
                exit()
    else:
        if args.download_model:
            snakemake_cmd += "all_downloads "
            config_l += "input_fasta=fake_input "
            find_download()
        else:
            parser.error('An input fasta file is required. '+
                        'Please use --i <input_fasta.fa>')


    # Define optional args effects

    ## Add custom parameters to the Snakemake command
    if args.input_fasta:
        if args.download_model:
            snakemake_cmd += "all_downloads "
            find_download()
    if args.dryrun:
        snakemake_cmd += "-n "
        print("\nExecuting the dryrun (getting the number of chr and/or "+
                "chunks of chr that will be created)...This may take a bit of"+
                " time for large genomes (ex: <1 min for the human genome).\n")
    if args.step_size:
        config_l += f"step_size={args.step_size} "
    if args.strand:
        config_l += f"strand={args.strand} "
    if args.chunks:
        config_l += f"chunks={args.chunks} "
    if args.chunk_size:
        if args.chunks:
            config_l += f"chunk_size={args.chunk_size} "
        else:
            config_l += "chunk_size=None "
    
    snakemake_cmd += config_l
    
    ## Run Snakemake
    sp.call(snakemake_cmd, shell=True)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        main(no_arg=True)
    else:
        main()
