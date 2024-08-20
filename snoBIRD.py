#!/usr/bin/python3
import sys
import subprocess as sp
import argparse



def main(no_arg=False):
    # Set up argument parser
    msg = "python3 snoBIRD.py [-options] <fasta_of_genome.fa>"
    parser = argparse.ArgumentParser(add_help=False, usage=msg,
                    description=(
    "SnoBIRD identifies and refines C/D box snoRNA genes in a given genomic sequence."))
    
    # Set up optional options
    optional_group = parser.add_argument_group('Available options')
    optional_group.add_argument('--help', '-h', action='store_true', 
                                help="Show this help message")
    optional_group.add_argument('--version', '-v', action='store_true', 
                                help="Show SnoBIRD version")
    optional_group.add_argument('--step_size', '-s', type=int, 
            help="Step size between windows for genome-wide predictions (default: 5)", 
                            default=5)
    # --no-chunks is an implicit available flag (--chunks returns True; --no-chunks returns None)
    optional_group.add_argument('--chunks', 
                                action=argparse.BooleanOptionalAction, 
                                default=True,
            help="Divide (--chunks) or not (--no-chunks) large fasta files into smaller chunks/fastas to predict more efficiently in parallel (default: --chunks)")
    optional_group.add_argument('--chunk_size', '-cs', type=int, 
                            help="Maximal chunk size in megabase (Mb) when --chunks True (default: 5000000)", 
                            default=5000000)
    #parser.add_argument('--run_rule', type=str, help="Rule to run")
    #parser.add_argument('--configfile', type=str, help="Configuration file")
    #parser.add_argument('--cores', type=int, help="Number of cores")
    
    # No options or input were given
    if no_arg == True:
        parser.print_help()
        exit()

    # Define args effects
    args = parser.parse_args()
    if args.help:
        parser.print_help()
    if args.version:
        print('SnoBIRD v0.1')

    ## Build Snakemake command
    snakemake_cmd = "cd workflow && snakemake --cores 1 -n "

    ## Add custom parameters to the Snakemake command
    config_l = "--config "
    if args.step_size:
        config_l += f"step_size={args.step_size} "
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
