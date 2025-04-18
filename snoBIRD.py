#!/usr/bin/python3
import sys
import subprocess as sp
import argparse
import os
import warnings
warnings.formatwarning = lambda msg, *args: f"\033[93m{msg}\033[0m\n"

snoBIRD_version = 'v0.1'

def find_download(dir_="workflow/data/references/", quiet=False):
    "Find if SnoBIRD models have already been downloaded."
    model_download = "will be downloaded!"
    model_download_bool = False
    if (os.path.exists(dir_+'models')) & (
        os.path.exists(dir_+'DNA_BERT_6_tokenizer')) & (
            os.path.exists(dir_+'DNA_BERT_6_pretrained_model') & (
                os.path.exists('workflow/envs/snoBIRD_env.tar.gz'))
        ):
        potential_downloads_model = os.listdir(dir_+'models')
        potential_downloads_bert = os.listdir(
                                        dir_+'DNA_BERT_6_pretrained_model')
        potential_downloads_tokenizer = os.listdir(
                                        dir_+'DNA_BERT_6_tokenizer')
        if ('snoBIRD_first_model.pt' in potential_downloads_model) & (
            'snoBIRD_second_model.pt' in potential_downloads_model) & (
                'model.safetensors' in potential_downloads_bert) & (
                    'tokenizer.json' in potential_downloads_tokenizer):
            model_download = ("have already been downloaded. You should now "+
                        "run SnoBIRD without the -d/--download_model option.")
            model_download_bool = True
    if quiet == False:
        print(f'SnoBIRD models and env {model_download}')
    return model_download_bool

def arg_value_range(arg_, arg_name, positive=True):
    """ Find if provided arg is within the allowed value range."""
    if positive == True:  # > 0
        if arg_ <= 0:
            raise ValueError(f'Provided {arg_name}: "{arg_}" must be > 0.')
    elif positive == "zero_or_more":  # >= 0
        if arg_ < 0:
            raise ValueError(f'Provided {arg_name}: "{arg_}" must be >= 0.')
    elif positive == "probability":  # must be between 0.5 and 1
        if (arg_ < 0.5) | (arg_ >= 1): 
            raise ValueError(
                f'Provided {arg_name}: "{arg_}" must be >= 0.5 and < 1.')
    elif positive == False: # < 0
        if arg_ > 0:
            raise ValueError(f'Provided {arg_name}: "{arg_}" must be <= 0.')

def is_sbatch_installed():
    """ Find if SLURM's sbatch command is installed to see if the user is on a 
        local computer or a SLURM cluster."""
    try:
        # Try to run sbatch --version to check if sbatch is installed
        result = sp.run(['sbatch', '--version'], stdout=sp.PIPE, 
                    stderr=sp.PIPE, text=True)
        
        # If the command was successful, return True
        if result.returncode == 0:
            return True
        else:
            return False
    except FileNotFoundError:
        # If sbatch is not found, return False
        return False

def verify_input_bed(fasta_path, bed_path):
    """ Verify that bed cols are separated by tab, that the file has at least 
        6 columns, that bed entries are <= 194 nt and that start < end cols."""

    # Check if bed is tab-separated
    tab_cmd = """awk -v FS='\t' '{if (NF < 2) {print $0}}' """ + bed_path + \
                """ | wc -l"""
    untab_lines = sp.run([tab_cmd], shell=True, 
                    capture_output=True, text=True).stdout.strip()
    if int(untab_lines) > 0:
        raise ValueError(f'Your input bed file contains "{untab_lines}" lines'+
            ' that are not tab-separated. Please ensure that all bed entries '+
            'are tab-separated. You can view which lines are tab-separated or'+
            ' not using:\n\tcat -t <input_bed.bed>')
    
    # Check if file has at least 6 columns
    col_cmd = """awk -v FS='\t' '{print NF}' """ + bed_path + \
                """ | sort | uniq """
    num_cols = sp.run([col_cmd], shell=True, 
                    capture_output=True, text=True).stdout.strip()
    num_col_list = list(num_cols.split('\n'))
    if num_col_list == ['']:
        raise ValueError('It seems like you provided an empty input bed file '+
            '(no tab-separated columns were found). Please provide an '+
            'accurate input bed file and make sure the bed path is right.')
    if len(num_col_list) > 1:
        raise ValueError('Your input bed file contains an inconsistant '+
            'number of columns across lines. Lines with the following number '+
            f'of columns were found: {num_col_list}. Please ensure that all '+
            'bed entries have the same number of tab-separated columns.')
    if (len(num_col_list) == 1) & (int(num_col_list[0]) < 6):
        num_i = int(num_col_list[0])
        raise ValueError(f"Your input bed file only contains {num_i} columns,"+
            " which is less than the 6 minimally required columns. Please "+
            "provide a tab-separated input bed file in the following format:"+
            "\n\tchr\tstart\tend\tgene_or_locus_id\tscore\tstrand\nIf you "+
            "don't know the score and/or strand, use '.' as a value for these"+
            " columns.")

    # Check if chr names in bed files match those in input_fasta
    chrbed_cmd = """awk -v OFS='\t' '{print $1}' """+bed_path+ \
                """ | sort | uniq"""
    chrfasta_cmd = f"grep '>' {fasta_path} | sort | uniq "
    bed_chr = sp.run([chrbed_cmd], shell=True, 
                    capture_output=True, text=True).stdout.strip()
    bed_chr = list(bed_chr.split('\n'))
    fa_chr = sp.run([chrfasta_cmd], shell=True, 
                    capture_output=True, text=True).stdout.strip()
    fa_chr = [i.split(' ')[0].replace('>', '') for i in fa_chr.split('\n')]
    missing_chr = []
    for b in bed_chr:
        if b not in fa_chr:
            missing_chr.append(b)
    if len(missing_chr) > 0:
        raise ValueError('Your input bed file contains entries for which the '+
            'chromosome name (1st column) is not present in the input fasta '+
            'file that you provided. Specifically, bed entries with the '+
            'following chromosomes:\n\t'+
            f'{missing_chr}\nare not present in the fasta chromosomes:'+
            f'\n\t{fa_chr}'+
            '\nPlease fix/remove these mismatch entries from your input bed '+
            'file (including the header line if present).')

    # Check if the start is < end for all entries
    loc_cmd = """awk -v FS='\t' '$2>=$3 {print $0}' """ + bed_path + \
                """ | wc -l"""
    loc_coord = sp.run([loc_cmd], shell=True, 
                    capture_output=True, text=True).stdout.strip()
    if int(loc_coord) > 0:
        raise ValueError(f"Your input bed file contains '{loc_coord}' entries"+
            " for which the start is >= to the end value (2nd column >= "+
            "3rd column). Please fix/remove these lines, as the start should "+
            "always be smaller than the end, even for genes on the minus "+
            "strand. You can view these lines with the following command:\n\t"+
            r"awk -v FS='\t' '$2>=$3 {print $0}' <input_bed.bed>")

    # Raise warning for bed entries that are > 194 nt in length
    len_cmd = """awk -v OFS='\t' '$3-$2+1>194 {print $0}' """ + bed_path + \
                """ | wc -l"""
    bigger_regions = sp.run([len_cmd], shell=True, 
                    capture_output=True, text=True).stdout.strip()
    if int(bigger_regions) > 0:
        warnings.warn("UserWarning: Your input bed file contains "+  
            f"'{bigger_regions}' entries of length > 194 nt. As SnoBIRD "+
            "can only predict on windows of size <= 194 nt, only the centered"+
            " window (of length 194 nt) within these larger entries will be "+
            "predicted on by SnoBIRD by default. Consider reducing the size "+
            "of these windows to be <= 194 nt or split them in overlapping "+
            "windows of size 194 nt if you really want them to be included in"+
            " the predictions (see README for more information). You can view"+
            " which entries exceed the window size limit using the following "+
            "command:\n\t"+
            r"awk -v OFS='\t' '$3-$2+1>194 {print $0}' <input_bed.bed>")
    




def main(no_arg=False):
    # Set up argument parser
    msg = ("python3 snoBIRD.py [-options] -i "+
            "</home/your_username/full_path/to/input_fasta.fa>")
    parser = argparse.ArgumentParser(add_help=False, usage=msg,
                    description=(
    "SnoBIRD identifies C/D box snoRNA genes in a given genomic sequence."))
    
    # Set up local and cluster profiles 
    # (cluster by default, unless -n or -d is used)
    profile_local = '--profile ../profile_local '
    profile_slurm = '--profile ../profile_slurm '

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
        help="Download the models that SnoBIRD uses (should be used only "+
        "once, before running SnoBIRD with an input fasta file)")
    optional_group.add_argument('--dryrun', '-n', action='store_true', 
        help="Run a snakemake dry run to just print a summary of the DAG of "+
        "jobs and verify job dependencies")
    optional_group.add_argument('--local_profile', '-L', action='store_true', 
        help="Force SnoBIRD to run locally, i.e. not on a HPC cluster. This "+
            "should ONLY be used for small input fasta files (<1Mb). By "+
            "default, SnoBIRD assumes that you use it on a HPC cluster (since"+
            " it was designed to run in parallel on several GPUs) and will "+
            "deal automatically with job dependencies. This local_profile "+
            "option is thus disabled by default, unless you use the -n and/or"+
            " -d options (which will run locally a dryrun or download_model, "+
            "respectively).")
    optional_group.add_argument('--cores', '-k', type=int, 
        help="Number of CPU cores to provide for parallelizing SnoBIRD's "+
            "steps across different cores (default: 1). WARNING: this "+
            "option only has an effect when SnoBIRD is run locally on "+
            "your computer (-L option), i.e. not on a HPC cluster. Although "+
            "the number of asked CPU and GPU cores is already optimized when "+
            "SnoBIRD runs on a HPC cluster, if you really want to override "+
            "these values for SnoBIRD's different steps on a HPC cluster, you "+
            "can modify them manually in the file "+
            "workflow/profile_slurm/cluster.yaml", default=1)
    optional_group.add_argument('--gpu_generation', '-G', type=str, 
        choices=['H100', 'A100', 'V100', 'P100', 'Unknown'], 
        help='GPU generation (from NVIDIA) that will be used to run SnoBIRD '+
        '(default: A100). WARNING: if Unknown is set, SnoBIRD identification '+
        'step (genome_prediction) will take longer to start running, as the '+
        'same generic longer time will be asked for SnoBIRD to run on each '+
        'chr/chunks. In addition, if Unknown is set, the shap_snoBIRD step '+
        'will also take longer to start running because a generic default '+
        'longer time is allocated for that step. However, one can '+
        'change these generic longer times by changing the time keyword in '+
        'the genome_prediction and/or shap_snoBIRD entries in the file '+
        'workflow/profile_slurm/cluster.yaml. To find which GPU generations '+
        'are available on your SLURM cluster, use the command '+
        'sinfo -o "%%N %%G" | grep "[pPaAvVhH]100"', default="A100")
    optional_group.add_argument('--first_model_only', '-f', action='store_true', 
        help="Run only the first SnoBIRD model and not the second (i.e. "+
        "predict only the presence of C/D snoRNA genes in the general sense; "+
        "no predictions if they are expressed snoRNAs or snoRNA pseudogenes "+
        "by the second model)")
    # --no-chunks is an implicit available flag 
    # (--chunks returns True; --no-chunks returns None)
    optional_group.add_argument('--chunks', 
                                action=argparse.BooleanOptionalAction, 
                                default=True,
        help="Divide (--chunks) or not (--no-chunks) large fasta files into"+
            " smaller chunks/fastas to predict more efficiently in parallel; "+
            "if --no-chunks, SnoBIRD will be less efficient "+
            "(default: --chunks)")
    optional_group.add_argument('--chunk_size', '-cs', type=int, 
        help="Maximal chunk size in megabases (Mb) when --chunks is chosen "+
            "(default: 5)", default=5)
    optional_group.add_argument('--step_size', '-s', type=int, 
        help="Step size between predicted windows (default: 5); where s=1 "+
        "will predict on ALL possible windows contained in a given sequence; "+
        "WARNING: setting s=1 will increase SnoBIRD's runtime by ~5 fold "+
        "compared to the default value", default=5)
    optional_group.add_argument('--strand', '-S', type=str, 
        choices=['both', 'plus', 'minus'], 
        help="Strand on which to predict (default: both)", default="both")
    optional_group.add_argument('--batch_size', '-b', type=int, 
        help='Number of sequences per batch that are processed in parallel '+
            'by SnoBIRD (default: 128); increasing this number increases '+
            'VRAM/RAM usage', default=128)
    optional_group.add_argument('--output_name', '-o', type=str, 
        help="Desired output file name (no extension needed) that will be "+
            "located in the workflow/results/final/ directory "+
            "(default: snoBIRD_complete_predictions)",
            default="snoBIRD_complete_predictions")
    optional_group.add_argument('--output_type', type=str, 
        choices=['tsv', 'fa', 'bed', 'gtf'],
        help="Desired output file type, i.e. either a tab-separated (.tsv), "+
            "bed (.bed) or fasta (.fa) file (default: tsv)", default="tsv")
    optional_group.add_argument('--input_bed', type=str, 
        help="COMPLETE (absolute) path to input bed file containing "+
        "specific genomic locations for SnoBIRD to predict on. This is "+
        "useful for example if you have identified peaks (regions) of "+
        "interest in an experiment (CLIP, RNA-Seq, etc.) and want to find if "+
        "these regions harbor a C/D box snoRNA without having to run SnoBIRD "+
        "on the entire genome sequence (thereby reducing significantly "+
        "SnoBIRD's runtime). Entries in the bed file must be <=194 nt and "+
        "you must still provide an input fasta file of the genome with "+
        "--input_fasta. The input bed should NOT include a header as the "+
        "first line. When using --input_bed, the following options are "+
        "ignored: --chunks/--no-chunks, --chunk_size, --step_size, --strand, "+
        "--consecutive_windows")
    optional_group.add_argument('--prob_first_model', '-p1', type=float, 
        help="Minimal prediction probability of a given window to be "+
            "considered as a C/D snoRNA gene by SnoBIRD's first model in the "+
            "identification step (default: 0.999, or 0.9 if --input_bed is "+
            "specified or if all input sequences are of length <= 194 nt); "+
            "WARNING: lowering this number can dramaticaly increase the "+
            "number of false positives, whereas increasing this number can "+
            "increase the number of false negatives ", default=0.999)
    optional_group.add_argument('--consecutive_windows', '-w', type=int, 
        help="Minimum number of consecutive windows predicted as a C/D snoRNA"+
            " gene by SnoBIRD's first model to be considered as a C/D snoRNA "+
            "gene in the identification step (default: 10); in other words, "+
            "one needs at least 10 positive consecutive windows overlapping a"+
            " candidate for it to be considered as a C/D snoRNA gene", 
            default=10)
    optional_group.add_argument('--prob_second_model', '-p2', type=float, 
        help="Minimal prediction probability of a given window to be "+
            "confidently considered as either an expressed C/D snoRNA or a "+
            "C/D snoRNA pseudogene by SnoBIRD's second model in the "+
            "refinement step (default: 0.999)", default=0.999)
    optional_group.add_argument('--box_score', '-B', type=int, 
        help="Box score threshold value used in the refinement step to filter"+
            " SnoBIRD's second model predictions; it represents the maximal "+
            "number of mutations across boxes tolerated in an expressed C/D "+
            "snoRNA (default: 5)", default=5)
    optional_group.add_argument('--score_c', '-C', type=int, 
        help="Threshold value used in the refinement step to filter"+
            " SnoBIRD's second model predictions; it represents the maximal "+
            "number of mutations within the C box tolerated in an expressed "+
            "C/D snoRNA (default: 2)", default=2)
    optional_group.add_argument('--score_d', '-D', type=int, 
        help="Threshold value used in the refinement step to filter"+
            " SnoBIRD's second model predictions; it represents the maximal "+
            "number of mutations within the D box tolerated in an expressed "+
            "C/D snoRNA (default: 1)", default=1)
    optional_group.add_argument('--terminal_stem_score', '-T', type=float, 
        help="Threshold value used in the refinement step to filter"+
            " SnoBIRD's second model predictions; it represents the maximal "+
            "terminal stem score tolerated in an expressed C/D snoRNA "+
            "(default: -25); lowering this threshold (e.g. -100) will return "+
            "more stringent predictions for expressed C/D snoRNA genes",
             default=-25.0)
    optional_group.add_argument('--normalized_sno_stability', '-N', type=float, 
        help="Threshold value used in the refinement step to filter"+
            " SnoBIRD's second model predictions; it represents the maximal "+
            "normalized secondary structure stability (in "+
            "kcal/mol/nucleotide) tolerated in an expressed C/D snoRNA "+
            "(default: -0.2); lowering this threshold (e.g. -0.5) will return"+
            " more stringent predictions for expressed C/D snoRNA genes",
             default=-0.2)
    optional_group.add_argument('--unlock', '-u', action='store_true', 
        help="Unlock working directory. This option will be necessary only "+
        "when a kill signal or a power loss happens while the jobs are "+
        "dispatched on the cluster (thereby creating a remaining lock that "+
        "interferes with snakemake). This should be used once (which will "+
        "unlock the working directory), then rerun the snoBIRD command "+
        "without this option")
             
    

    
    args = parser.parse_args()

    # Define help and version args effects
    if args.help:
        parser.print_help()
        exit()
    if args.version:
        print(f'SnoBIRD {snoBIRD_version}')
        exit()

    # No options or input were given
    if no_arg == True:
        parser.print_help()
        exit()

    # Verify if custom arg is in the right value range
    arg_value_range(args.cores, "cores")
    arg_value_range(args.chunk_size, "chunk_size")
    arg_value_range(args.step_size, "step_size")
    arg_value_range(args.batch_size, "batch_size")
    arg_value_range(args.consecutive_windows, "consecutive_windows")
    arg_value_range(args.box_score, "box_score", "zero_or_more")
    arg_value_range(args.score_c, "score_c", "zero_or_more")
    arg_value_range(args.score_d, "score_d", "zero_or_more")
    arg_value_range(args.prob_first_model, "prob_first_model", "probability")
    arg_value_range(args.prob_second_model, "prob_second_model", "probability")
    arg_value_range(args.terminal_stem_score, "terminal_stem_score", False)
    arg_value_range(args.normalized_sno_stability, "normalized_sno_stability", 
                    False)
    # Verify if output_name is <= 100 characters 
    # (this is for shap_snoBIRD time reallocation to work properly)
    if len(str(args.output_name)) > 100:
        raise ValueError(f'Provided output_name: "{args.output_name}" length '+
        'cannot be > 100 characters. Please reduce your output_name length to'+
        ' <= 100 characters.')


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
                print('SnoBIRD models and env need to be downloaded first '+
                        'with:\npython3 snoBIRD.py --download_model')
                exit()
        if args.input_bed:
            config_l += f"input_bed={args.input_bed} "
            if '/' not in args.input_bed:
                raise argparse.ArgumentTypeError(f'You must provide the '+ 
                   f'complete (absolute) path to {args.input_bed}; for '+
                    'example, </home/your_username/full_path/to/input_bed.bed>'
                    )
    else:
        if args.download_model:
            snakemake_cmd = "cd workflow && snakemake "
            snakemake_cmd += "all_downloads " + profile_local
            config_l += "input_fasta=fake_input "
            config_l += "download_model=True "
            config_l += "profile=local "
            find_download()
        elif args.input_bed:
            parser.error('An input fasta file of your genome is also required'+
                        ' when using the --input_bed option. '+
                        'Please use -i <input_fasta.fa>')
        else:
            parser.error('An input fasta file is required. '+
                        'Please use -i <input_fasta.fa>')


    # Define optional args effects

    ## Add custom parameters to the Snakemake command
    if args.local_profile:
        snakemake_cmd = "cd workflow && snakemake " + profile_local +\
                        f"--cores {args.cores} "
        config_l += "profile=local "
        
    if args.input_fasta:
        if args.download_model:
            snakemake_cmd = "cd workflow && snakemake "
            snakemake_cmd += "all_downloads " + profile_local
            config_l += "profile=local "
            find_download()
    if args.dryrun:
        snakemake_cmd = "cd workflow && snakemake "+ profile_local
        snakemake_cmd += "-n "
        config_l += "dryrun=True " 
        if not args.download_model:
            if not args.input_bed:
                print("\nExecuting the dryrun (getting the number of chr "+
                "and/or chunks of chr that will be created)...This may take "+
                "a bit of time for large genomes (ex: <1 min for the human "+
                "genome).\n")
            else:
                print("\nExecuting the dryrun (verifying your input bed and "+
                "getting the steps that will run)...\n")
                verify_input_bed(args.input_fasta, args.input_bed)
        else:
            snakemake_cmd = "cd workflow && snakemake all_downloads -n "+ \
                            profile_local 
    if not args.local_profile:  # default value
        if not args.download_model and not args.dryrun:
            if is_sbatch_installed():
                snakemake_cmd = "cd workflow && snakemake " + \
                                profile_slurm
            else:
                raise ModuleNotFoundError("It seems that you are not on a "+
                "SLURM cluster, as sbatch is not installed. If you want to "+
                "run SnoBIRD locally, add the --local_profile option to your "+
                "command.")

    if args.step_size:
        config_l += f"step_size={args.step_size} "
    if args.strand:
        config_l += f"strand='{args.strand}' "
    if args.batch_size:
        config_l += f"batch_size={args.batch_size} "
    if args.chunks:
        config_l += f"chunks={args.chunks} "
    if args.output_name:
        config_l += f"output_name={args.output_name} "
    if args.output_type:
        config_l += f"output_type={args.output_type} "
    if args.gpu_generation:
        config_l += f"gpu_generation={args.gpu_generation} "
    if args.prob_first_model:
        if args.input_bed:
            if args.prob_first_model != 0.999:  # user-defined
                config_l += (
                f"min_probability_threshold_first_model={args.prob_first_model} ")
            else:  # less stringent threshold as less windows are predicted on 
                config_l += (
                f"min_probability_threshold_first_model={0.9} ")
        else:
            config_l += (
                f"min_probability_threshold_first_model={args.prob_first_model} ")
    if args.consecutive_windows:
        config_l += (
            f"min_consecutive_windows_threshold={args.consecutive_windows} ")
    if args.prob_second_model:
        config_l += (
        f"min_probability_threshold_second_model={args.prob_second_model} ")
    if args.box_score:
        config_l += f"box_score_threshold={args.box_score} "
    if args.score_c:
        config_l += f"score_c_threshold={args.score_c} "
    if args.score_d:
        config_l += f"score_d_threshold={args.score_d} "
    if args.terminal_stem_score:
        config_l += (
                f"terminal_stem_score_threshold={args.terminal_stem_score} ")
    if args.normalized_sno_stability:
        config_l += (
        f"normalized_sno_stability_threshold={args.normalized_sno_stability} ")
    if args.chunk_size:
        if args.chunks:
            # Convert chunk_size from Megabase to base
            config_l += f"chunk_size={args.chunk_size * 1000000} "
        else:
            config_l += "chunk_size=None "
    if args.first_model_only:
        config_l += "first_model_only=True "
    else:
        config_l += "first_model_only=False "
    
    snakemake_cmd += config_l

    if args.unlock:
        snakemake_cmd += "--unlock "
    
    ## Run Snakemake
    sp.call(snakemake_cmd, shell=True)

    ## Return the actual command line that was run for reproducibility
    if (args.download_model == False) & (args.unlock == False) & (
        args.dryrun == False):
        no_return_var = ['help', 'version', 'download_model', 'dryrun', 
                        'unlock']  # args that are not reproducibility-related
        if args.first_model_only:  # don't return args w/r to 2nd model
            no_return_var.extend(['prob_second_model', 'box_score', 'score_c', 
                'score_d', 'terminal_stem_score', 'normalized_sno_stability'])
        if args.input_bed:  # no args w/r to split_chr and merge_filter_windows
            no_return_var.extend(['chunks', 'no-chunks', 'chunk_size', 
                    'step_size', 'consecutive_windows', 'strand'])

        if not args.local_profile:  # don't return 'cores' arg if cluster usage
            no_return_var.append('cores')

        cmd_line = "python3 snoBIRD.py "
        dict_var = vars(args)
        output_n = dict_var['output_name']
        output_t = dict_var['output_type']
        for k,v in dict_var.items():
            if k not in no_return_var:
                cmd_line += f"--{k} {v} "
        return cmd_line, output_n, output_t

if __name__ == "__main__":
    if len(sys.argv) == 1:
        main(no_arg=True)
    else:
        sp.call('mkdir -p workflow/logs/', shell=True)
        result = main()
        if result is not None:
            arguments, output_final, ext_ = result
            sp.call('echo "####SnoBIRD report usage" >> snoBIRD_usage.log', 
                                                                    shell=True)
            sp.call('date >> snoBIRD_usage.log', shell=True)
            sp.call('echo "##[Used command line]:" >> snoBIRD_usage.log', 
                                                                    shell=True)
            sp.call(f'echo {arguments} >> snoBIRD_usage.log', shell=True)
            sp.call('echo "##[Desired output file location]:" >> '+
                    'snoBIRD_usage.log', shell=True)
            sp.call(f'echo workflow/results/final/{output_final}.{ext_} >> '+
                    'snoBIRD_usage.log', shell=True)
            sp.call('echo "##[SnoBIRD version]:" >> snoBIRD_usage.log', 
                                                                    shell=True)
            sp.call(f'echo {snoBIRD_version} >> snoBIRD_usage.log', shell=True)
            sp.call('echo "##[Used packages versions]:" >> snoBIRD_usage.log', 
                                                                    shell=True)
            if is_sbatch_installed():
                sp.call(
                    'workflow/snoBIRD_env/bin/pip list >> snoBIRD_usage.log', 
                    shell=True)
            else:
                sp.call(
                    'conda list -p workflow/snoBIRD_env/ >> snoBIRD_usage.log', 
                    shell=True)
            sp.call('echo "\n\n" >> snoBIRD_usage.log', shell=True)
        
