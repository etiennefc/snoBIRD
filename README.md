<div style="text-align: justify">  

# SnoBIRD
SnoBIRD (**B**ERT-based **I**dentification and **R**efinement of C/**D** box **sno**RNAs) is a Snakemake-based tool to predict C/D box snoRNA genes across genomic sequences. 

## Tested with the following dependencies

- **Python** version == 3.10.14<br>
- **Conda** version == 23.7.2<br>
- **Mamba** version == 1.4.9<br>
- **Snakemake** version == 7.18.2<br>

## Installation
1 - Conda (Miniconda3) needs to be installed (https://docs.conda.io/en/latest/miniconda.html) via the following commands:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Answer `yes` to `Do you wish the installer to initialize Miniconda3?`

2 - Mamba needs to be installed via conda (mamba greatly speeds up environment creation, but conda can still be used instead of mamba):
```bash
conda install -n base -c conda-forge mamba
```

3 - To create the snakemake environment that SnoBIRD uses, run the following commands:

```bash
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake=7.18
```

4 - Activate snakemake environment needed for SnoBIRD to run:
```bash
conda activate snakemake
```

5 - Download the SnoBIRD repository:
```bash
git clone https://github.com/etiennefc/snoBIRD.git &&
cd snoBIRD/
```

6 - Finally, download the models and tokenizer that constitute the core of SnoBIRD, as well as create the environment in which SnoBIRD will run:
```bash
python3 snoBIRD.py --download_model
```

## Usage
The most basic use of SnoBIRD is as follows, where [-options] are optional flags to modulate SnoBIRD usage, and `input_fasta.fa` is your mandatory input sequence (in FASTA format) for which you want SnoBIRD to predict on and for which you must provide the **full path**:
```bash
python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa>
```
By default, this command will run **both** the first and second model of SnoBIRD, which will identify in your input sequence C/D box snoRNA genes (first model) and refine these predicted C/D box snoRNA genes by predicting if they are expressed snoRNAs or snoRNA pseudogenes (second model). SnoBIRD also assumes, by default, that it is run on High-Performance Computing (HPC) cluster. 

Your final predictions will be located in the file `workflow/results/final/snoBIRD_complete_predictions.tsv`, a tab-separated file (.tsv) by default.

You can view as follows all the available options to modulate SnoBIRD usage:
```bash
python3 snoBIRD.py --help
```

## Tips 
- ***Dry Run***: You can add the `-n/--dryrun` option to visualize which steps will be run
by SnoBIRD. This goes through the whole pipeline without actually running it and 
summarizes which steps will be executed. 
    ```bash
    python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> -n
    ```
- ***Output***: You can specify which output type you want using the `--output_type` option (default: tsv). The default output format is a tab-separated (.tsv) file containing all predicted snoRNAs and other relevant information (genomic location, box positions, structural feature values, prediction probability, snoRNA sequence, etc.). The other available output formats are FASTA and BED files, which also contain all the aforementioned snoRNA information. In addition, you can also modify the default output file name `snoBIRD_complete_predictions` using the `-o/--output_name` option (no need to provide the file extension in the name, SnoBIRD automatically adds it depending on the chosen `--output_type` option). Of note, the default output file directory (`workflow/results/final/`) cannot be changed, only the file name can.
    ```bash
    python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> --output_type <tsv|fa|bed> -o <your_favorite_file_name>
    ```
- ***Run First Model Only***: If you want to run only the first model of SnoBIRD, i.e. that you are interested in knowing only if there are C/D box snoRNA genes in your input sequence and therefore you don't want to differentiate these predictions between if they are expressed or snoRNA pseudogenes, you should run the following command:
    ```bash
    python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> -f
    ```
    Of note, this will only slightly reduce SnoBIRD's overall runtime, as the most time-consuming step occurs during the first model's prediction over the entire input sequence. 
- ***GPU Usage***: There are different generations of GPUs that might be available on your SLURM cluster (H100, A100, V100, P100, etc.). You can view (highlight) which GPU generation is available as follows:
    ```bash
    sinfo -o "%N %G" | grep "[pPaAvVhH]100"
    ```
    The latest GPU type is usually far more efficient/faster than its predecessor (H100 > A100 > V100 > P100). For SnoBIRD usage, we recommend using minimally V100 GPUs (although A100 and H100 GPUs will be far more efficient, if they're ever available to you). One could use P100 GPUs, but they will show a significant increase in terms of running time and should only be used as a last resort and for smaller genomes. By default, SnoBIRD assumes that it will be run using NVIDIA's A100 GPUs, but you can specify which GPU generation you will use (in order to optimize the time that SnoBIRD spends using the GPU) as follows:
    ```bash
    python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> -G <H100|A100|V100|P100|Unknown>
    ```
- ***Cluster Resources Allowance***: The resources asked for each step (rule) in the SnoBIRD pipeline are defined by default in the `profile_slurm/cluster.yaml` file. For example, the allowed memory, number of CPU/GPU and time are defined for each task. You can change these parameters directly in the file if necessary (e.g. if your job runs out of time, consider increasing the allowed time for the given rule). The only parameter that cannot be overriden directly is the `time` parameter in the `genome_prediction` rule, which predicts the presence of C/D box snoRNA genes (in the general sense) with SnoBIRD's first model on a given chromosome or chunk of chromosome. By default, the allowed time per chromosome or chunk of chrosomome for this rule is already optimized based on the type of GPU (`-G`) as well as the chunk size (`-cs`) and step size (`-s`) that you provided. If you want to override this behavior, you can change the time value for `genome_prediction` in the `profile_slurm/cluster.yaml` file (default: 11 hours) and you must add the following `-G` option:
```bash
python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> -G Unknown
```  

## Notes 
While SnoBIRD can technically run on a local computer without any GPU, its architecture (BERT models and parallelization steps) is designed to be run a HPC cluster using GPUs in order to speed up dramatically its runtime. Therefore, we recommend its use on HPC clusters with GPUs (SnoBIRD assumes by default that it will be run on a HPC cluster). However, you can realistically run SnoBIRD on a local computer if you have a small input FASTA sequence (<1Mb) (or if you have a lot of spare time ahead of you for sequences of greater size). In that case, you should add the `-L/--local_profile` option to your command so that SnoBIRD works properly on a local computer. When using the `-L/--local_profile` option, you can also specify the number of CPU cores that you want to provide locally to increase parallelism between SnoBIRD's steps (and therefore decrease runtime) using the `-k/--cores` option (default: 1).
```bash
python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> -L -k <number_of_cores_to_use>
```


</div>

