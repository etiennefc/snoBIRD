
![SnoBIRD_logo_small](https://github.com/user-attachments/assets/05f17c18-3b79-4eab-916d-743affdb5f7e)

<div style="text-align: justify">  

# SnoBIRD
SnoBIRD (**B**ERT-based **I**dentification and **R**efinement of C/**D** box **sno**RNAs) is a Snakemake-based tool to predict C/D box snoRNA genes across genomic sequences. 

## Tested with the following dependencies

- **Python** version == 3.10.14<br>
- **Conda** version == 23.7.2<br>
- **Mamba** version == 1.4.9<br>
- **Snakemake** version == 7.18.2<br>

## Installation
We recommend installing SnoBIRD on the scratch space of your High-Performance Computing (HPC) cluster. This storage space usually contains more memory than your regular project space (which fits well SnoBIRD's storage space requirements as it creates large intermediate files), although it should only be used temporarily (usually in terms of weeks). Therefore, once SnoBIRD is installed and has run successfully, you should move your output predictions to a long-term storage space disk.

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
You can also provide a BED file in addition to your input fasta. This is useful if you already suspect regions in your genome of interest to contain snoRNAs (e.g. you get a BED file after peak calling following CLIP- or RNA-Seq experiments). Using the `--input_bed` option as follows and providing the **full path** to the BED file, SnoBIRD will only predict on the regions present in your BED file rather than on the whole input fasta, thereby dramatically decreasing SnoBIRD's runtime:
```bash
python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> --input_bed </home/your_username/full_path/to/input_bed.bed>
```

By default, these commands will run **both** the first and second model of SnoBIRD, which will identify in your input sequence C/D box snoRNA genes (first model) and refine these predicted C/D box snoRNA genes by predicting if they are expressed snoRNAs or snoRNA pseudogenes (second model). SnoBIRD also assumes, by default, that it is run on a HPC cluster. Therefore, these previous command lines should be run **directly on the login node** of your HPC cluster (SnoBIRD will automatically submit your jobs on computing nodes with pre-optimized parameters for each job (e.g. memory, number of CPU/GPU, time, etc.)).

Your final predictions will be located in the file `workflow/results/final/snoBIRD_complete_predictions.tsv`, a tab-separated file (.tsv) by default. All your log files will be located in the directory `workflow/logs/`.

You can view as follows all the available options to modulate SnoBIRD usage:
```bash
python3 snoBIRD.py --help
```

We <span style="color:green;">**STRONGLY RECOMMEND**</span> that you read at least the *Dry Run* and *GPU Usage* subsections in the **Tips** section below, as it will likely impact SnoBIRD performance and usage. Feel free to read the rest of the section to optimize your use of SnoBIRD.

## Tips 
- ***Dry Run***: We **strongly recommend** to use the dryrun option before actually running SnoBIRD. To do so, you can add the `-n/--dryrun` option which will help you visualize which steps will be run
by SnoBIRD. 
    ```bash
    python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> -n
    ```
    This goes through the whole pipeline without actually running it and 
summarizes which steps will be executed. For instance, you will see how many `genome_prediction` jobs will run and if you want to adjust that number by playing with the chunk size (`-cs`) and/or step size (`-s`) parameters, depending on your GPU availabilities. This also shows you how many fasta entries are present in your input sequence and can help you spot unwanted sequences before running SnoBIRD on them unnecessarily (e.g. there are usually a lot of scaffold sequences at the end of genome sequences that you probably don't want to predict on).
- ***GPU Usage***: There are different generations of GPUs that might be available on your SLURM cluster (H100, A100, V100, P100, etc.). You can view (highlight) which GPU generation is available as follows:
    ```bash
    sinfo -o "%N %G" | grep "[pPaAvVhH]100"
    ```
    The latest GPU type is usually far more efficient/faster than its predecessor (H100 > A100 > V100 > P100). For SnoBIRD usage, we recommend using minimally V100 GPUs (although A100 and H100 GPUs will be far more efficient, if they're ever available to you). One could use P100 GPUs, but they will show a significant increase in terms of running time and should only be used as a last resort and for smaller genomes. By default, <span style="color:red;">**SnoBIRD assumes that it will be run using NVIDIA's A100 GPUs**</span>, but you can specify which GPU generation you will use (in order to optimize the time that SnoBIRD spends using the GPU) as follows:
    ```bash
    python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> -G <H100|A100|V100|P100|Unknown>
    ```
    Some HPC clusters can harbor different GPU types (e.g. A100 or V100 depending on the node) or can offer different GPU configurations (e.g. multi-instance GPU that split a GPU between users). To specify which GPU type/configuration to use in these cases, you can modify or add specific sbatch options in the `profile_slurm/cluster.yaml` file. Be sure to change it for all jobs that use a GPU, i.e. `merge_filter_windows`, `genome_prediction`, `predict_and_filter_bed_windows`, `shap_snoBIRD` and `sno_pseudo_prediction`.
- ***Input BED***: When using the `--input_bed` option, you must also provide an input fasta of your genome using `--input_fasta`. This serves the purpose of retrieving the sequence of the entries in your input BED file. By default, only BED entries with a length <= 194 nt will be predicted on by SnoBIRD (as SnoBIRD was trained on sequences of this size only). If entries in your BED file exceed 194 nt in length, SnoBIRD will, by default, predict only on the centered window (of 194 nt in length) of that larger entry. For example, if the initial entry has a length of 214 nt, SnoBIRD will only predict on the 194 nt window starting after the 10th nucleotide and ending at the 10th before last nucleotide. To override this behavior, you can remove your initial long BED entry, convert it to as many overlapping windows of length 194 nt as you want and add these windows of appropriate size back to your initial BED file (e.g. convert a 214 nt window into 21 overlapping windows of length 194 nt using a step size of 1). In addition, if no strand is provided in the BED file, SnoBIRD will predict by default on both + and - strand at the given genomic location, thereby duplicating that given BED entry. 
- ***Output***: You can specify which output type you want using the `--output_type` option (default: tsv). The default output format is a tab-separated (.tsv) file containing all predicted snoRNAs and other relevant information (genomic location, box positions, structural feature values, prediction probability, snoRNA sequence, etc.). The other available output formats are FASTA, BED and GTF files, which also contain all the aforementioned snoRNA information. In addition, you can also modify the default output file name `snoBIRD_complete_predictions` using the `-o/--output_name` option (no need to provide the file extension in the name, SnoBIRD automatically adds it depending on the chosen `--output_type` option). Of note, the default output file directory (`workflow/results/final/`) cannot be changed, only the file name can.
    ```bash
    python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> --output_type <tsv|fa|bed> -o <your_favorite_file_name>
    ```
- ***Run First Model Only***: If you want to run only the first model of SnoBIRD, i.e. that you are interested in knowing only if there are C/D box snoRNA genes in your input sequence and therefore you don't want to differentiate these predictions between if they are expressed or snoRNA pseudogenes, you should run the following command:
    ```bash
    python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> -f
    ```
    Of note, this will only slightly reduce SnoBIRD's overall runtime, as the most time-consuming step occurs during the first model's prediction over the entire input sequence. 

- ***Cluster Resources Allocation***: The resources asked for each step (rule) in the SnoBIRD pipeline are defined by default in the `profile_slurm/cluster.yaml` file. For example, the allocated memory, number of CPU/GPU and time are defined for each task. You can change these parameters directly in the file if necessary (e.g. if your job runs out of time, consider increasing the allocated time for the given rule). The only parameters that cannot be overriden directly are the `time` parameter in the `genome_prediction` and `predict_and_filter_bed_windows` rules, which predict the presence of C/D box snoRNA genes (in the general sense) with SnoBIRD's first model on a given chromosome/chunk of chromosome or on specified bed intervals respectively, as well as the `time` parameter in the `shap_snoBIRD` rule, which computes SHAP values for all predicted snoRNAs. By default, the allocated time for these rules is already optimized based on the following parameters: type of GPU (`-G`), chunk size (`-cs`), step size (`-s`), total number of bed intervals (for the `predict_and_filter_bed_windows` step only) or total number of predicted snoRNAs (for the `shap_snoBIRD` step only). If you want to override these behaviors, you can change the time value for `genome_prediction`, `predict_and_filter_bed_windows` and/or `shap_snoBIRD` in the `profile_slurm/cluster.yaml` file (default: 11, 3 and 15 hours respectively) and you must add the following `-G Unknown` option:
    ```bash
    python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> -G Unknown
    ```  
- ***HPC Account Usage***: Some users might have more than one account (linked with a given resource allocation) associated to their username on a HPC cluster. By default, SnoBIRD assumes that you only have one default account and this might result in the following sbatch error when trying to use SnoBIRD:
    ```bash
    sbatch: error: You are associated with multiple _cpu allocations...
    sbatch: error: Please specify one of the following accounts to submit this job:
    sbatch: error:   RAS default accounts: <def_account1>, <def_account2>, 
    sbatch: error: Use the parameter --account=desired_account when submitting your job
    sbatch: error: Batch job submission failed: Unspecified error
    ```
    To fix this error, modify the script `workflow/scripts/python/slurmSubmit.py` by uncommenting and modifying the following line accordingly:
    ```bash
    #cmdline = "sbatch --account=[def-your_account] "
    ```

## Notes 
While SnoBIRD can technically run on a local computer without any GPU, its architecture (BERT models and parallelization steps) is designed to be run a HPC cluster using GPUs in order to speed up dramatically its runtime. Therefore, we recommend its use on HPC clusters with GPUs (SnoBIRD assumes by default that it will be run on a HPC cluster). However, you can realistically run SnoBIRD on a local computer if you have a small input FASTA sequence (<1Mb) (or if you have a lot of spare time ahead of you for sequences of greater size). In that case, you should add the `-L/--local_profile` option to your command so that SnoBIRD works properly on a local computer. When using the `-L/--local_profile` option, you can also specify the number of CPU cores that you want to provide locally to increase parallelism between SnoBIRD's steps (and therefore decrease runtime) using the `-k/--cores` option (default: 1).
```bash
python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> -L -k <number_of_cores_to_use>
```


</div>

