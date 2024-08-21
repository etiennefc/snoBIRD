# SnoBIRD
SnoBIRD (**B**ERT-based **I**dentification and **R**efinement of C/**D** box **sno**RNAs) is a Snakemake-based tool to predict C/D box snoRNA genes across genomic sequences. 

## Tested with the following dependencies

- **Python** version == 3.11.4<br>
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

6 - Finally, download the two models that constitute the core of SnoBIRD:
```bash
python3 snoBIRD.py --download_model
```

## Usage
The most basic use of SnoBIRD is as follows, where [-options] are optional flags to modulate SnoBIRD usage, and input_fasta.fa is your mandatory input sequence (in FASTA format) for which you want SnoBIRD to predict on:
```bash
python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa>
```

You can view as follows all the available options to modulate SnoBIRD usage:
```bash
python3 snoBIRD.py --help
```

**Tips**: You can add the -n/--dryrun option to visualize which steps will be run
by SnoBIRD. This goes through the whole pipeline without actually running it and 
summarizes which steps will be executed. 
```bash
python3 snoBIRD.py [-options] -i </home/your_username/full_path/to/input_fasta.fa> -n
```

**Note**: 
While SnoBIRD can technically run on a local computer without any GPU, its architecture (BERT models and parallelization steps) is designed to be run a HPC cluster using GPUs in order to speed up dramatically its runtime. Therefore, we recommend its use on HPC clusters with GPUs. However, you can realistically run SnoBIRD locally if you have a small input FASTA sequence (<1Mb) (or if you have a lot of spare time ahead of you for sequences of greater size).





