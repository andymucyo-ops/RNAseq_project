# RNAseq FS 2025

## PART 1.

This part has been conducted on the IBU cluster, where data have been retrieved from and all the following scripts have been run, to asses reads quality, map them to the reference genome and count the reads obtain to produce the inputs for the Differential expression analysis in R.


The data used for this project were retireved from the following repository 
    
    /data/courses/rnaseq_course/toxoplasma_de/reads_Blood

### 1) Quality Check of the data

A frist script has been run to run the fastqc tool on all raw data

    sbatch ./scripts/QC/run_fastqc.slurm

Then a second scirpt using multiqc for applied on the fastqc output for better visualisation of the quality check of all reads 
    
    sbatch ./scripts/QC/run_multiqc.slurm

### 2) MAP reads to the reference genome

The reference genome and associated annotations were downloaded form the following website: https://www.ensembl.org/info/data/ftp/index.html 
Using the following commands: 

For reference genome: 

    wget ftp://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz ./data/reference_genome

For associated annotaions:

    wget ftp://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz ./data/reference_genome

The reference genome is then unziped with `bash gzip -d ./data/reference_genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz` and the name is then changed to `bash mv ./data/reference_genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa ./data/reference_genome/genome.fa`

Then run the following script to produce the indexes for mapping:

    sbatch ./scripts/alignment/Hist2_index.slurm 

to obtain the sam file for the aligement run:

    ./scripts/alignment/mapping.sh

that will run the `bash mapping.slurm` script on every raw fastq data file

Finally run 

    ./scripts/alignment/samtools.sh

that will run the `bash samtools.slurm` on the .sam files. This scripts consits of three actions 
    1) conversion of .sam file, to .bam file (compresssed form containing the informations)
    2) sorting the .bam files 
    3) indexing the sorted .bam files

### 3) Read counting

run the following script ont the sorted bam files:

    sbatch scripts/count_reads/count_reads.slurm

## PART 2.

### 4) Differential expression analysis 

first clone the repo on your local machine in a repo of your choice `git clone https://github.com/andymucyo-ops/RNAseq_project`

and enter directory `cd RNAseq_project`

Then the analysis can be run in two ways:

#### 1) conda environment + R Studio
- first install miniforge (conda/mamba installer) following the steps on official [Mamba documentation](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

- Then, create conda environment from `docker-build/environment.yml`, and activate it with the following commands:
    
    ```{bash}
    (mamba create -f docker-build/environment.yml)
    ```
    ```{bash}
    conda activate rnaseq-r
    ```
    
- once the environment is activated you can open R Studio and execute the DES_data_analysis.R script
- ⚠️ When using R Studio make sure that it is using R from rnaser-r environment, run the following command in R Studio console: `R.home()` 

#### 2) Docker container

- pull docker image: 
    ```{bash}
    docker pull andymucyo/rnaseq_r:1.0
    ```

- then run the following command from the project root (RNAseq_project/):
    ```{bash}
    docker run --rm -v $(pwd):/project -w /project andymucyo/rnaseq_r:1.0 Rscript scripts/R_scripts/DES_data_analysis.R
    ```

- outputs will be found in `RNAseq_project/results/R_plots` and `RNAseq_project/results/tables` 
