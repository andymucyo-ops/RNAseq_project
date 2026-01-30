# RNAseq FS 2025

This project analyzes RNA-seq data from mouse blood samples infected with Toxoplasma gondii to identify differentially expressed genes (DEGs) between wild-type (WT) and double-knockout (DKO) genotypes under case (infected) and control (uninfected) conditions. The workflow includes quality control, read alignment, quantification, and comprehensive differential expression analysis using DESeq2, followed by Gene Ontology enrichment analysis and interferon module profiling.

## PART 1.

This part has been conducted on the IBU cluster, where data have been retrieved from and all the following scripts have been run, to asses reads quality, map them to the reference genome and count the reads obtain to produce the inputs for the Differential expression analysis in R.


| Tool | Description | Documentation |
|------|-------------|---------------|
| FastQC | Quality control for raw sequencing reads | [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| MultiQC | Aggregate QC reports across samples | [MultiQC Documentation](https://multiqc.info/) |
| HISAT2 | Spliced alignment of RNA-seq reads | [HISAT2 Manual](http://daehwankimlab.github.io/hisat2/) |
| Samtools | SAM/BAM file manipulation and sorting | [Samtools Documentation](http://www.htslib.org/doc/samtools.html) |
| featureCounts | Read quantification at gene level | [featureCounts Documentation](https://subread.sourceforge.net/) |


The data used for this project were retireved from the following repository on the IBU cluster 
    
    /data/courses/rnaseq_course/toxoplasma_de/reads_Blood
But can also be fond on the Gene Expression Omnibus, with the accession number GSE119855 ([GEO link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119855)) 
### 1) Quality Check of the data

A frist script has been run to run the fastqc tool on all raw data (`*.fastq.gz`)

    sbatch ./scripts/QC/run_fastqc.slurm

Then a second scirpt using multiqc for applied on the fastqc output for better visualisation of the quality check of all reads 
    
    sbatch ./scripts/QC/run_multiqc.slurm

### 2) MAP reads to the reference genome

The reference genome and associated annotations were downloaded form the following website: https://www.ensembl.org/info/data/ftp/index.html 
Using the following commands: 

For reference genome: 

    wget https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz ./data/reference_genome

For associated annotaions:

    wget https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz ./data/reference_genome

The reference genome is then unziped with 
```{bash}
gzip -d ./data/reference_genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
```
and the name is then changed to 
```{bash} 
mv ./data/reference_genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa ./data/reference_genome/genome.fa
```

Then run the following script to produce the indexes for mapping:

    sbatch ./scripts/alignment/Hist2_index.slurm 

to obtain the sam file for the aligement run:

    ./scripts/alignment/mapping.sh

that will run the `mapping.slurm` script on every raw fastq data file

Finally run 

    ./scripts/alignment/samtools.sh

that will run the `samtools.slurm` on the .sam files. This scripts consits of three actions 
- 1) conversion of .sam file, to .bam file (compresssed form containing the informations)
- 2) sorting the .bam files 
- 3) indexing the sorted .bam files

### 3) Read counting

run the following script on the sorted bam files:

    sbatch scripts/count_reads/count_reads.slurm

## PART 2.

This part is conducted on your local machine, for that you'll need to clone this repo 

```{bash}
git clone https://github.com/andymucyo-ops/RNAseq_project
```

- The requirment to conduct the analysis are listed in the [environment.yml file](./docker-build/environment.yml)

- Other requirments are either Mamba (installed using miniforge or miniconda) or Docker, installed on your local machine

### 4) Differential expression analysis 

- first clone the repo on your local machine in a repo of your choice 
```{bash}
git clone https://github.com/andymucyo-ops/RNAseq_project
```

and enter directory `cd RNAseq_project`

Then the analysis can be run in two ways:

#### 1) conda environment + R Studio
- first install miniforge (conda/mamba installer) following the steps on official [Mamba documentation](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

- Then, create conda environment from `docker-build/environment.yml`, and activate it with the following commands:
    
    ```{bash}
    mamba env create -f docker-build/environment.yml
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

##### Outputs

**Tables** (`results/tables/`):
- `DEG_*.csv`: Complete differential expression results for each comparison
- `DEG_summary.csv`: Summary of total DEGs, upregulated, and downregulated genes
- `*_GO_complete.csv`: Full GO enrichment results
- `*_GO_top15.csv`: Top 15 enriched GO terms

**Plots** (`results/R_plots/`):
- `PCAplot.png`: Principal component analysis
- `heatmap_*.png`: Expression heatmaps for each comparison
- `volcano_plot_*.png`: Volcano plots for each comparison
- `dotplot_enrichGO_*.png`: GO enrichment visualization
- `heatmap_B11_module.png`: B11 interferon module expression
- `heatmap_B14_module.png`: B14 interferon module expression

