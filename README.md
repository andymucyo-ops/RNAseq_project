# RNAseq FS 2025

The data used for this project were retireved from the following repository 

    /data/courses/rnaseq_course/toxoplasma_de/reads_Blood

## 1) Quality Check of the data

A frist script has been run to run the fastqc tool on all raw data

    sbatch ./scripts/QC/run_fastqc.slurm

Then a second scirpt using multiqc for applied on the fastqc output for better visualisation of the quality check of all reads 
    
    sbatch ./scripts/QC/run_multiqc.slurm

## 2) MAP reads to the reference genome

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

## 3) Read counting

run the following script ont the sorted bam files:

    sbatch scripts/count_reads/count_reads.slurm
    
## 4) Exploratory analysis

### Step 1) reformat the data as expected by DESeq2
  
  apply the following command on the featurecount output file to create a new file, without the first line and the unwanted column

    tail -n +2 corrected_count.txt | cut -f 7-21 > count.txt

  and same for the README in the reads_blood directory, to create a new file containing the sample information

    tail -n +10 README > sample_info.txt


### Step 2) Differential expression analysis 
  
  run the following script
    
    DES_data_analysis.R 





