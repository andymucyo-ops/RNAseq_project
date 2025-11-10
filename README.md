# RNAseq FS 2025

The data used for this project were retireved from the following repository 

    /data/courses/rnaseq_course/toxoplasma_de/reads_Blood

1) Quality Check of the data

A frist script has been run to run the fastqc tool on all data

    ./scripts/run_fastqc.slurm

Then a second scirpt using multiqc for applied on the fastqc output for better visualisation of all 
reads 
    
    ./scripts/run_multiqc.slurm

2) MAP reads to the reference genome

The reference genome and associated annotations were downloaded form the following website: https://www.ensembl.org/info/data/ftp/index.html 
Using the following commands: 

For reference genome: 

    wget ftp://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz

For associated annotaions:

    wget ftp://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz

The reference genome is then unziped with `bash gzip -d ./data/reference_genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz` and the name is then changed to `bash mv ./data/reference_genome/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa ./data/reference_genome/genome.fa`

Then produce the required index running the `bash sbatch ./scripts/Hist2_index.slurm` script
