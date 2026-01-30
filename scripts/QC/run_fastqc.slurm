#!/usr/bin/env bash

#SBATCH --job-name=Blood_sample
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --partition=pibu_el8
#SBATCH --output=./logs/out/fastq_%j.out
#SBATCH --error=./logs/err/fastq_%j.err
#SBATCH --time=02:10:00

#set paths
CONTAINER="/containers/apptainer/fastqc-0.12.1.sif" 
INPUT_DIR="./data/reads_Blood"
OUTPUT_DIR="./results/QC_data"

# Run FastQC using Apptainer container
apptainer exec ${CONTAINER} fastqc \
	${INPUT_DIR}/*fastq.gz \
	-o ${OUTPUT_DIR} \
	-t $SLURM_CPUS_PER_TASK

echo 'FastQC analysis completed!'

