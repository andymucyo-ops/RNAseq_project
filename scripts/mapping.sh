#!/usr/bin/env bash

# set paths

INPUT_DIR="/data/users/ankunzimana/RNAseq_project/data/reads_Blood"
SCRIPT_DIR="/data/users/ankunzimana/RNAseq_project/scripts"

for file_1 in `ls -1 ${INPUT_DIR}/*_1.fastq.gz`; do 
    file_2=${file_1%_1.fastq.gz}_2.fastq.gz
    sample_name=$(basename ${file_1%_1.fastq.gz})
    sbatch ${SCRIPT_DIR}/mapping.slurm ${file_1} ${file_2} ${sample_name}
done
