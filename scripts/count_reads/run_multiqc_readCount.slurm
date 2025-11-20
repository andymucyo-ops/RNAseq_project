#!/usr/bin/env bash

#SBATCH --job-name=multi_qc_RNAseq
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/ankunzimana/RNAseq_project/logs/out/multiqc_%j.out
#SBATCH --error=/data/users/ankunzimana/RNAseq_project/logs/err/multiqc_%j.err
#SBATCH --time=02:10:00
#SBATCH --mail-user=andy.nkunzimana@students.unibe.ch
#SBATCH --mail-type=end,fail

#define the pathways 
CONTAINER="/containers/apptainer/multiqc-1.19.sif"
INPUT_DIR="/data/users/ankunzimana/RNAseq_project/results/count_reads/corrected_count.txt.summary"
OUTPUT_DIR="/data/users/ankunzimana/RNAseq_project/results/count_reads" 

# run scirpt
apptainer exec ${CONTAINER} multiqc \
    ${INPUT_DIR}\
    -o ${OUTPUT_DIR}

echo 'multiqc analysis completed !'