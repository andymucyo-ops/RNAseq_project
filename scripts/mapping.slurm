#!/usr/bin/env bash

#SBATCH --job-name=Hisat2_index_RNAseq
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/ankunzimana/RNAseq_project/logs/out/Hisat2_index_%j.out
#SBATCH --error=/data/users/ankunzimana/RNAseq_project/logs/err/Hisat2_index_%j.err
#SBATCH --time=02:00:00
#SBATCH --mail-user=andy.nkunzimana@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

#set paths
CONTAINER="containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
INDEX_DIR="/data/users/ankunzimana/RNAseq_project/results/hisat2_indexes"
FILE_1= $1
FILE_2= $2
SAMPLE_NAME= $3
OUTPUT_DIR="/data/users/ankunzimana/RNAseq_project/results/alignment/sam"

# mapping script 

apptainer exec ${CONTAINER} hisat2 -p $SLURM_CPUS_PER_TASK\
    -x ${INDEX_DIR} \
    -1 ${FILE_1} \
    -2 ${FILE_2} \
    -S ${OUTPUT_DIR}/${SAMPLE_NAME}.sam \
    --summary-file ${OUT_DIR}/${SAMPLE_NAME}_summary.txt\
    --rna-strandness RF 

echo'mapping completed!'