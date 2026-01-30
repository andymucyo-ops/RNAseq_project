#!/usr/bin/env bash

#SBATCH --job-name=Hisat2_align_RNAseq
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=pibu_el8
#SBATCH --output=./logs/out/Hisat2_align_%j.out
#SBATCH --error=./logs/err/Hisat2_align_%j.err
#SBATCH --time=01:00:00
#SBATCH --mail-user= <your email>
#SBATCH --mail-type=end,fail

#set paths
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
INDEX_BASENAME="./results/hisat2_indexes/index"
FILE_1=${1}
FILE_2=${2}
SAMPLE_NAME=${3}
OUTPUT_DIR="./results/alignment/sam"

# mapping script 

apptainer exec ${CONTAINER} hisat2 -p $SLURM_CPUS_PER_TASK\
    -x ${INDEX_BASENAME} \
    -1 ${FILE_1} -2 ${FILE_2} \
    -S ${OUTPUT_DIR}/${SAMPLE_NAME}.sam \
    --summary-file ${OUTPUT_DIR}/summary/${SAMPLE_NAME}_summary.txt\
    --rna-strandness RF 

echo 'mapping completed!'
