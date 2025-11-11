#!/usr/bin/env bash

#SBATCH --job-name=Hisat2_index_RNAseq
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/ankunzimana/RNAseq_project/logs/out/Hisat2_index_%j.out
#SBATCH --error=/data/users/ankunzimana/RNAseq_project/logs/err/Hisat2_index_%j.err
#SBATCH --time=01:00:00
#SBATCH --mail-user=andy.nkunzimana@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

#set pathways 
INPUT_DIR="/data/users/ankunzimana/RNAseq_project/data/reference_genome/*.fa"
OUTPUT_DIR="/data/users/ankunzimana/RNAseq_project/results/hisat2_indexes/"
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"



# run hisat2-built

apptainer exec ${CONTAINER} hisat2-build \
    -f ${INPUT_DIR} \
    ${OUTPUT_DIR} \
    -p $SLURM_CPUS_PER_TASK

echo 'Hisat2 index has been created!'
