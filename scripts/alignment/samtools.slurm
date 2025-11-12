#!/usr/bin/env bash

#SBATCH --job-name=samtools_RNAseq
#SBATCH --mem=32G
#SBATCH --cpus-per-task=6
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/ankunzimana/RNAseq_project/logs/out/samtools_%j.out
#SBATCH --error=/data/users/ankunzimana/RNAseq_project/logs/err/samtools_%j.err
#SBATCH --time=02:00:00
#SBATCH --mail-user=andy.nkunzimana@students.unibe.ch
#SBATCH --mail-type=end,fail


#set paths

SAM_FILE=${1}

OUTPUT_FILE=${2}

# SAMPLE_NAME=${3}

# SORTED_BAM=${4}


container="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

#conversion from sam to bam
apptainer exec ${container} samtools view -bo ${OUTPUT_FILE}.bam ${SAM_FILE} 


# sort bam files directly from sam file
apptainer exec ${container} samtools sort ${OUTPUT_FILE}.bam ${OUTPUT_FILE}.sorted.bam \
    -m ${SLURM_MEM_PER_CPU} -@ ${SLURM_CPUS_PER_TASK}

# index sorted bam
apptainer exec ${container} samtools index -@ ${SLURM_CPUS_PER_TASK} ${OUTPUT_FILE}.sorted.bam 

