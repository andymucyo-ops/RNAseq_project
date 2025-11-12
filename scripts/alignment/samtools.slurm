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

BAM_FILE=${2}

SORTED_BAM=${3}



container="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

#conversion from sam to bam
apptainer exec ${container} samtools view -hbo ${BAM_FILE} ${SAM_FILE}
echo "$(basename ${SAM_FILE}) converted to bam"

# sort bam files directly from sam file
apptainer exec ${container} samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${SORTED_BAM} -T temp ${BAM_FILE} 
echo "$(basename ${BAM_FILE}) has been sorted"

# index sorted bam
apptainer exec ${container} samtools index ${SORTED_BAM} 
echo "$(basename ${SORTED_BAM}) has been indexed"
