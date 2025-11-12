#!/usr/env/bin bash

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

sam_file= ${1}

bam_file= ${2}

sorted_bam= ${3}

output_dir= "/data/users/ankunzimana/RNAseq_project/results/alignment/bam"

container="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

#conversion from sam to bam
apptainer exec ${container} samtools view -hbS ${sam_file} > ${bam_file}

# sort bam files 
apptainer exec ${container} samtools sort ${bam_file} -o ${sorted_bam} \
    -m ${SLURM_MEM_PER_CPU} -@ ${SLURM_CPUS_PER_TASK}

# index sorted bam
apptainer exec ${container} smatools index ${sorted_bam}

