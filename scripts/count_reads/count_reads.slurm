#!/usr/bin/env bash 

#SBATCH --job-name=samtools_RNAseq
#SBATCH --mem=1G
#SBATCH --cpus-per-task=4
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/ankunzimana/RNAseq_project/logs/out/counts_%j.out
#SBATCH --error=/data/users/ankunzimana/RNAseq_project/logs/err/counts_%j.err
#SBATCH --time=02:00:00
#SBATCH --mail-user=andy.nkunzimana@students.unibe.ch
#SBATCH --mail-type=end,fail

#set work directories
input_files="/data/users/ankunzimana/RNAseq_project/results/alignment/bam/*sorted.bam"
annotation_file="/data/users/ankunzimana/RNAseq_project/data/reference_genome/Mus_musculus.GRCm39.115.gtf"
container="/containers/apptainer/subread_2.0.6.sif"
output_dir="/data/users/ankunzimana/RNAseq_project/results/count_reads"

echo "${input_dir}"

apptainer exec ${container} featureCounts -p --countReadPairs\
    -a ${annotation_file}\
    -o ${output_dir}/counts.txt\
    ${input_files}\
    -g gene_id -t exon\
    -T $SLURM_CPUS_PER_TASK \
    -s 2

sed -E -i 's#[^[:space:]]*/(SRR[0-9]+)\.sorted\.bam#\1#g'



