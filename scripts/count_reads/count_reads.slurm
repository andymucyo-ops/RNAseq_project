#!/usr/bin/env bash 

#SBATCH --job-name=samtools_RNAseq
#SBATCH --mem=1G
#SBATCH --cpus-per-task=4
#SBATCH --partition=pibu_el8
#SBATCH --output=./logs/out/counts_%j.out
#SBATCH --error=./logs/err/counts_%j.err
#SBATCH --time=02:00:00
#SBATCH --mail-user= <your email>
#SBATCH --mail-type=end,fail

#set work directories
input_files="./results/alignment/bam/*sorted.bam"
annotation_file="./data/reference_genome/Mus_musculus.GRCm39.115.gtf"
container="/containers/apptainer/subread_2.0.6.sif"
output_dir="./results/count_reads.txt"

apptainer exec ${container} featureCounts -p --countReadPairs\
    -a ${annotation_file}\
    -o ${output_dir}/counts.txt\
    ${input_files}\
    -g gene_id -t exon\
    -T $SLURM_CPUS_PER_TASK \
    -s 2

sed -E -i 's#[^[:space:]]*/(SRR[0-9]+)\.sorted\.bam#\1#g' ${output_dir}/counts.txt
tail -n +2 ${output_dir}/counts.txt | cut -f 1,7-21 > ${output_dir}"/counts.txt"
sed -E -i 's#[^[:space:]]*/(SRR[0-9]+)\.sorted\.bam#\1#g' ${output_dir}/counts.txt.summary


