#!/usr/bin/env bash

# setting pathways

input_dir="/data/users/ankunzimana/RNAseq_project/results/alignment/sam"

script_dir="/data/users/ankunzimana/RNAseq_project/scripts/alignment"

output_dir="/data/users/ankunzimana/RNAseq_project/results/alignment/bam"

#apply samtools.slurm on every sam file 
for sam_file in `ls -1 ${input_dir}/*.sam`; do

    sample_name=$(basename ${sam_file%.sam})

    bam_file="${output_dir}/${sample_name}.bam"

    touch ${output_dir}/${sample_name}.sorted.bam

    sorted_bam="${output_dir}/${sample_name}.sorted.bam"
    # echo "${sam_file} ${bam_file} ${sorted_bam}"
    sbatch ${script_dir}/samtools.slurm ${sam_file} ${bam_file} ${sorted_bam}
done