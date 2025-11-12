#!/usr/bin/env bash

input_dir="/data/users/ankunzimana/RNAseq_project/results/alignment/sam"

script_dir="/data/users/ankunzimana/RNAseq_project/scripts/alignment"

output_dir="/data/users/ankunzimana/RNAseq_project/results/alignment/bam"

for sam_file in `ls -1 ${input_dir}/*.sam`; do

    sample_name=$(basename ${sam_file%.sam})
    # echo "${sam_file} ${output_dir}/${sample_name}"
    sbatch ${script_dir}/samtools.slurm ${sam_file} "${output_dir}/${sample_name}"
done