#!/usr/bin/env bash

input_dir="/data/users/ankunzimana/RNAseq_project/results/alignment"

script_dir="/data/users/ankunzimana/RNAseq_project/scripts/alignment"

for sam_file in `ls -1 ${input_dir}/sam/*.sam`; do
    bam_file= "${input_dir}/bam/$(basename ${sam_file%.sam}).bam"
    sorted_bam= "${input_dir}/bam/$(basename ${sam_file%.sam}).sorted.bam"

    sbatch ${script_dir}/samtools.slurm ${sam_file} ${bam_file} ${sorted_bam}
done