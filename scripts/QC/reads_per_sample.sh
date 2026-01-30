#!/usr/bin/env bash

input_dir="./data/reads_Blood/*fastq.gz"

output_file="./results/QC_data/QC_info/reads_per_sample.txt" 

for file in `ls -1 ${input_dir}`; do 

    # read_count= zcat ${file} | grep '@' | wc -l

    echo "$(zcat ${file} | grep '@' | wc -l) in $(basename ${file%_fastq.gz})" >> ${output_file}
done
