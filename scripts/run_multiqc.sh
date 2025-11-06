#!/usr/bin/env bash

mkdir /results/multiQC
#SBATCH --job-name=Blood_sample
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/ankunzimana/RNAseq_project/logs/out/multiqc_%j.out
#SBATCH --error=/data/users/ankunzimana/RNAseq_project/logs/err/multiqc_%j.err
#SBATCH --time=02:10:00