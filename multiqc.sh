#!/bin/bash

#SBATCH --job-name=trimming
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=14:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mjw85@student.le.ac.uk
#SBATCH --account=evo-epi

### Include when testing
# SBATCH --partition=devel

### Run script in the working directory it was submitted in
cd $SLURM_SUBMIT_DIR

# Load modules
module load py-multiqc/1.14-ateq76h
module load fastqc/0.12.1-hkgpcde

## Trimming:
module load trimmomatic/0.39-tfj436w
for file in $(ls *_1.fastq.gz)
do
	base=$(basename ${file} "_1.fastq.gz")
	trimmomatic PE -threads 16 ${base}_1.fastq.gz ${base}_2.fastq.gz \
    ${base}_trim_1.fq.gz ${base}_unpaired_1.fq.gz \
	${base}_trim_2.fq.gz ${base}_unpaired_2.fq.gz \
    ILLUMINACLIP:illumina_adapters_new.fa:2:30:10 MINLEN:50 HEADCROP:10
done 

for file in $(ls *trim*.fq.gz)
do
	fastqc -t 24 ${file}
done

multiqc ./
