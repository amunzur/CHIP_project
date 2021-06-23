#!/bin/bash

#SBATCH --job-name=Mutect2
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH --cpus-per-task=31
#SBATCH --mem 200000 # memory pool for all cores
#SBATCH -t 08:30:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log

bam_file=$1
output_dir="/groups/wyattgrp/users/amunzur/chip_project/mutect2_results/GU_finland_download"

printf "\n"
printf "*******************************\n"
printf "FIXMATE - $bam_file\n"
printf "*******************************\n"

picard MarkDuplicates I=$bam_file O=rmdup_$bam_file M=/groups/wyattgrp/users/amunzur/chip_project/finland_bams/MarkDuplication_metrics/rmdup_$bam_file REMOVE_DUPLICATES=true

picard FixMateInformation I=rmdup_$bam_file O=fixed_mate_$bam_file ADD_MATE_CIGAR=true

picard AddOrReplaceReadGroups \
    I=fixed_mate_$bam_file \
    O=RG_$bam_file \
    RGID=1 \
    RGLB=library \
    RGPL=ILLUMINA \
    RGPU=unit \
    RGSM=sample

samtools index RG_$bam_file

printf "\n"
printf "******************************\n*"
printf "VARIANT CALLING\n"
printf "*******************************\n"

/home/amunzur/gatk-4.2.0.0/gatk Mutect2 \
-R /groups/wyattgrp/users/amunzur/chip_project/references/hg38.fa \
-I $RG_bam_file \
-O $output_dir$bam_vcf.gz

