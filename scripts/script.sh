#!/bin/bash
bam_file=$1
output_dir="/groups/wyattgrp/users/amunzur/chip_project/mutect2_results"
conda activate mito

printf "\n"
printf "*******************************\n"
printf "FIXMATE - $bam_file\n"
printf "*******************************\n"

picard MarkDuplicates I=bam_file O='rmdup_$bam_file' M='/groups/wyattgrp/users/amunzur/chip_project/finland_bams/MarkDuplication_metrics/rmdup_$bam_file' --REMOVE_DUPLICATES=true

# picard FixMateInformation I='rmdup_$bam_file' O='fixed_mate_$bam_file' ADD_MATE_CIGAR=true

# picard AddOrReplaceReadGroups \
#    I="fixed_mate_$bam_file" \
#    O="RG_$bam_file" \
#    RGID=1 \
#    RGLB=library \
#    RGPL=ILLUMINA \
#    RGPU=unit \
#    RGSM=sample

# samtools index "RG_$bam_file"

# printf "\n"
# printf "******************************\n*"
# printf "VARIANT CALLING\n"
# printf "*******************************\n"

# /home/amunzur/gatk-4.2.0.0/gatk Mutect2 \
# -R /groups/wyattgrp/reference/hg38/hg38.fa \
# -I "RG_$bam_file" \
# -O "$output_dir$bam_vcf.gz"
