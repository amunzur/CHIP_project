#!/bin/bash

#SBATCH --job-name=Mutect2
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH --cpus-per-task=31
#SBATCH --mem 200000 # memory pool for all cores
#SBATCH -t 08:30:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log

vcf_file=$1
type=$2

if [[ "$type" == "ctDNA" ]]
then
		output_dir="/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/ctDNA_prognosis/"
else
		output_dir="/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/GU_finland_download/"
fi

printf "\n"
printf "*******************************\n"
printf "FILTER MUTECT CALLS - $bam_file\n"
printf "*******************************\n"

/home/amunzur/gatk-4.2.0.0/gatk FilterMutectCalls \
-R /groups/wyattgrp/users/amunzur/chip_project/references/hg38.fa \
-V $vcf_file \
-O "${output_dir}${vcf_file}_FILTERED_vcf"

