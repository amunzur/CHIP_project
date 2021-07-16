#!/bin/bash
configfile="/groups/wyattgrp/users/amunzur/chip_project/scripts/Mutect2/config.txt"
conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh"
source ${configfile};
source ${conda_profile_path};
conda activate mito 

###########################################
# option 1: run the code on all bam files in the dir given below
###########################################
# cd /groups/wyattgrp/users/amunzur/chip_project/finland_bams/GU_finland_download_ORIGINAL

# for bam_file in GU-*.bam 
# do
# 	printf "bash ${script_dir}/run_Mutect2.sh ${bam_file}"
# 	sbatch --exclude=cn[01-05] /groups/wyattgrp/users/amunzur/chip_project/scripts/run_Mutect2.sh ${bam_file};
# done

###########################################
# option 2: run the code on the bam files given in a text file, ids only, not absolute paths.
###########################################
# or can also give a list of paths to bam files to run, instead of running it on every bam file in the given dir.
path_to_bamList="/groups/wyattgrp/users/amunzur/chip_project/scripts/Mutect2/bamList.txt"
cd /groups/wyattgrp/users/amunzur/chip_project/finland_bams/new_finland_download # I'll fix this later but for now, we need to cd into the dir with bams, even if we use a file with bam ids.

cat $path_to_bamList | while read bam_file || [[ -n $bam_file ]];
do
	printf "bash ${script_dir}/run_Mutect2.sh ${bam_file} "
	sbatch --exclude=cn[01-05] ${script_dir}/run_Mutect2.sh ${bam_file};
done

