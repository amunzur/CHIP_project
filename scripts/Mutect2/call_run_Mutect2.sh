#!/bin/bash
configfile="/groups/wyattgrp/users/amunzur/chip_project/scripts/Mutect2/config.txt"
conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh"
source ${configfile};
source ${conda_profile_path};
conda activate mito 

cd /groups/wyattgrp/users/amunzur/chip_project/finland_bams/GU_finland_download_ORIGINAL

for bam_file in GU-*.bam 
do
	printf "bash ${script_dir}/run_Mutect2.sh ${bam_file}"
	sbatch --exclude=cn[01-05] /groups/wyattgrp/users/amunzur/chip_project/scripts/run_Mutect2.sh ${bam_file};
done
