#!/bin/bash
configfile="/groups/wyattgrp/users/amunzur/chip_project/scripts/Mutect2/config.txt"
conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh"
source ${configfile};
source ${conda_profile_path};
conda activate mito 

type="WBC"

if [[ "$type" == "ctDNA" ]]
then
		cd /groups/wyattgrp/users/amunzur/chip_project/mutect2_results/ctDNA_prognosis
		printf "TUMOR DNA\n"
else
		cd /groups/wyattgrp/users/amunzur/chip_project/mutect2_results/GU_finland_download
		printf "WBC DNA\n"
fi

###########################################
# option 1: run the code on all vcf files in the dir given below
###########################################
for vcf_file in GU-*.bam_vcf 
do
	printf "bash ${script_dir}/filter_Mutect2.sh ${vcf_file} ${type} "
	sbatch --exclude=cn[01-05] /groups/wyattgrp/users/amunzur/chip_project/scripts/Mutect2/filter_Mutect2.sh ${vcf_file} ${type};
done

###########################################
# option 2: run the code on some samples only
# This option uses the same list of samples as the run_Mutect2.sh script.
###########################################
