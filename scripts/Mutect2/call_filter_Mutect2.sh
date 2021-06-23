#!/bin/bash
configfile="/groups/wyattgrp/users/amunzur/chip_project/scripts/Mutect2/config.txt"
conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh"
source ${configfile};
source ${conda_profile_path};
conda activate mito 

type=$1

if [[ "$type" == "ctDNA" ]]
then
		cd /groups/wyattgrp/users/amunzur/chip_project/mutect2_results/ctDNA_prognosis
		printf "TUMOR DNA\n"
else
		cd /groups/wyattgrp/users/amunzur/chip_project/mutect2_results/GU_finland_download
		printf "WBC DNA\n"
fi

for vcf_file in GU-*.bam_vcf 
do
	printf "bash ${script_dir}/filter_Mutect2.sh ${vcf_file} ${type} "
	sbatch --exclude=cn[01-05] /groups/wyattgrp/users/amunzur/chip_project/scripts/Mutect2/filter_Mutect2.sh ${vcf_file} ${type};
done
