#!/bin/bash
configfile="/groups/wyattgrp/users/amunzur/chip_project/scripts/Mutect2/config.txt"
conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh"
source ${configfile};
source ${conda_profile_path};

type=$1

if [[ "$type" == "ctDNA" ]]
then
		cd /groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/ctDNA_prognosis

else
		cd /groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/GU_finland_download
fi

for vcf_file in GU-*.bam_vcf_FILTERED_vcf
do
        printf "bash /groups/wyattgrp/users/amunzur/chip_project/scripts/annovar/run_annovar.sh ${vcf_file} ${type} "
        sbatch --exclude=cn[01-05] /groups/wyattgrp/users/amunzur/chip_project/scripts/annovar/run_annovar.sh ${vcf_file} ${type};
done

