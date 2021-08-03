bam_file=$1

output_dir="/groups/wyattgrp/users/amunzur/chip_project/mutect2_results/kidney_samples/"
output_dir_filtered="/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/kidney_samples/"

echo $bam_file

/home/amunzur/gatk-4.2.0.0/gatk FilterMutectCalls \
-R /groups/wyattgrp/users/amunzur/chip_project/references/hg38.fa \
-V "${output_dir}${bam_file}_vcf.gz" \
-O "${output_dir_filtered}${bam_file}_FILTERED_vcf"


# cd /groups/wyattgrp/users/amunzur/chip_project/finland_bams/kidney_samples # I'll fix this later but for now, we need to cd into the dir with bams, even if we use a file with bam ids.
# ls GU*WBC*.bam | parallel -j8 -k bash /groups/wyattgrp/users/amunzur/chip_project/scripts/Mutect2/temp/filter_Mutect_kidney_WBC.sh
