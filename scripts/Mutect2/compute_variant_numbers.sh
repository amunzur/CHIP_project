#!/bin/bash

# To compute the number of variants that passed filtering in each sample
path_to_mutect_FILTERED="/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/new_finland_download"
path_to_metrics="/groups/wyattgrp/users/amunzur/chip_project/metrics/mutect_variant_numbers/new_finland_download_VARIANTS"

cd $path_to_mutect_FILTERED
touch $path_to_metrics

for vcf in *.bam_FILTERED_vcf
do
	variant_number=$(grep -v "#" $vcf | wc -l)
	echo $vcf $variant_number
	echo $vcf $variant_number >> $path_to_metrics
done

# To compute the number of common variants in samples based on merged csv files
path_to_common_variants="/groups/wyattgrp/users/amunzur/chip_project/common_variants/new_finland_download"
path_to_metrics="/groups/wyattgrp/users/amunzur/chip_project/metrics/mutect_variant_numbers/new_finland_download_COMMON"
cd $path_to_common_variants

for csv in *
do 
	variant_number=$(wc -l $csv | awk '{print $1-1}') # -1 because the first line is the header
	echo $csv $variant_number
	echo $csv $variant_number >> $path_to_metrics
done

