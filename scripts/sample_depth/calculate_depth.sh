#!/bin/bash 
# This script will calculate coverage for a list of bam files, and should be run with GNU parallel to minimize the time. 
# run this in the dir with the bams of interest:
# cd /groups/wyattgrp/users/amunzur/chip_project/finland_bams/ctDNA_prognosis
# ls -d -1 "$PWD/"filtered_RG_*.bam | parallel -j10 -k bash /groups/wyattgrp/users/amunzur/chip_project/scripts/sample_depth/calculate_depth.sh

# ca variant_calling
bam=$1 # supplied by GNU parallel in the pipe
path_to_bed="/groups/wyattgrp/users/amunzur/chip_project/references/baits_hg38.bed" 
dir_to_read_depth="/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/" # the path where read depth information will be saved, OUTPUT
# output_path="/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/GU_finland_download_coverage.txt"

output_file_txt="${dir_to_read_depth}new_finland_download_FILTERED.txt"
output_file_csv="${dir_to_read_depth}new_finland_download_FILTERED.csv"

touch $output_file_txt # initialize an empty file to append coverage information

read_num=$(samtools view $bam | wc -l)
if (( $read_num > 0 )); then # if we have any reads
	
	depth=`samtools depth -b $path_to_bed $bam |  awk '{sum+=$3} END { print sum/NR}'`
	echo $bam $depth >> $output_file_txt
	echo $bam
fi

# convert from txt to csv
# cat $output_file_txt | tr -s '[:blank:]' ',' > $output_file_csv 
# rm $output_file_txt
