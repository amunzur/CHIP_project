#!/bin/bash 
# ls filtered_RG_*.bam | parallel -j10 -k bash /groups/wyattgrp/users/amunzur/chip_project/scripts/sample_depth/calculate_depth_individual.sh

bam=$1 # given by GNU parallel
dir_to_read_depth="/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/individual_samples/new_finland_download/" # the dir where read depth information will be saved, end with /
dir_to_bams="/groups/wyattgrp/users/amunzur/chip_project/finland_bams/new_finland_download/"
path_to_bed="/groups/wyattgrp/users/amunzur/chip_project/references/baits_hg38.bed"

# note that the read depth was calculated on samples after removing reads with low MQ and too short reads
cd $dir_to_bams # dir with all bams

path_to_metrics="${dir_to_read_depth}${bam}" # we keep the "bam" extension
echo $bam
samtools depth -b $path_to_bed $bam > $path_to_metrics

