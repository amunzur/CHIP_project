#!/bin/bash
# cd /groups/wyattgrp/users/amunzur/chip_project/finland_bams/kidney_samples
# ls GU*cfDNA*.bam | parallel -j10 -k bash /groups/wyattgrp/users/amunzur/chip_project/scripts/fix_finland_bams/get_read_names.sh

# This script, run with GNU Parallel, saves the read names from a bam file to a text file. 
bam=$1
path_readnames_all="/groups/wyattgrp/users/amunzur/chip_project/metrics/read_names/read_names_all/readNumbers_$bam"

echo $bam
samtools view $bam | cut -f 1 | sort > path_readnames_all