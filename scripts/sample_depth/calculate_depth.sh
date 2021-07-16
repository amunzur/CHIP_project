#!/bin/bash 
# This script will calculate coverage for a list of bam files, and should be run with GNU parallel to minimize the time. You can either supply a text file where 
# each line contains the absolute path to the bam files, or give a path to the directory of bams. The script assumes that the list of files are piped into it. 

# ls -d -1 "$PWD/"NEW_filtered*.bam | parallel -j6 -k bash /groups/wyattgrp/users/amunzur/chip_project/scripts/sample_depth/calculate_depth.sh

bam=$1
dir_to_read_depth="/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/" # the path where read depth information will be saved, OUTPUT
# output_path="/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/GU_finland_download_coverage.txt"

temp1="${dir_to_read_depth}temp1" # the file with sample names
temp2="${dir_to_read_depth}temp2" # the file with depth information, will be concatted with temp1 later on
combined="${dir_to_read_depth}ctDNA_prognosis_coverage.txt"

cd $dir_to_read_depth # dir where the text files with depth information are saved
rm $temp1 $temp2
touch $temp1 $temp2 # initialize an empty text file in the read depth directory

read_num=$(samtools view $bam | wc -l)
if (( $read_num > 0 )); then # if we have any reads
    echo $bam # print the sample name to terminal 
    echo $bam >> $temp1; # sample name
    samtools depth $bam |  awk '{sum+=$3} END { print sum/NR}' >> $temp2; # average read depth from one sample
fi

a=`wc -l $temp1`
b=`wc -l $temp2`

if [[ $a -eq $b ]]; then # if both have the same number of lines
    paste $temp1 $temp2 > $combined;
fi