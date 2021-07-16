#!/bin/bash
# ls GU*.bam | parallel -j6 -k bash ../process_samples.sh

sample=$1
if [ ! -f "NEW_filtered_${sample}" ]
then
	printf "Started sample ${sample}.\n"
	samtools sort -n -o "name_sorted_${sample}" $sample # name sort for fixmate
	samtools fixmate -m -O bam "name_sorted_${sample}" "fixmate_${sample}" # run fixmate and save output
	samtools sort -@4 -o "sorted_${sample}" "fixmate_${sample}" # sort by coordinate
	samtools markdup -r "sorted_${sample}" "rmdup_${sample}" # mark and remove duplicates

	samtools view -bq 20 -F 4 "rmdup_${sample}" > "NEW_filtered_${sample}" # remove unmapped and reads with MQ less than 20 

else
	echo "SAMPLE ALREADY PROCESSED."

fi