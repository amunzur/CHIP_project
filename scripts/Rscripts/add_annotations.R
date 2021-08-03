#!/bin/bash

# Given annovar annotations, this script adds annotations in any given df. 

library(tidyverse)
library(stringr)
library(tidyr)

# FUNCTIONS
add_annotations <- function(some_row) {
	# some row is taken from the df to annotate 
	# DIR_annovar is the path to the directory that contains annovar results for samples

	DIR_annovar_df <- "/groups/wyattgrp/users/amunzur/chip_project/annovar/annovar_results/new_finland_download" # annovar annotation results 

	some_row <- as.list(some_row)
	message(c(some_row$SAMPLE_ID, " ", some_row$CHROM, " ", some_row$POSITION))

	annovar_df <- as.data.frame(read.delim(file.path(DIR_annovar_df, paste0(some_row$SAMPLE_ID, ".bam_FILTERED_vcf.ANNOTATED.hg38_multianno.txt"))))
	annovar_df <- annovar_df %>% 
		filter(Chr == some_row$CHROM, Start == some_row$POSITION) %>% # filter to keep the variant of interest
		select(Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene, ExAC_ALL) # exonic or intronic, which gene, type of variation, which aa change, exac score

	print(dim(some_row))
	print(annovar_df)


	some_row <- as.data.frame(cbind(some_row, annovar_df))

	return(some_row)

}

# DECLARE VARIABLES
bam_file <- "GU-19-331-WBC-2019Aug14.bam"
DIR_df <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/curated_muts_panel.csv" # whatever file you are working with 

replace <- FALSE # replace the original file with a new one with annotations

snv_list <- apply(df, 1, add_annotations) # saves annotations in a list, which we will make a df out of
snv_df <- as.data.frame(do.call(rbind, snv_list))


snv_df %>% separate(V1, c(CHROM, POSITION, SAMPLE_ID, Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene, ExAC_ALL), " ")


str_split_fixed(snv_df$V1, "\t", 9)
