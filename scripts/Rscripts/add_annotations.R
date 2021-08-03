#!/bin/bash

# Given annovar annotations, this script adds annotations in any given df. 

library(tidyverse)
library(stringr)
library(tidyr)

# FUNCTIONS
add_annotations <- function(path_to_curated_muts_df, dir_to_annovar, output_path) {
	# some row is taken from the df to annotate 
	# dir_to_annovar_df is the path to the directory that contains annovar results for samples

	
	# now go through the muts in the df - will include an apply function instead in the future
	df <- as.data.frame(read_csv(path_to_curated_muts_df))

	row_list <- list()
	i <- 1 
	while (i <= dim(df)[1]){ # go through each row 
		
		some_row <- df[i, ] 
		message(c(some_row$SAMPLE_ID, " ", some_row$CHROM, " ", some_row$POSITION))

		# filter annovar results based on which sample and variant we are at
		# read in annovar results
		annovar_df <- as.data.frame(read.delim(file.path(dir_to_annovar, paste0(some_row$SAMPLE_ID, ".bam_FILTERED_vcf.ANNOTATED.hg38_multianno.txt"))))
		annovar_df <- annovar_df %>% 
			filter(Chr == some_row$CHROM, Start == some_row$POSITION) %>% # filter to keep the variant of interest
			select(Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene, ExAC_ALL) # exonic or intronic, which gene, type of variation, which aa change, exac score

		some_row <- as.data.frame(cbind(as.data.frame(some_row), annovar_df)) # combine the called variant with annovar annotations
		row_list <- append(row_list, list(some_row)) # add the modified row to the list

		i <- i + 1

	} # end of while loop

	updated_df <- as.data.frame(do.call(rbind, row_list)) # updated df with the annotations
	write_csv(updated_df, output_path)

	return(updated_df)

}

path_to_curated_muts_df <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/curated_muts_panel.csv" # whatever file you are working with 
dir_to_annovar <- "/groups/wyattgrp/users/amunzur/chip_project/annovar/annovar_results/new_finland_download"
output_path <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/curated_muts_panel_annotated.csv"

add_annotations(path_to_curated_muts_df, dir_to_annovar)