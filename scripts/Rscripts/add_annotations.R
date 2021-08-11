#!/bin/bash

# Given annovar annotations, this script adds annotations in any given df. 

library(tidyverse)
library(stringr)
library(tidyr)

cohort <- "new_finland_download"
germline_threshold <- 0.05 # only keep variants if they are below this threshold

# FUNCTIONS
add_annotations <- function(path_to_muts_df, dir_to_annovar, germline_threshold, output_path, output_path_germline) {
	# path_to_muts_df can be the curated df based on snapshots, or any other df with mutation information 
	# some row is taken from the df to annotate 
	# dir_to_annovar_df is the path to the directory that contains annovar results for samples

	df <- as.data.frame(read_csv(path_to_muts_df)) 	# now go through the muts in the df - will include an apply function instead in the future
	# this bit helps accomodate different methods of naming 
	if (length(which(names(df) == "SAMPLE_ID")) == 0) {names(df)[which(names(df) == "sample_t")] <- "SAMPLE_ID"}
	if (length(which(names(df) == "POSITION")) == 0) {names(df)[which(names(df) == "POS")] <- "POSITION"}

	row_list <- list()
	i <- 1 
	while (i <= dim(df)[1]){ # go through each row 
		
		some_row <- df[i, ] 
		message(c(some_row$SAMPLE_ID, " ", some_row$CHROM, " ", some_row$POSITION))

		# filter annovar results based on which sample and variant we are at
		if (cohort == "first_batch") {annovar_suffix <- ".bam_vcf.ANNOTATED.hg38_multianno.txt"} else {".bam_FILTERED_vcf.ANNOTATED.hg38_multianno.txt"}
		annovar_df <- as.data.frame(read.delim(file.path(dir_to_annovar, paste0(some_row$SAMPLE_ID, annovar_suffix)))) # read annovar results 
		annovar_df <- annovar_df %>% 
			filter(Chr == some_row$CHROM, Start == some_row$POSITION) %>% # filter to keep the variant of interest
			select(Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene, ExAC_ALL) # exonic or intronic, which gene, type of variation, which aa change, exac score

		some_row <- as.data.frame(cbind(as.data.frame(some_row), annovar_df)) # combine the called variant with annovar annotations
		row_list <- append(row_list, list(some_row)) # add the modified row to the list

		i <- i + 1

	} # end of while loop

	updated_df <- as.data.frame(do.call(rbind, row_list)) # updated df with the annotations

	# germline filtering
	updated_df$ExAC_ALL[grep("^\\.$", updated_df$ExAC_ALL)] <- 0 # if no exac_all annots, add 0 to retain them
	updated_df_germline <- filter(updated_df, as.numeric(as.vector(updated_df$ExAC_ALL)) < germline_threshold) 

	write_csv(updated_df, output_path)
	write_csv(updated_df_germline, output_path_germline)

	return(updated_df)

}

path_to_muts_df <- file.path("/groups/wyattgrp/users/amunzur/chip_project/subsetted", cohort, "curated_muts_panel.csv") # input
dir_to_annovar <- file.path("/groups/wyattgrp/users/amunzur/chip_project/annovar/annovar_results", cohort)
output_path <- file.path("/groups/wyattgrp/users/amunzur/chip_project/subsetted", cohort, "curated_muts_panel_annotated.csv")
output_path_germline <- file.path("/groups/wyattgrp/users/amunzur/chip_project/subsetted", cohort, "curated_muts_panel_annotated_germFiltered.csv")

updated_df <- add_annotations(path_to_muts_df, dir_to_annovar, germline_threshold, output_path, output_path_germline)