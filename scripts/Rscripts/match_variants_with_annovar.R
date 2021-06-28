library(tidyverse)

filter_variants_frequency <- function(r){

	# Given a row that has variant information, add a PASS or FAIL status based on AF. 
	
	idx <- which(unlist(r) > 0 & unlist(r) < 1) # identify all elements that are between 0 and 1, they are the allele frequencies 
	
	if (length(idx) == 0 | any(r[idx] > 0.05)) {
		return("BAD")
	} else if (all(r[idx] < 0.05)){
		return("GOOD")
	}
}

filter_variants_bed <- function(df, bed){

	# Given a bed file and a csv file that has variant information, filter the df to retain the regions in the bed file only. 

	i = 1
	bed_status <- list()
	while (i <= nrow(df)){
		
		# extract chrom and position from the df 
		chrom_df <- as.vector(df[i, ]$CHROM) 
		pos_df <- as.vector(df[i, ]$POS)

		bed_chrom <- bed %>% filter(chrom == chrom_df) # filter the bed file to retain the chrom of interest
		if (any(pos_df >= bed$start & pos_df <= bed$stop)) {bed_status <- append(bed_status, "keep")} else {bed_status <- append(bed_status, "remove")}

		i <- i + 1 
	}

	idx <- which(bed_status == "keep")
	df_filtered <- df[idx, ]

	return(df)
}

annovar_path <- "/groups/wyattgrp/users/amunzur/chip_project/annovar_results/ctDNA_prognosis/"
common_variants_path <- "/groups/wyattgrp/users/amunzur/chip_project/common_variants/"
bed_path <- "/groups/wyattgrp/users/echen/baits_hg38.bed"

# load the dfs with common variants
for (file_name in list.files(common_variants_path)) {

	sample_name <- strsplit(file_name, ".csv")[[1]][[1]]
	print(c("Started sample", sample_name))
	idx <- grep(paste0("^", sample_name, ".+_vcf_FILTERED_vcf.ANNOTATED.hg38_multianno.txt$"), 
	                  list.files(annovar_path))
	anno_name <- list.files(annovar_path)[idx]

	# load both files
	sample <- read.csv(paste0(common_variants_path, file_name))
	sample$X <- NULL
	anno <- as.data.frame(read.delim(paste0(annovar_path, anno_name), sep = "\t"))
	names(anno)[1:5] <- c("CHROM", "POS", "END", "REF", "ALT")
	
	# KEEP ONLY COMMON ONES
	combined <- merge(x = sample, y = anno, by = c("CHROM", "POS", "ALT", "REF")) # merge!
	combined[combined == "."] <- NA # replace all dots with NA, makes it easier to filter later on
	combined <- Filter(function(x) !all(is.na(x)), combined) # drop cols if they are all filled with NA

	# REMOVE GERM LINE BASED ON DATABASE INFORMATION
	status_list <- apply(combined, 1, filter_variants_frequency) # apply the function to every row to make of list of pass or fail status
	idx <- which(status_list == "GOOD")
	combined_filtered <- combined[idx, ]

	# FILTERING BASED ON BED FILE
	bed <- read.delim(bed_path, header = TRUE)
	combined_filtered <- filter_variants_bed(combined_filtered, bed)

	path <- paste0("/groups/wyattgrp/users/amunzur/chip_project/bed_variants/", sample_name, ".tsv")
	write_delim(combined_filtered, path, delim = "\t")

}