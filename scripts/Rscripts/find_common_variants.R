library(dplyr)

return_sample <- function(ids_list, str_separator){

	# Extract the sample id from a given list of ids. ids list contains the full name of the bam files, here we extract the sample name only and disregard the rest of the strings. 
	sample_names <- lapply(ids_list, function(x) strsplit(x, str_separator)[[1]][[1]])

	return(sample_names)

}

subset_list <- function(ids_list, samples_list, str_separator){
	# Given a list of strings, subset it based on another list of strings. 

	# subset both to common elements 
  ids_list_stripped <- lapply(ids_list, function(x) strsplit(x, str_separator)[[1]][[1]])
  idx <- match(intersect(ids_list_stripped, samples_list), ids_list_stripped)
	ids_list <- ids_list[idx]

	return(ids_list)

}

identify_vcf_files <- function(sample_type, sample_list){
	# sample type: given as a string, either "cfDNA" or "WBC". 
	# sample_list: Sample ids we are interested in retrieving as a vcf file, since not all files have a tumor and wbc match. 

	# load vcf based on the sample type provided by user
	if (sample_type == "cfDNA") {path_to_vcf <- paste("/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/ctDNA_prognosis")}
	if (sample_type == "WBC") {path_to_vcf <- paste("/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/GU_finland_download")}
	
	vcf_list <- list() # intialize an empty list for the loop below 
	for (sample_name in samples_list) {

		# use regex to catch files starting witht the sample id and ending with the appropriate suffix
		idx <- grep(paste0("^", sample_name, ".+bam_vcf_FILTERED_vcf$"), list.files(path_to_vcf), ignore.case = TRUE) # idx of the file we are interested in based on grep results, we need to grep once more though to remove the stats files we dont want
		vcf_file <- list.files(path_to_vcf)[idx]
		
		vcf_list <- append(vcf_list, vcf_file) # add the file to the list 

	} # end of for loop

	return(as.list(sort(unlist(vcf_list))))
}

find_common_variants <- function(tumor_vcf, wbc_vcf)
	# this function assumes that these files are the same length, and they have the same ids in the same location in both

# get the sample ids and sort 
tumor_ids <- as.list(sort(grep("^[GU]", list.files("/groups/wyattgrp/users/amunzur/chip_project/finland_bams/ctDNA_prognosis_ORIGINAL"), value = TRUE)))
wbc_ids <- as.list(sort(grep("^[GU]", list.files("/groups/wyattgrp/users/amunzur/chip_project/finland_bams/GU_finland_download_ORIGINAL"), value = TRUE)))

# subset
tumor_sample <- return_sample(tumor_ids, "-cfDNA")
wbc_sample <- return_sample(wbc_ids, "-WBC")

if (identical(tumor_sample, wbc_sample)) {print("Good to go!")} else {
  print("Something is seriously wrong, working on it.")
  
  samples_list <- intersect(tumor_sample, wbc_sample)
  print("All good now.")
    }

tumor_ids <- subset_list(tumor_ids, samples_list, "-cfDNA")
wbc_ids <- subset_list(wbc_ids, samples_list, "-WBC")

# identify the vcf files
tumor_vcf_paths <- as.list(identify_vcf_files("cfDNA", samples_list))
wbc_vcf_paths <- as.list(identify_vcf_files("WBC", samples_list))

# load vcf files 
tumor_vcf <- 



