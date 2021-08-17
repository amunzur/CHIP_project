ptm <- proc.time() # start timer
noquote(format(Sys.time(), "%a %b %d %X %Y"))

# make the log file and create a file connection to write to it.
# log_file_name <- paste0(length(list.files("/groups/wyattgrp/users/amunzur/log/")), ".log")
# file.create(log_file_name)
# log_con <- file(log_file_name, open = "a") # open in the append more so that we can add stuff cumulatively

library(tidyverse)
library(vcfR)

# the path where the common variants will be saved as csv
path_to_common_variants <- "/groups/wyattgrp/users/amunzur/chip_project/common_variants/kidney_samples/" # make sure this ends with /

return_sample <- function(ids_list, str_separator){

	# Extract the sample id from a given list of ids. ids list contains the full name of the bam files, here we extract the sample name only and disregard the rest of the strings. 
	# This is the only way to match the tumor and wbc samples.
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

get_sample_ids <- function(sample_keyword, sample_type) {

	# Super simple function that calls grep to get the sample ids, NOT PATHS. Having a function for it made it more flexible for different cohorts.
	# These ids are sorted so that control and tumor samples match. Later on they are used to load the relevant vcf files.

	# sample_keyword is the dir name that contains the bams. options are: 
	# ctDNA_prognosis_ORIGINAL
	# GU_finland_download
	# new_finland_download
	# kidney_samples

	# sample_type is either ctDNA or WBC.
	
	main_path_to_bams <- "/groups/wyattgrp/users/amunzur/chip_project/finland_bams/" # path that stays the same regardless of the samples we work with
	path_to_bams <- file.path(main_path_to_bams, sample_keyword)

	# check later to see if these file paths are correct
	if (sample_keyword == "ctDNA_prognosis_ORIGINAL" | sample_keyword == "GU_finland_download") { 
		ids <- as.list(sort(grep("^[GU]", list.files(path_to_bams, pattern = ".bam$"), value = TRUE, ignore.case = TRUE))) 
	
	} else if (sample_type == "ctDNA") { 
		ids <- as.list(sort(grep("Baseline|cfDNA", list.files(path_to_bams, pattern = ".bam$"), value = TRUE, ignore.case = TRUE)))
	
	} else if (sample_type == "WBC") {
		ids <- as.list(sort(grep("WBC", list.files(path_to_bams, pattern = ".bam$"), value = TRUE, ignore.case = TRUE)))

	} # end of if loop

	# some further filtering based on sample_keyword
	if (sample_keyword == "new_finland_download") {ids <- ids[!grepl("filtered_RG_", ids)] # if we are dealing with new_finland_download, remove the duplicated names
	} else if (sample_keyword == "kidney_samples") {ids <- ids[grepl("^[GU]", ids)]}  

	return(ids)

} # end of function

identify_vcf_files <- function(sample_type, samples_list){
	# Given the sample type and a list of samples, load the relevant vcf files. 

	# Sample type: given as a string, choose one from the following. This helps locate the dir that has the relevant vcf files. 
	# "cfDNA" -> initial cohort, tumor samples
	# "WBC" -> initial cohort, wbc samples
	# "new_samples" -> second group of samples, wbc and tumor in the same dir
	# samples_list: Sample ids we are interested in retrieving as a vcf file, since not all files have a tumor and wbc match. 

	# load vcf based on the sample type provided by user
	if (sample_type == "new_samples") {path_to_vcf <- paste("/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/new_finland_download")}
	if (sample_type == "cfDNA") {path_to_vcf <- paste("/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/ctDNA_prognosis")}
	if (sample_type == "WBC") {path_to_vcf <- paste("/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/GU_finland_download")}
	if (sample_type == "kidney_samples") {path_to_vcf <- paste("/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/kidney_samples")}
	
	vcf_list <- list() # intialize an empty list for the loop below 
	for (sample_name in samples_list) {

		# use regex to catch files starting witht the sample id and ending with the appropriate suffix
		if (sample_type != "new_samples" & sample_type != "kidney_samples") {vcf_file <- file.path(path_to_vcf, paste0(sample_name, "_vcf_FILTERED_vcf")) # initial cohort
		} else { vcf_file <- file.path(path_to_vcf, paste0(sample_name, "_FILTERED_vcf")) } # second cohort

		if (file.exists(vcf_file)) { vcf_list <- append(vcf_list, vcf_file) # add the file to the list if it exists
		} else { print(c("File doesn't exist:", vcf_file)) }

	} # end of for loop

	return(as.list(sort(unlist(vcf_list))))
}

find_common_variants <- function(samples_list, tumor_vcf_paths, wbc_vcf_paths, path_to_common, save){

	i <- 1
	while (i <= length(samples_list)){

		sample <- samples_list[i]
		tumor_path <- tumor_vcf_paths[grep(sample, tumor_vcf_paths)]
		wbc_path <- wbc_vcf_paths[grep(sample, wbc_vcf_paths)]

		if (length(tumor_path) > 0 & length(wbc_path) > 0){
			
			message(c("Found both files, started sample ", unlist(sample)))
			# these are actually lists, we will extract the actual path if we can identify both the ctDNA and WBC sample
			tumor_path <- tumor_path[[1]]
			wbc_path <- wbc_path[[1]]

			# load the vcf files 
			message("Started loading files.")
			tumor <- read.vcfR(tumor_path, verbose = FALSE)
			wbc <- read.vcfR(wbc_path, verbose = FALSE)

			# only get what we need from the files
			tumor <- as.data.frame(tumor@fix)[, c("CHROM",  "POS", "REF", "ALT", "FILTER")]
			wbc <- as.data.frame(wbc@fix)[, c("CHROM",  "POS", "REF", "ALT", "FILTER")]

			# remove the variants that failed
			message("Filtering.")
			tumor <- filter(tumor, FILTER == "PASS")
			wbc <- filter(wbc, FILTER == "PASS")

			# now merge to keep the variants if they exist in both files
			message("Merging.")
			combined <- merge(tumor, wbc)
			path_combined <- paste0(path_to_common, sample, ".csv")
			
			if (dim(combined)[1] != 0) {
				message(c("Writing to csv. Found ", dim(combined)[1], " common variants."))
				#cat(c("Writing to csv. Found", dim(combined)[1], "common variants."), , file = log_con)
				if (save == TRUE) {
					write.csv(combined, path_combined)
					message(c("Saved to ", path_combined))
					} else {"No common variants."}}

			message("\n")

		} else {
			message(c("Skipped sample: ", sample))

		}
		i = i + 1

	} # end of while loop

} # end of function 
 
####################################
# FOR GU COHORT, OUR INITIAL SAMPLES:
####################################

# tumor_ids <- get_sample_ids("cfDNA")
# wbc_ids <- get_sample_ids("WBC")

# tumor_sample <- return_sample(tumor_ids, "-cfDNA") # subset to make sure we have the same samples in both groups
# wbc_sample <- return_sample(wbc_ids, "-WBC")

# tumor_vcf_paths <- as.list(identify_vcf_files("cfDNA", samples_list))
# wbc_vcf_paths <- as.list(identify_vcf_files("WBC", samples_list))

# add the prefixes to both to include the absolute path 
# tumor_vcf_paths <- lapply(tumor_vcf_paths, function(some_string) paste0("/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/ctDNA_prognosis/", some_string))
# wbc_vcf_paths <- lapply(wbc_vcf_paths, function(some_string) paste0("/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/GU_finland_download/", some_string))

####################################
# FOR THE NEW BATCH OF SAMPLES FROM FINLAND - new_finland_download
####################################
# tumor_ids <- get_sample_ids("new_samples_cfDNA") # needs to be updated with the updated get_sample_ids function
# wbc_ids <- get_sample_ids("new_samples_WBC")

# #################################### this part is just a sanity check
# tumor_sample <- return_sample(tumor_ids, "-cfDNA|_cfDNA|-Baseline|_Baseline") # subset to make sure we have the same samples in both groups
# wbc_sample <- return_sample(wbc_ids, "-WBC|_WBC")

# if (identical(tumor_sample, wbc_sample)) {noquote("Good to go!")} else {
#   noquote("Something is seriously wrong, working on it.")
  
#   samples_list <- intersect(tumor_sample, wbc_sample)
#   noquote("FANTASTIC JOB! The ids match.")
#     }

# # subset to remove any samples that may not have both tumor and wbc samples
# tumor_ids <- subset_list(tumor_ids, samples_list, "-cfDNA|_cfDNA|-Baseline|_Baseline")
# wbc_ids <- subset_list(wbc_ids, samples_list, "-WBC|_WBC")

# # identify the vcf files
# tumor_vcf_paths <- as.list(identify_vcf_files("new_samples", tumor_ids))
# wbc_vcf_paths <- as.list(identify_vcf_files("new_samples", wbc_ids))

# # load vcf files, add the necessary prefix to the paths
# find_common_variants(samples_list, 
# 	tumor_vcf_paths, 
# 	wbc_vcf_paths, 
# 	path_to_common = path_to_common_variants,
# 	save = TRUE)

####################################
# FOR THE KIDNEY SAMPLES - third batch
####################################
tumor_ids <- get_sample_ids("kidney_samples", "ctDNA")
wbc_ids <- get_sample_ids("kidney_samples", "WBC")

#################################### this part is just a sanity check
tumor_sample <- return_sample(tumor_ids, "-cfDNA|_cfDNA|-Baseline|_Baseline") # subset to make sure we have the same samples in both groups
wbc_sample <- return_sample(wbc_ids, "-WBC|_WBC")

if (identical(tumor_sample, wbc_sample)) {noquote("Good to go!")} else {
  noquote("Something is seriously wrong, working on it.")
    }

samples_list <- intersect(tumor_sample, wbc_sample)

# subset to remove any samples that may not have both tumor and wbc samples
tumor_ids <- subset_list(tumor_ids, samples_list, "-cfDNA|_cfDNA|-Baseline|_Baseline")
wbc_ids <- subset_list(wbc_ids, samples_list, "-WBC|_WBC")

# identify the vcf files
tumor_vcf_paths <- as.list(identify_vcf_files("kidney_samples", tumor_ids))
wbc_vcf_paths <- as.list(identify_vcf_files("kidney_samples", wbc_ids))

# load vcf files, add the necessary prefix to the paths
find_common_variants(samples_list, 
	tumor_vcf_paths, 
	wbc_vcf_paths, 
	path_to_common = path_to_common_variants,
	save = TRUE)


proc.time() - ptm # end timer
noquote(format(Sys.time(), "%a %b %d %X %Y"))

