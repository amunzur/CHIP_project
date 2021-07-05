ptm <- proc.time() # start timer
noquote(format(Sys.time(), "%a %b %d %X %Y"))

# make the log file and create a file connection to write to it.
log_file_name <- paste0(length(list.files("/groups/wyattgrp/users/amunzur/log/")), ".log")
file.create(log_file_name)
log_con <- file(log_file_name, open = "a") # open in the append more so that we can add stuff cumulatively

library(dplyr)
library(vcfR)

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
	# Given the smaple type and a list of samples,load the relevant vcf files. 

	# sample type: given as a string, either "cfDNA" or "WBC". This helps locate the dirname thathas the files.
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

load_vcf_from_list <- function(vcf_paths_list) {

	# super simple function that loads the relevant vcf files given a list of paths, and saves them in a list as dataframes
	vcf_list <- list()

	for (vcf_path in vcf_paths_list){

		vcf_file <- read.vcfR(vcf_path, verbose = FALSE)
		vcf_file <- vcf_file@fix # only grab the portion we are interested in

		vcf_list <- append(vcf_list, vcf_file)

	} # end of for loop

	return(vcf_list)

}

find_common_variants <- function(samples_list, tumor_vcf_paths, wbc_vcf_paths, path_to_common, save){

	i <- 1
	while (i <= length(samples_list)){

		sample <- samples_list[i]
		tumor_path <- tumor_vcf_paths[grep(sample, tumor_vcf_paths)]
		wbc_path <- wbc_vcf_paths[grep(sample, wbc_vcf_paths)]

		if (length(tumor_path) > 0 & length(wbc_path) > 0){
			
			print(c("Found both files, started sample", sample))
			# these are actually lists, we will extract the actual path if we can identify both the ctDNA and WBC sample
			tumor_path <- tumor_path[[1]]
			wbc_path <- wbc_path[[1]]

			# load the vcf files 
			print("Started loading files.")
			tumor <- read.vcfR(tumor_path, verbose = FALSE)
			wbc <- read.vcfR(wbc_path, verbose = FALSE)

			# only get what we need from the files
			tumor <- as.data.frame(tumor@fix)[, c("CHROM",  "POS", "REF", "ALT", "FILTER")]
			wbc <- as.data.frame(wbc@fix)[, c("CHROM",  "POS", "REF", "ALT", "FILTER")]

			# remove the variants that failed
			print("Filtering.")
			tumor <- filter(tumor, FILTER == "PASS")
			wbc <- filter(wbc, FILTER == "PASS")

			# now merge to keep the variants if they exist in both files
			print("Merging.")
			combined <- merge(tumor, wbc)
			path_combined <- paste0(path_to_common, sample, ".csv")
			
			if (dim(combined)[1] != 0) {
				print(c("Writing to csv. Found", dim(combined)[1], "common variants."))
				#cat(c("Writing to csv. Found", dim(combined)[1], "common variants."), , file = log_con)
				if (save == TRUE) {write.csv(combined, path_combined)} else {"No common variants."}}

			print("\n")

		} else {
			print(c("Skipped sample:", sample))
			#cat(c("Skipped sample:", sample), file = log_con)

		}
		i = i + 1

	} # end of while loop

} # end of function 

# get the sample ids and sort 
tumor_ids <- as.list(sort(grep("^[GU]", list.files("/groups/wyattgrp/users/amunzur/chip_project/finland_bams/ctDNA_prognosis_ORIGINAL", pattern = ".bam$"), value = TRUE)))
wbc_ids <- as.list(sort(grep("^[GU]", list.files("/groups/wyattgrp/users/amunzur/chip_project/finland_bams/GU_finland_download_ORIGINAL", pattern = ".bam$"), value = TRUE)))

# subset
tumor_sample <- return_sample(tumor_ids, "-cfDNA")
wbc_sample <- return_sample(wbc_ids, "-WBC")

if (identical(tumor_sample, wbc_sample)) {noquote("Good to go!")} else {
  noquote("Something is seriously wrong, working on it.")
  
  samples_list <- intersect(tumor_sample, wbc_sample)
  noquote("All good now.")
    }

tumor_ids <- subset_list(tumor_ids, samples_list, "-cfDNA")
wbc_ids <- subset_list(wbc_ids, samples_list, "-WBC")

# identify the vcf files
tumor_vcf_paths <- as.list(identify_vcf_files("cfDNA", samples_list))
wbc_vcf_paths <- as.list(identify_vcf_files("WBC", samples_list))

# add the prefixes to both to include the absolute path 
tumor_vcf_paths <- lapply(tumor_vcf_paths, function(some_string) paste0("/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/ctDNA_prognosis/", some_string))
wbc_vcf_paths <- lapply(wbc_vcf_paths, function(some_string) paste0("/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/GU_finland_download/", some_string))


# load vcf files, add the necessary prefix to the paths
find_common_variants(samples_list, 
	tumor_vcf_paths, 
	wbc_vcf_paths, 
	path_to_common = "/groups/wyattgrp/users/amunzur/chip_project/common_variants_new/",
	save = TRUE)

proc.time() - ptm # end timer
noquote(format(Sys.time(), "%a %b %d %X %Y"))

