library(tidyverse)

get_sample_name <- function(path_to_dir, key) {

	path_list <- readLines(path_to_dir)
	path_list <- lapply(path_list, basename) # get the file name from the files list 
	samples_list <- sort(unlist(lapply(path_list, function(some_name) strsplit(some_name, key)[[1]][[1]])))

	# for cfDNA some samples require a bit more processing
	if (key == "-cfDNA") { 
		samples_list <- sort(unlist(lapply(samples_list, function(some_name) strsplit(some_name, "-1st")[[1]][[1]])))
		samples_list <- sort(unlist(lapply(samples_list, function(some_name) strsplit(some_name, "-16Aug2017")[[1]][[1]]))) }

	return(samples_list)
}

find_missing_samples <- function(ctdna_samples, wbc_samples, path_to_ctdna_original, path_to_wbc_original) {

	# index of missing wbc samples
	idx <- wbc_samples %in% ctdna_samples
	idx <- !idx
	missing_wbc <- wbc_samples[idx]

	idx <- ctdna_samples %in% wbc_samples  
	idx <- !idx
	missing_ctdna <- ctdna_samples[idx]

	print(c("Missing ctDNA samples are:", missing_ctdna))
	print(c("Missing WBC samples are:", missing_wbc))

	# some string cleaning 
	missing_ctdna <- lapply(missing_ctdna, function(some_name) strsplit(some_name, "_")[[1]][[2]])
	missing_wbc <- lapply(missing_wbc, function(some_name) strsplit(some_name, "_")[[1]][[2]])

	# reattach rest of the strings - CTDNA
	path_list_ctdna <- readLines(path_to_ctdna_original)
	wbc_matches <- unique(grep(paste(missing_ctdna, collapse="|"), path_list_ctdna, value=TRUE))

	# reattach rest of the strings - WBC
	path_list_wbc <- readLines(path_to_wbc_original)
	ctdna_matches <- unique(grep(paste(missing_wbc, collapse="|"), path_list_wbc, value=TRUE))

}

remove_samples <- function(to_remove, path_to_file, overwrite){

	# overwrite is boolean
	# only provide the sample id like 18-246 and 18-329
	path_list <- as.list(readLines(path_to_file))
	
	for (x in to_remove) {
		idx <- grep(x, path_list)
		path_list[idx] <- NULL # remove the unwanted sample
	}

	# function will overwrite the file with removed samples if the user specifies
	if (overwrite) {

		file.remove(path_to_file) # delete the old file 

		fileConn <- file(path_to_file) # establish a new file connection, essentially the same path
		writeLines(unlist(path_list), fileConn)
		close(fileConn)

	} # end of if loop

} # end of function

replace_sample_with_full_names <- function(path_to_full_names, samples_list) {

	full_names <- unlist(lapply(as.list(readLines(path_to_full_names)), basename))

	for (sample in samples_list) {

		idx <- which(samples_list == sample) # this will help replace the sample with the full_sample name later on
		sample <- strsplit(sample, "RG_")[[1]][[2]] # remove the RG_string at teh beginning of the sample name
		
		sample_full_name <- grep(sample, full_names, value = TRUE) # grep the full name of the sample
		samples_list[idx] <- sample_full_name # make inplace replacement ot avoid making another variable 

	} #end of for loop 

	return(samples_list)

} # end of function

path_to_ctdna <- "/groups/wyattgrp/users/amunzur/chip_project/finland_bams/bamslist/ctDNA_bams"
path_to_wbc <- "/groups/wyattgrp/users/amunzur/chip_project/finland_bams/bamslist/wbc_bams"

path_to_ctdna_original <- "/groups/wyattgrp/users/amunzur/chip_project/finland_bams/bamslist/ctDNA_original"
path_to_wbc_original <- "/groups/wyattgrp/users/amunzur/chip_project/finland_bams/bamslist/wbc_original"

# remove the samples we arent interested in
remove_samples(path_to_ctdna, c("18-246", "18-329"), FALSE)
remove_samples(path_to_wbc, c("18-246", "18-329"), FALSE)

# this just helps to cleanup the strings and prepare the sample names for matching
ctdna_samples <- get_sample_name(path_to_ctdna, "-cfDNA")
wbc_samples <- get_sample_name(path_to_wbc, "-WBC")

# this prints out the missing samples we need to work on before moving onto the next step
find_missing_samples(ctdna_samples, wbc_samples, path_to_ctdna_original, path_to_wbc_original)

###################################################
# PREPARE BAMSLIST FOR TNVSTATS
###################################################
# bamslist is a delim file with sample names and paths to those samples
# this step assumes that we have a ctdna and wbc match for every sample
ctdna_samples <- sort(ctdna_samples)
wbc_samples <- sort(wbc_samples)

# check if these two lists are identical - if all good to go, continue.
if (identical(ctdna_samples, wbc_samples)) {print("The samples are identical. Good to go!")} else {"Something is seriously wrong."}

# now find the paths associated with these samples and save them in another vector 
ctdna_paths_list <- list()
wbc_paths_list <- list()

for (x in ctdna_samples) {
	
	ctdna_file <- readLines(path_to_ctdna) # load the file with the ctdna paths
	file_path <- grep(x, ctdna_file, value = TRUE, ignore.case = TRUE)
	ctdna_paths_list <- append(ctdna_paths_list, file_path)}

for (x in wbc_samples) {
	
	wbc_file <- readLines(path_to_wbc) # load the file with the wbc paths 
	file_path <- grep(x, wbc_file, value = TRUE, ignore.case = TRUE)
	wbc_paths_list <- append(wbc_paths_list, file_path)}

# so far the samples we have been working with don't have the full sample name, only the patient id. 
# this piece of code below completes the id with the rest of it, not just the patient id. 
ctdna_samples <- replace_sample_with_full_names(path_to_ctdna_original, ctdna_samples)
wbc_samples <- replace_sample_with_full_names(path_to_wbc_original, wbc_samples)

# now save everything to a delim file 
DF <- data.frame(
	ctdna_samples = ctdna_samples, 
	ctdna_paths = unlist(ctdna_paths_list), 
	wbc_samples = wbc_samples, 
	wbc_paths = unlist(wbc_paths_list))

write.csv(DF, "/groups/wyattgrp/users/amunzur/chip_project/finland_bams/bamslist/tnvstats_bamList.csv", row.names = FALSE)

