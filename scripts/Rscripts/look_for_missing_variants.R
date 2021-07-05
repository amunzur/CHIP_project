library(tidyverse)

# Goal of this dscript is to identify the variants not picked by in house variant caller, and figure out why.
# We will first start with 4 interesting variants, we go back to tnvstats for these variants and extract a bunch of useful stats.
path_to_tnvstats <- "/groups/wyattgrp/users/amunzur/chip_project/tnvstats_new"
path_to_curated_combined <-"/groups/wyattgrp/users/amunzur/chip_project/tnvstats_mutect_compared/combined.csv" # all muts are here, inbdicating whether they are found in only one pipeline or both
path_to_final <- "/groups/wyattgrp/users/amunzur/chip_project/tnvstats_mutect_compared/only_mutect_samples_info.csv" # this is the pah to the df that contains tnvstat information for samples that only mutect identified
combined <- as.data.frame(read_csv(path_to_curated_combined))

mutect <- combined %>% filter(STATUS == "mutect_only") # get mutect only

tumor_ids <- lapply(as.list(mutect$SAMPLE_ID), function(some_name) paste(some_name, "bam", sep = ".")) # add the ".bam" to tumor ids 
tumor_ids <- lapply(tumor_ids, function(some_name) file.path(path_to_tnvstats, some_name, paste0(some_name, ".tnvstat")))

tnvstats <- lapply(tumor_ids, read.delim) # load all tnvstats at once, and filter based on the locations we are interested in

i <- 1
filtered_list <- list() # an empty list to append the bits and pieces from the filtered tnvstats
while (i <= dim(mutect)[1]){
	
	tnvstat <-  as.data.frame(tnvstats[i]) # load the tnv stat for the related sample
	filtered <- as.data.frame(unlist(tnvstat %>% filter(chrom == mutect$CHROM[i], pos == mutect$POSITION[i])))
	filtered_list <- append(filtered_list, filtered)
	
	i <- i + 1
	
	print(i)
	print(to_filter)
	print(filtered)
}

df <- as.data.frame(do.call(rbind, filtered_list)) # combine all pieces into one large df
rownames(df) <- NULL # no need for rownames
names(df) <- names(tnvstat)


















