# The point of this script is to get sample names from png snapshots, manually curated, and filter the main mutext2 results file based on the names we curate 

library(tidyverse)

# FILES TO UPLOAD 
path_to_snapshots <- "/wyattgrp/Asli/chip_project/snapshot/mutect_results/CHIP_MUTS" # dir that contains igv snapshots, after manual curation
path_to_mutect <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted_new/ALL_MERGED_FILTERED_0.2.csv" # file (not dir) that contains variants, after merging and filtering based on read support and vaf 
path_to_curated_elie <- "/groups/wyattgrp/users/echen/CHIP_project/results/curated2_bg_error_10.tsv"
path_to_mutect_and_elie <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted_new/MUTECT_ELIE_COMBINED.csv" # file that contains elies calls and mutect calls, output of this script

# FILES TO WRITE TO CSV 
path_to_curated_mutect <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted_new/curated_muts.csv" # file that will contain the curated muts, output of this script
path_to_both <- "/groups/wyattgrp/users/amunzur/chip_project/tnvstats_mutect_compared/both.csv"
path_to_mutect_only <- "/groups/wyattgrp/users/amunzur/chip_project/tnvstats_mutect_compared/mutect_only.csv"
path_to_mutect_only_filtered <- "/groups/wyattgrp/users/amunzur/chip_project/tnvstats_mutect_compared/mutect_only_filtered.csv"
path_to_tnvstats_only <- "/groups/wyattgrp/users/amunzur/chip_project/tnvstats_mutect_compared/elie_only.csv"
path_to_curated_combined <-"/groups/wyattgrp/users/amunzur/chip_project/tnvstats_mutect_compared/combined.csv" # all muts are here, inbdicating whether they are found in only one pipeline or both

cleanup_mutect_results <- function(mutect_combined){

	i <- 1
	
	tumor_vafs_list <- list()
	wbc_vafs_list <- list()

	tumor_rs_list <- list()
	wbc_rs_list <- list()

	while (i <= dim(mutect_combined)[1]){

		alt <- as.vector(mutect_combined[i, "ALT"]) # the variant
	    
	    # identify
	    tumor_vaf <- as.vector(mutect_combined[i, paste0(alt, "AF_t")])
	    normal_vaf <- as.vector(mutect_combined[i, paste0(alt, "AF_n")])
	    tumor_rs <- as.vector(mutect_combined[i, paste0(alt, "_t")])
        normal_rs <- as.vector(mutect_combined[i, paste0(alt, "_n")])

        # add to lists 
	    tumor_vafs_list <- append(tumor_vafs_list, tumor_vaf)
		wbc_vafs_list <- append(wbc_vafs_list, normal_vaf)

		tumor_rs_list <- append(tumor_rs_list, tumor_rs)
		wbc_rs_list <- append(wbc_rs_list, normal_rs)

		i <- i + 1

	} # end of while loop
	
	# df with the vafs
	df1 <- data.frame("VAF" = unlist(tumor_vafs_list), 
		"WBC_VAF" = unlist(wbc_vafs_list), 
		"Mutant_Reads" = unlist(tumor_rs_list), 
		"Mutant_Reads_WBC" = unlist(wbc_rs_list))

	# df with misc information, tnvstats has these too so i tried to make them look similar
	df2 <- mutect_combined %>% select(CHROM, POSITION, SAMPLE_ID, sample_n, ALT, REF)

	# and now combine them all 
	df3 <- cbind(df1, df2)

	return(df3)

} # end of function

tnvstats <- read.delim(path_to_curated_elie) # read elie's results
tnvstats$WBC_ID <- lapply(as.list(as.vector(tnvstats$WBC_ID)), function(some_name) gsub("^.*?GU","GU", some_name)) # remove everything before GU so that the ids match across dfs

mutect <- read.csv(path_to_mutect) # read mutect results

# get the bam id, chr and location from the snapshot names
snapshot_names <- unlist(lapply(as.list(list.files(path_to_snapshots)), function(some_name) strsplit(some_name, "_")[[1]]))
SAMPLE_ID <- snapshot_names[seq(1, length(snapshot_names), 3)] # starting from 1, every 3rd element is the sample id
CHROM <- snapshot_names[seq(2, length(snapshot_names), 3)] # starting from 2, every 3rd elements is the chrom
POSITION <- unlist(lapply(as.list(snapshot_names[seq(3, length(snapshot_names), 3)]), function(some_name) strsplit(some_name, ".", fixed = TRUE)[[1]][[1]])) # starting from 3, every 3rd elements is the position

df <- data.frame(SAMPLE_ID, CHROM, POSITION) # this df has the chrom and pos info about the manually curated variants
rownames(df) <- NULL # no need to rownames

####################################################
# FILTER THE ALL_MERGED_FILTERED_0.2.csv FILE BASED ON MANUAL CURATION
####################################################

idx <- na.omit(match(names(mutect), c("sample_t", "CHROM", "POS"))) # making sure the names match
mutect$POS <- as.character(mutect$POS)
df$POSITION <- as.character(df$POSITION)
tnvstats$POSITION <- as.character(tnvstats$POSITION)

mutect_combined <- inner_join(df, mutect, by = c("CHROM" = "CHROM", "POSITION" = "POS", "SAMPLE_ID" = "sample_t"))
mutect_combined <- cleanup_mutect_results(mutect_combined)
mutect_combined <- mutect_combined %>% select(CHROM, POSITION, SAMPLE_ID, VAF, WBC_VAF, Mutant_Reads, ALT, REF, sample_n) #reorder

# some cleaning up for the same names
mutect_combined$SAMPLE_ID <- unlist(lapply(as.list(as.vector(mutect_combined$SAMPLE_ID)), function(some_name) strsplit(some_name, ".", fixed = TRUE)[[1]][[1]]))
mutect_combined$sample_n <- unlist(lapply(as.list(as.vector(mutect_combined$sample_n)), function(some_name) strsplit(some_name, ".", fixed = TRUE)[[1]][[1]]))

write_csv(mutect_combined, path_to_curated_mutect)

####################################################
# COMPARE ELIE TO CURATED MUTECT
####################################################
both <- inner_join(mutect_combined, tnvstats, by = c("CHROM", "POSITION", "SAMPLE_ID"))
both$STATUS <- "both"
both <- as.data.frame(apply(both, 2, unlist)) # we have to do this otherwise it complains while saving

mutect_only <- dplyr::setdiff(mutect_combined[, c("CHROM", "POSITION", "SAMPLE_ID")], tnvstats[, c("CHROM", "POSITION", "SAMPLE_ID")])
mutect_only$STATUS <- "mutect_only"
mutect_only <- inner_join(mutect_only, mutect_combined, by = c("CHROM", "POSITION", "SAMPLE_ID")) # this helps add the missing cols after set diff
mutect_only <- as.data.frame(apply(mutect_only, 2, unlist)) # we have to do this otherwise it complains while saving

tnvstats_only <- dplyr::setdiff(tnvstats[, c("CHROM", "POSITION", "SAMPLE_ID")], mutect_combined[, c("CHROM", "POSITION", "SAMPLE_ID")])
tnvstats_only$STATUS <- "tnvstats_only"
tnvstats_only <- inner_join(tnvstats_only, tnvstats, by = c("CHROM", "POSITION", "SAMPLE_ID"))
tnvstats_only <- as.data.frame(apply(tnvstats_only, 2, unlist)) # we have to do this otherwise it complains while saving

# drop some from mutect_only, these are likely to be wrong
idx <- which(duplicated(mutect_only[, 1:2]))
mutect_only_filtered <- mutect_only[-idx, ]

# write these files to csv
write_csv(both, path_to_both)
write_csv(mutect_only, path_to_mutect_only)
write_csv(mutect_only_filtered, path_to_mutect_only_filtered)
write_csv(tnvstats_only, path_to_tnvstats_only)

# MAKE ONE BIG DF FOR PLOTTING 
both <- both %>% 
	select(CHROM, POSITION, SAMPLE_ID, WBC_ID, VAF.x, WBC_VAF.x, REF.x, ALT.x, STATUS)

mutect_only_filtered <- mutect_only_filtered %>% 
	select(CHROM, POSITION, SAMPLE_ID, sample_n, VAF, WBC_VAF, REF, ALT, STATUS)

tnvstats_only <- tnvstats_only %>% 
	select(CHROM, POSITION, SAMPLE_ID, WBC_ID, VAF, WBC_VAF, REF, ALT, STATUS)

# make sure they have the same name
names(both) <- c('CHROM', 'POSITION', 'SAMPLE_ID', 'WBC_ID', 'VAF', 'WBC_VAF', 'REF', 'ALT', 'STATUS')
names(mutect_only_filtered) <- names(both)
names(tnvstats_only) <- names(both)

combined <- rbind(both, mutect_only_filtered, tnvstats_only)
write_csv(combined, path_to_curated_combined)