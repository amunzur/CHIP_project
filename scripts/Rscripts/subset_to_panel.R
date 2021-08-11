library(tidyverse)

path_to_bed <- "/groups/wyattgrp/users/amunzur/chip_project/references/baits_hg38.bed" # bed file with locations of interest

path_to_df <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/ALL_MERGED.csv" # df to be subsetted
path_to_df_subsetted <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/ALL_MERGED_panel.csv"

path_to_df <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/ALL_MERGED_FILTERED_0.2.csv"
path_to_df_subsetted <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/ALL_MERGED_FILTERED_0.2_panel.csv"

path_to_df <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/ALL_MERGED_FILTERED_INDEL.csv"
path_to_df_subsetted <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/ALL_MERGED_FILTERED_INDEL_panel.csv"

path_to_df <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/curated_muts.csv"
path_to_df_subsetted <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/curated_muts_panel.csv"

# kidney samples
path_to_bed <- "/groups/wyattgrp/users/amunzur/chip_project/misc_files/HumanOncologyPanel_Design_Files/bed_file_capture.bed" # bed file with locations of interest
path_to_df <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/kidney_samples_second/curated_muts.csv"
path_to_df_subsetted <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/kidney_samples_second/curated_muts_panel.csv"


subset_to_panel <- function(path_to_bed, path_to_df, path_to_df_subsetted) {
	bed <- as.data.frame(read.delim(path_to_bed, header = FALSE))
	if (ncol(bed) == 4) {names(bed) <- c("chrom", "start", "stop", "gene")} else {names(bed) <- c("chrom", "start", "stop")}
	df <- as.data.frame(read_csv(path_to_df))

	to_keep <- list() # list of positions to remove

	i <- 1 
	print(i)
	while (i <= dim(df)[1]){
		chrom_subsetted <- df[i, 1] # pick the chrom we are at 
		location <- df[i, 2] # pick the location we are at 
		bed_subsetted <- bed %>% filter(chrom == chrom_subsetted) # subset the bed by chrom
		message(c("position:", i))

		j <- 1
		while(j <= dim(bed_subsetted)[1]) {
			start_pos <- bed_subsetted[j, 2]
			end_pos <- bed_subsetted[j, 3]

			if (all(location >= start_pos, location <= end_pos)) {
				to_keep <- append(to_keep, i) # saving the row index of muts we are keeping
				break 
				print("Position in panel.")
				} else {j <- j + 1} # if the location isn't in the panel, go check out the next position in the bed file.
				# print(c(j, "j"))

		} # end of inner while loop

	i <- i + 1 # next identified variant

	} # end of outer while loop - looping through identified variants

	df <- df[unlist(to_keep), ]
	write_csv(df, path_to_df_subsetted)

	return(df)

} # end of function

df <- subset_to_panel(path_to_bed, path_to_df, path_to_df_subsetted)


