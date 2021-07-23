library(tidyverse)

snv_df_path <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/curated_muts.csv" # mutations after manual curation 
coverage_files_path <- "/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/individual_samples/new_finland_download" # the individual coverage info files
MAPQ_files_path <- "/groups/wyattgrp/users/amunzur/chip_project/metrics/MAPQ_information/individual_samples/new_finland_download"
snv_df <- as.data.frame(read_csv(snv_df_path))

###################################################
# COVERAGE 
###################################################
i <- 1
tumor_coverage_list <- list()
wbc_coverage_list <- list()

while (i <= dim(snv_df)[1]){

	i_row <- snv_df[i, ]

	# load the samples
	tumor <- as.data.frame(read.delim(file.path(coverage_files_path, grep(i_row$SAMPLE_ID, list.files(coverage_files_path), value = TRUE)), header = FALSE)) # path to coverage
	wbc <- as.data.frame(read.delim(file.path(coverage_files_path, grep(i_row$sample_n, list.files(coverage_files_path), value = TRUE)), header = FALSE)) # path to_coverage

	names(tumor) <- c("CHROM", "POSITION", "COVERAGE")
	names(wbc) <- c("CHROM", "POSITION", "COVERAGE")

	tumor_coverage <- filter(tumor, POSITION == i_row$POSITION, CHROM == i_row$CHROM)[, 3]
	wbc_coverage <- filter(wbc, POSITION == i_row$POSITION, CHROM == i_row$CHROM)[, 3]

	# if the position isn't found in the coverage files, that means it isn't in the gene panel. append NA. 
	if (length(tumor_coverage) == 0) {tumor_coverage_list <- append(NA, tumor_coverage_list)} else {tumor_coverage_list <- append(tumor_coverage, tumor_coverage_list)}
	if (length(wbc_coverage) == 0) {wbc_coverage_list <- append(NA, wbc_coverage_list)} else {wbc_coverage_list <- append(wbc_coverage, wbc_coverage_list)}

	i <- i + 1
	print(i)

}

###################################################
# MAPPING QUALITY 
###################################################
i <- 1
tumor_mq_list <- list()
wbc_mq_list <- list()

while (i <= dim(snv_df)[1]){

	i_row <- snv_df[i, ]

	# load the samples, this time we have the header
	tumor <- as.data.frame(read.delim(file.path(MAPQ_files_path, grep(i_row$SAMPLE_ID, list.files(MAPQ_files_path), value = TRUE))))
	wbc <- as.data.frame(read.delim(file.path(MAPQ_files_path, grep(i_row$sample_n, list.files(MAPQ_files_path), value = TRUE))))

	tumor_mq <- filter(tumor, POSITION == i_row$POSITION, CHROM == i_row$CHROM)[, 3]
	wbc_mq <- filter(wbc, POSITION == i_row$POSITION, CHROM == i_row$CHROM)[, 3]

	# if the position isn't found in the coverage files, that means it isn't in the gene panel. append NA. 
	if (length(tumor_mq) == 0) {tumor_mq_list <- append(NA, tumor_mq_list)} else {tumor_mq_list <- append(tumor_mq, tumor_coverage_list)}
	if (length(wbc_mq) == 0) {wbc_mq_list <- append(NA, wbc_mq_list)} else {wbc_mq_list <- append(wbc_coverage, wbc_coverage_list)}

	i <- i + 1
	print(i)

}

###################################################
# ADD COVERAGE AND MAPQ TO THE CURATED_MUTS.CSV 
###################################################
tumor_coverage_list <- list()
wbc_coverage_list <- list()
tumor_mq_list <- list()
wbc_mq_list <- list()

snv_df$coverage_t <- unlist(tumor_coverage_list)
snv_df$coverage_n <- unlist(wbc_coverage_list)
snv_df$MAPQ_t <- unlist(tumor_mq_list)
snv_df$MAPQ_n <- unlist(wbc_mq_list)

# save to csv 
snv_df_path <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/curated_muts_new.csv" # mutations after manual curation 
write_csv(snv_df, snv_df_path)

