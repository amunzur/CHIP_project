# This script uses the coverage metrics to find the coverage at a given SNV. 

snv_df_path <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/curated_muts.csv" # mutations after manual curation 
coverage_files_path <- "/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/individual_samples/new_finland_download" # the individual coverage info files
snv_df <- as.data.frame(read_csv(snv_df_path))

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

	if (length(tumor_coverage) == 0) {tumor_coverage_list <- append(0, tumor_coverage_list)} else {tumor_coverage_list <- append(tumor_coverage, tumor_coverage_list)}
	if (length(wbc_coverage) == 0) {wbc_coverage_list <- append(0, wbc_coverage_list)} else {wbc_coverage_list <- append(wbc_coverage, wbc_coverage_list)}

	i <- i + 1
	print(i)

}