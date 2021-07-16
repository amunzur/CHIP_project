# this script is for transfering WBC cell samples from finland servers to ours. 
# we only copy samples in the PATH_samples_df.

# load the file with sample list 
PATH_samples_df <- "/groups/wyattgrp/users/amunzur/chip_project/data/CHIP_list_of_samples.txt"
samples_df <- read.delim(PATH_samples_df, sep = "\t", header = TRUE)
WBC_samples <- as.list(as.vector(samples_df$Sample_sequencing_name_server))
WBC_samples_index <- lapply(WBC_samples, function(some_name) paste(some_name, "bai", sep = "."))

source <- "rsync -aP vpc@tambio.uta.fi:/home/annalam/datasets/gu_cohort/72_gene/alignments"
dest <- "/groups/wyattgrp/data/bam/ctDNA_prognosis"

source_bams <- lapply(WBC_samples, function(some_name) paste(source, some_name, sep = "/"))
source_bai <- lapply(WBC_samples_index, function(some_name) paste(source, some_name, sep = "/"))

# concatenate 
FINAL_bams <- paste(source_bams, rep(dest, length(source_bams)))
FINAL_bai <- paste(source_bai, rep(dest, length(source_bams)))

FINAL_list <- c(FINAL_bams, FINAL_bai)

# write to file
writeLines(unlist(lapply(FINAL_list, paste, collapse=" ")), "/groups/wyattgrp/users/amunzur/chip_project/scripts/transfer_ctDNA.txt")
