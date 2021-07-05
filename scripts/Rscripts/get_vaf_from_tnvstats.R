library(tidyverse)

use_annovar <- FALSE
vaf_threshold <- 0.10
normal_read_support <- 10
tumor_read_support <- 8

dir_to_tnvstats <- "/groups/wyattgrp/users/amunzur/chip_project/tnvstats_new"
dir_to_mutect_results <- "/groups/wyattgrp/users/amunzur/chip_project/common_variants"
output_path <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted_new"

# tiny function to get the file extension
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

compare_tnvstats_and_mutect <- function(dir_to_tnvstats, dir_to_mutect_results, output_path){
  
  # Check if the merged file (mutect and tnvstats) has been computed before. If yes, skip sample.
  for (file_name in list.files(dir_to_tnvstats)) {

    sample_name <- strsplit(file_name, "-cfDNA")[[1]][[1]]
  
    path_to_mutect_results <- file.path(dir_to_mutect_results, paste0(sample_name, ".csv"))
    path_to_tnvstats <- file.path(dir_to_tnvstats, file_name, paste0(file_name, ".tnvstat"))

    if (file.exists(file.path(output_path, paste0(sample_name, ".csv")))) { print(c("Merge has been done. Skipping", sample_name))
  
  } else {   # if we havent computed the merged file for the group of samples

    # if there is a match in the variant files for tnvstats read both files
    if ((file.exists(path_to_mutect_results)) & (file.exists(path_to_tnvstats))) {
      
      print(file_name)
      
      mutect <- as.data.frame(read.csv(path_to_mutect_results)) %>% select(CHROM, POS, ALT, REF)
      tnvstats <- read.delim(path_to_tnvstats)
      merged <-inner_join(mutect, tnvstats, by=c("CHROM" = "chrom", "POS" = "pos"))
      
      write_csv(merged, file.path(output_path, paste0(sample_name, ".csv"))) # write as csv

  }  else {print(c("Tnvstats missing.", sample_name))}
        
      } 
      
    } # end of for loop 
    
  } # end of function

filter_variants <- function(myMergedData, vaf_threshold, normal_read_support, tumor_read_support, consider_indels, save){

  # consider_indels need to be boolean 
  if (consider_indels == FALSE){
    snv_df <- subset(myMergedData, nchar(as.character(ALT)) == 1)

    i <- 1 # start the counter to go through the rows 
    toremove <- list() # we'll save the idx for removing rows later on
    
    while (i <= dim(snv_df)[1]){ # iterate through rows
      
      alt <- as.vector(snv_df[i, "ALT"]) # tumor allele, the variant we identified
          
      # calculate vafs
      tumor_vaf <- as.vector(snv_df[i, paste0(alt, "AF_t")])
      normal_vaf <- as.vector(snv_df[i, paste0(alt, "AF_n")]) # for WBC vaf, look into the same variant allele, not the basenormal.

      # calculate read support
      tumor_rs <- as.vector(snv_df[i, paste0(alt, "_t")])
      normal_rs <- as.vector(snv_df[i, paste0(alt, "_n")])
      
      if (tumor_vaf < vaf_threshold & tumor_rs >= tumor_read_support) {print(c("TUMOR:", round(tumor_vaf, 2), tumor_rs))} # print to console
      if (normal_vaf < vaf_threshold & normal_rs >= normal_read_support) {print(c("WBC:", round(normal_vaf, 2), normal_rs))}

      if (tumor_vaf > vaf_threshold | normal_vaf > vaf_threshold | tumor_rs < tumor_read_support | normal_rs < normal_read_support) {toremove <- append(toremove, i)} # if either one of them is higher than the threshold, add to the list of removed idx
      
      i <- i + 1 # check the next row
      
    } # end of while loop

  } # end of main if loop

  # after going through all the rows, subset to the muts we will keep
  if (save == TRUE) {

    filtered <- snv_df[-unlist(toremove), ]
    write_csv(filtered, file.path(output_path, paste0("ALL_MERGED_FILTERED_", vaf_threshold, ".csv")))

  }

  return(filtered)

} # end of function

# running tnvstats helps exclude the genomic locations we don't want, and save the vaf information for the variants, from the tnvstats.
compare_tnvstats_and_mutect(dir_to_tnvstats, dir_to_mutect_results, output_path)

# now we merge all the individual tnvstats-mutect files we wrote.
setwd(output_path) # makes it easier to locate and upload the files
myMergedData <- do.call(rbind, lapply(list.files(path = output_path, pattern = "^GU-"), read.csv))
write_csv(myMergedData, file.path(output_path, "ALL_MERGED.csv")) # final file with the variant information from all samples

# TAKING A CLOSER LOOK AT THE MERGED DATA - identifying the variants with vafs within our range
df <- filter_variants(myMergedData, 
  vaf_threshold = 0.20, 
  normal_read_support = 10,
  tumor_read_support = 8, 
  consider_indels = FALSE, 
  save = TRUE)
