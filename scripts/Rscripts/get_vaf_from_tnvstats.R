library(tidyverse)

use_annovar <- FALSE
vaf_threshold <- 0.10
normal_read_support <- 10 # variants with less than 10 rs are EXCLUDED
tumor_read_support <- 8

#INPUTS
dir_to_tnvstats <- "/groups/wyattgrp/users/amunzur/chip_project/tnvstats/new_finland_download"
dir_to_mutect_results <- "/groups/wyattgrp/users/amunzur/chip_project/common_variants/new_finland_download"
path_to_bamlist <- "/groups/wyattgrp/users/amunzur/chip_project/finland_bams/bamslist/tnvstats_bamList_new_finland_download.txt" # bam list used to run tnvstats

# OUTPUTS
path_to_metrics <- "/groups/wyattgrp/users/amunzur/chip_project/metrics/file_situation/new_finland_download"
output_path <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download"

# tiny function to get the file extension
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# helper function to evaluate if we have all the files we need for merging
am_i_good_to_go <- function(dir_to_tnvstats, dir_to_mutect_results, path_to_bamlist, path_to_metrics, extra_safe){

  tnvstats_exists <- list()
  tnvstats_full <- list()

  mutect2_exists <- list()
  mutect2_full <- list()

  samples_df <- as.data.frame(read.delim(path_to_bamlist, header = FALSE))
  names(samples_df) <- c("name_tumor", "path_tumor", "name_wbc", "path_wbc")
  samples <- as.character(samples_df$name_tumor) # list of all tumor samples

  for (sample in samples) {

    print(sample)

    # CHECK TNVSTATS
    path_to_tnvstats <- paste0(file.path(dir_to_tnvstats, sample), "/", sample, ".tnvstat")
    
    # does tnvstats exist?
    if (file.exists(path_to_tnvstats)) {
      tnvstats_exists <- append(tnvstats_exists, TRUE)
      tnvstats <- read.delim(path_to_tnvstats)

      # does tnvstats have data?
      if (dim(tnvstats)[1] == 0) {tnvstats_full <- append(tnvstats_full, FALSE)} else {tnvstats_full <- append(tnvstats_full, TRUE)}

      } else {append(tnvstats_exists, FALSE)}
    
    if (file.exists(paste0(file.path(dir_to_tnvstats, sample), ".tnvstat"))) {tnvstats_exists <- append(tnvstats_exists, TRUE)} else {append(tnvstats_exists, FALSE)}

    # CHECK MUTECT2
    sample_name <- strsplit(sample, "-cfDNA|_cfDNA|-Baseline|_Baseline")[[1]][[1]] # just like samples, but the name here is cleaned up
    path_to_mutect_results <- file.path(dir_to_mutect_results, paste0(sample_name, ".csv")) # path to common variants

    if (file.exists(path_to_mutect_results)) {
      mutect2_exists <- append(mutect2_exists, TRUE)
      mutect2 <- as.data.frame(read.csv(path_to_mutect_results))

      # does mutect have data?
      if (dim(mutect2)[1] == 0) {mutect2_full <- append(mutect2_full, FALSE)} else {mutect2_full <- append(mutect2_full, TRUE)}

  } else {append(mutect2_exists, FALSE)}

  } # end of for loop

  # since function takes a long time to run, we might wanna be extra safe. this would have been much faster if written in bash. oh well.
  if (extra_safe == TRUE){

    write_delim(as.data.frame(tnvstats_exists), paste(path_to_metrics, "tnvstats_exists", sep = "_"))
    write_delim(as.data.frame(tnvstats_full), paste(path_to_metrics, "tnvstats_full", sep = "_"))
    write_delim(as.data.frame(mutect2_exists), paste(path_to_metrics, "mutect2_exists", sep = "_"))
    write_delim(as.data.frame(mutect2_full), paste(path_to_metrics, "mutect2_full", sep = "_"))

  }

  # make a df with all the info we collected 
  metrics_df <- as.data.frame(cbind(samples, tnvstats_exists, tnvstats_full, mutect2_exists, mutect2_full))
  # names(metrics_df) <- c('samples', 'tnvstats_exists', 'tnvstats_full', 'mutect2_exists', 'mutect2_full')

  # and write to file
  write_csv(metrics_df, path_to_metrics)

  # we need to honor the name of the function 
  print("YOU ARE GOOD TO GO!")

  return(metrics_df)

} # end of function


compare_tnvstats_and_mutect <- function(dir_to_tnvstats, dir_to_mutect_results, output_path){

  # Check if the merged file (mutect and tnvstats) has been computed before. If yes, skip sample.
  for (file_name in list.files(dir_to_tnvstats)) {

    sample_name <- strsplit(file_name, "-cfDNA|_cfDNA|-Baseline|_Baseline")[[1]][[1]]
  
    path_to_mutect_results <- file.path(dir_to_mutect_results, paste0(sample_name, ".csv"))
    path_to_tnvstats <- file.path(dir_to_tnvstats, file_name, paste0(file_name, ".tnvstat"))

    if (file.exists(file.path(output_path, paste0(sample_name, ".csv")))) { message(c("Merged already, skipping ", sample_name)) # so that we dont waste time merging files that have been merged already
  
  } else {   # if we havent computed the merged file for the group of samples

    # if there is a match in the variant files for tnvstats read both files
    if ((file.exists(path_to_mutect_results)) & (file.exists(path_to_tnvstats))) {
            
      mutect <- as.data.frame(read.csv(path_to_mutect_results)) %>% select(CHROM, POS, ALT, REF)
      tnvstats <- read.delim(path_to_tnvstats)
      
      if(dim(tnvstats)[1] == 0) {message(c("Tnvstats empty, failed merge ", sample_name))
      } else {
        message(c("MERGING ", sample_name))
        merged <-inner_join(mutect, tnvstats, by=c("CHROM" = "chrom", "POS" = "pos"))
        write_csv(merged, file.path(output_path, paste0(sample_name, ".csv")))}

  }  else {message(c("Missing file ", sample_name))

    }
        
      } 
      
    } # end of for loop 
    
  } # end of function

# filter variants based on vaf, read support and indel status. 
filter_variants <- function(myMergedData, vaf_threshold, normal_read_support, tumor_read_support, consider_indels, save){

  if(consider_indels == FALSE) {variant_df <- subset(myMergedData, nchar(as.character(ALT)) == 1)
  
  i <- 1 # start the counter to go through the rows 
  toremove <- list() # we'll save the idx for removing rows later on
  
  while (i <= dim(variant_df)[1]){ # iterate through rows
    
    alt <- as.vector(variant_df[i, "ALT"]) # tumor allele, the variant we identified
        
    # calculate vafs
    tumor_vaf <- as.vector(variant_df[i, paste0(alt, "AF_t")])
    normal_vaf <- as.vector(variant_df[i, paste0(alt, "AF_n")]) # for WBC vaf, look into the same variant allele, not the basenormal.

    # calculate read support
    tumor_rs <- as.vector(variant_df[i, paste0(alt, "_t")])
    normal_rs <- as.vector(variant_df[i, paste0(alt, "_n")])
    
    if (tumor_vaf < vaf_threshold & tumor_rs >= tumor_read_support) {print(c("TUMOR:", round(tumor_vaf, 2), tumor_rs))} # print to console
    if (normal_vaf < vaf_threshold & normal_rs >= normal_read_support) {print(c("WBC:", round(normal_vaf, 2), normal_rs))}

    if (tumor_vaf > vaf_threshold | normal_vaf > vaf_threshold | tumor_rs < tumor_read_support | normal_rs < normal_read_support) {toremove <- append(toremove, i)} # if either one of them is higher than the threshold, add to the list of removed idx
    
    i <- i + 1 # check the next row

    if (save == TRUE) {
    filtered <- variant_df[-unlist(toremove), ]
    write_csv(filtered, file.path(output_path, paste0("ALL_MERGED_FILTERED", vaf_threshold, ".csv")))}
    
  } # end of while loop

  # consider indels
  } else { 

    variant_df <- subset(myMergedData, nchar(as.character(ALT)) > 1)
    if (save == TRUE) {
    filtered <- variant_df
    write_csv(filtered, file.path(output_path, paste0("ALL_MERGED_FILTERED_INDEL", ".csv")))}

  } # end of if loop

  # after going through all the rows, subset to the muts we will keep

  return(filtered)

} # end of function

# running tnvstats helps exclude the genomic locations we don't want, and save the vaf information for the variants, from the tnvstats.
# metrics_df <- am_i_good_to_go(dir_to_tnvstats, dir_to_mutect_results, path_to_bamlist, path_to_metrics, extra_safe = TRUE)
compare_tnvstats_and_mutect(dir_to_tnvstats, dir_to_mutect_results, output_path)

# now we merge all the individual tnvstats-mutect files we wrote.
setwd(output_path) # makes it easier to locate and upload the files
myMergedData <- do.call(rbind, lapply(list.files(path = output_path, pattern = "^GU-"), read.csv)) # pattern here allows us to merge only certain files
write_csv(myMergedData, file.path(output_path, "ALL_MERGED.csv")) # final file with the variant information from all samples

# TAKING A CLOSER LOOK AT THE MERGED DATA - identifying the variants with vafs within our range
df <- filter_variants(myMergedData, 
  vaf_threshold = 0.10, 
  normal_read_support = 10,
  tumor_read_support = 8, 
  consider_indels = FALSE, 
  save = TRUE)

df <- filter_variants(myMergedData, 
  vaf_threshold = 0.10, 
  normal_read_support = 10,
  tumor_read_support = 8, 
  consider_indels = TRUE, 
  save = TRUE)

