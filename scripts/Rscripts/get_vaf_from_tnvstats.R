library(tidyverse)

use_annovar <- FALSE
vaf_threshold <- 0.20
normal_read_support <- 10 # variants with less than 10 rs are EXCLUDED
tumor_read_support <- 8

#INPUTS
cohort <- "kidney_samples" # only update this

dir_to_tnvstats <- file.path("/groups/wyattgrp/users/amunzur/chip_project/tnvstats", cohort) # dir where the tnvstats are saved to
dir_to_mutect_results <- file.path("/groups/wyattgrp/users/amunzur/chip_project/common_variants", cohort) # dir where common vars between tumor and wbc saved to
path_to_bamlist <- file.path("/groups/wyattgrp/users/amunzur/chip_project/tnvstats/kidney_samples/bamList.txt") # bam list used to run tnvstats

# OUTPUTS
path_to_metrics <- file.path("/groups/wyattgrp/users/amunzur/chip_project/metrics/file_situation", cohort)
output_path <- file.path("/groups/wyattgrp/users/amunzur/chip_project/subsetted", cohort)

# tiny function to get the file extension
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# helper function to evaluate if we have all the files we need for merging
am_i_good_to_go <- function(dir_to_tnvstats, dir_to_mutect_results, path_to_bamlist, path_to_metrics, extra_safe){

  tnvstats_exists <- list() # does tnvstats exist? 
  tnvstats_full <- list() # does tnvstats have data, i.e. at least one column?

  mutect2_exists <- list()
  mutect2_full <- list()

  samples_df <- as.data.frame(read.delim(path_to_bamlist, header = FALSE)) # read in the bam list made to run tunvstats
  names(samples_df) <- c("name_tumor", "path_tumor", "name_wbc", "path_wbc") # name the cols
  samples <- as.character(samples_df$name_tumor) # list of all tumor samples

  for (sample in samples) { # go through each tumor sample

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

  # since function takes a long time to run, we might wanna be extra safe. this would have been much faster if written in bash. oh well. whatevs.
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
  
  } else { # if we havent computed the merged file for the group of samples

    # if there is a match in the variant files for tnvstats read both files
    if ((file.exists(path_to_mutect_results)) & (file.exists(path_to_tnvstats))) {
            
      message(c("Found both files, started merge.", sample_name))
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
modify_variant_cols <- function(df, which_info) {
  # given a df, this function identifies which col contains the vaf and rs for a given alt allele, and creates a new col with the vafs. 
  # for vaf, which_info="vaf", for rs, which_info="rs"
  
  if (which_info == "vaf") {

    tumor_vaf_names <- paste0(df$ALT, "AF_t") # only colnames that contain the vaf for the alt allele
    wbc_vaf_names <- paste0(df$ALT, "AF_n")

  } else {

    tumor_vaf_names <- paste0(df$ALT, "_t") # only colnames that contain the vaf for the alt allele
    wbc_vaf_names <- paste0(df$ALT, "_n")

  }

  # get the col index and subset, then add as new cols
  idx <- as.numeric(sapply(tumor_vaf_names, function(x) grep(x, colnames(df)))) 
  df$tumor_info <- as.numeric(df[cbind(seq_along(tumor_vaf_names), idx)])

  idx <- as.numeric(sapply(wbc_vaf_names, function(x) grep(x, colnames(df)))) 
  df$wbc_info <- as.numeric(df[cbind(seq_along(wbc_vaf_names), idx)])

  # rename, based on whether we added rs or vaf
  names(df)[names(df) == 'tumor_info'] <- paste0("tumor_", which_info)
  names(df)[names(df) == 'wbc_info'] <- paste0("wbc_", which_info)

  return(df)

} 

filter_variants <- function(myMergedData, vaf_threshold, normal_read_support, tumor_read_support, consider_indels, save){

  if(consider_indels == FALSE) {variant_df <- subset(myMergedData, nchar(as.character(ALT)) == 1)

    if (save == TRUE) {write_csv(variant_df, file.path(output_path, paste0("ALL_SNVs", ".csv")))}

    variant_df <- modify_variant_cols(variant_df, "vaf")
    variant_df <- modify_variant_cols(variant_df, "rs")

    # filter based on vaf and rs
    filtered <- variant_df %>%
      filter(tumor_vaf <= vaf_threshold, 
             wbc_vaf <= vaf_threshold, 
             tumor_rs >= tumor_read_support, 
             wbc_rs >= normal_read_support) 

    if (save == TRUE) {write_csv(filtered, file.path(output_path, paste0("ALL_MERGED_FILTERED_", vaf_threshold, ".csv")))}
  
  # consider indels
  } else { 

    variant_df <- subset(myMergedData, nchar(as.character(ALT)) > 1)
    if (save == TRUE) {
    filtered <- variant_df
    write_csv(filtered, file.path(output_path, paste0("ALL_MERGED_FILTERED_INDEL", ".csv")))}

  } # end of if to consider indels loop

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
  vaf_threshold = 0.20, 
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

