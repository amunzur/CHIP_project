library(tidyverse)
use_annovar <- FALSE

path_to_tnvstats <- "/groups/wyattgrp/users/amunzur/chip_project/tnvstats"
path_to_mutect_results <- "/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered"
output_path <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted2"

# tiny function to get the file extension
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


compare_tnvstats_and_mutect <- function(path_to_tnvstats, path_to_mutect_results, output_path, sample_extension){}
  
  for (file_name in list.files(path_to_tnvstats)) {
    
    sample_name <- strsplit(file_name, "-cfDNA")[[1]][[1]]	
    sample_name <- strsplit(sample_name, "RG_")[[1]][[2]]
    sample_name <- paste0(sample_name, substrRight(sample_name, 3)) # this last line helps add the right extension 
    
    # if there is a match in the variant files for tnvstats read both files
    if (file.exists(file.path(path_to_mutect_results, sample_name))) { 
      
      # read.delim or csv, based on file type
      if (substrRight(sample_name, 3) == "tsv") {mutect <- as.data.frame(read.delim(file.path(path_to_mutect_results, sample_name))) %>% select(CHROM, POS, ALT, REF)}
      if (substrRight(sample_name, 3) == "csv") {mutect <- as.data.frame(read.csv(file.path(path_to_mutect_results, sample_name))) %>% select(CHROM, POS, ALT, REF)}
      
      if (file.exists(file.path(path_to_tnvstats, file_name, paste0(file_name, ".tnvstat")))) {
        
        print(file_name)
        tnvstats <- read.delim(file.path(path_to_tnvstats, file_name, paste0(file_name, ".tnvstat")))
        merged <-inner_join(mutect, tnvstats, by=c("CHROM" = "chrom", "POS" = "pos"))
        write_delim(merged, file.path(output_path, str_replace(sample_name, "tsv", "csv")))
        
      } else {print("tnvstats missing.")}
      
    } # end of if loop 
    
  } # end of for loop
  









if (use_annovar){
  
  path_to_tnvstats <- "/groups/wyattgrp/users/amunzur/chip_project/tnvstats"
  path_to_mutect_results <- "/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered"
  
  for (file_name in list.files(path_to_tnvstats)) {
    
    sample_name <- strsplit(file_name, "-cfDNA")[[1]][[1]]	
    sample_name <- strsplit(sample_name, "RG_")[[1]][[2]]
    sample_name <- paste0(sample_name, ".tsv") # this last line helps add the right extension 
    
    # if there is a match in the variant files for tnvstats read both files
    if (file.exists(file.path(path_to_mutect_results, sample_name))) { 
      mutect <- as.data.frame(read.delim(file.path(path_to_mutect_results, sample_name))) %>% select(CHROM, POS, ALT, REF)
      
      if (file.exists(file.path(path_to_tnvstats, file_name, paste0(file_name, ".tnvstat")))) {
        
        tnvstats <- read.delim(file.path(path_to_tnvstats, file_name, paste0(file_name, ".tnvstat")))
        print(file_name)
        merged <-inner_join(mutect, tnvstats, by=c("CHROM" = "chrom", "POS" = "pos"))
        write_delim(merged, file.path("/groups/wyattgrp/users/amunzur/chip_project/subsetted", str_replace(sample_name, "tsv", "csv")))
        
        
      } else print("tnvstats missing.")
      
    } # end of if loop 
    
  } # end of for loop
  
  # rbind all samples here
  setwd("/groups/wyattgrp/users/amunzur/chip_project/subsetted")
  myMergedData <- 
    do.call(rbind,
            lapply(list.files(path = "/groups/wyattgrp/users/amunzur/chip_project/subsetted"), read.csv))
  
  write_csv(myMergedData, "/groups/wyattgrp/users/amunzur/chip_project/subsetted/ALL_MERGED.csv")
  
  i <- 1 
  while (i <= dim(myMergedData)[1]){
    
    alt <- as.vector(myMergedData[i, "ALT"]) # tumor allele, the variant we identified
    base_normal <- as.vector(myMergedData[i, "base_normal"])
    
    if (nchar(alt) == 1) {
      tumor_vaf <- as.vector(myMergedData[i, paste0(alt, "AF_t")])
      if (tumor_vaf < 0.30) {print(c("TUMOR:", tumor_vaf))}}
    
    normal_vaf <- as.vector(myMergedData[i, paste0(base_normal, "AF_n")])
    if (normal_vaf < 0.50) {print(c("WBC:", normal_vaf))}
    # print(c("WBC:", normal_vaf))
    
    i <- i + 1
    
  } # end of while loop
  
} else {
  
  path_to_tnvstats <- "/groups/wyattgrp/users/amunzur/chip_project/tnvstats"
  path_to_mutect_results <- "/groups/wyattgrp/users/amunzur/chip_project/common_variants"
  
  for (file_name in list.files(path_to_tnvstats)) {
    
    sample_name <- strsplit(file_name, "-cfDNA")[[1]][[1]]	
    sample_name <- strsplit(sample_name, "RG_")[[1]][[2]]
    sample_name <- paste0(sample_name, ".csv") # this last line helps add the right extension 
  
    # if there is a match in the variant files for tnvstats read both files
    if (file.exists(file.path(path_to_mutect_results, sample_name))) { 
      
      mutect <- as.data.frame(read.csv(file.path(path_to_mutect_results, sample_name))) %>% dplyr::select(CHROM, POS, ALT, REF)
      
      if (file.exists(file.path(path_to_tnvstats, file_name, paste0(file_name, ".tnvstat")))) {
        
        tnvstats <- read.delim(file.path(path_to_tnvstats, file_name, paste0(file_name, ".tnvstat")))
        print(file_name)
        merged <-inner_join(mutect, tnvstats, by=c("CHROM" = "chrom", "POS" = "pos"))
        write_delim(merged, file.path("/groups/wyattgrp/users/amunzur/chip_project/subsetted2", str_replace(sample_name, "tsv", "csv")))
        
        
      } else print("tnvstats missing.")
      
    } # end of if loop 
    
  } # end of for loop
  
  # rbind all samples here
  setwd("/groups/wyattgrp/users/amunzur/chip_project/subsetted")
  myMergedData <- 
    do.call(rbind,
            lapply(list.files(path = "/groups/wyattgrp/users/amunzur/chip_project/subsetted"), read.csv))
  
  write_csv(myMergedData, "/groups/wyattgrp/users/amunzur/chip_project/subsetted/ALL_MERGED.csv")
  
  i <- 1 
  while (i <= dim(myMergedData)[1]){
    
    alt <- as.vector(myMergedData[i, "ALT"]) # tumor allele, the variant we identified
    base_normal <- as.vector(myMergedData[i, "base_normal"])
    
    if (nchar(alt) == 1) {
      tumor_vaf <- as.vector(myMergedData[i, paste0(alt, "AF_t")])
      if (tumor_vaf < 0.30) {print(c("TUMOR:", tumor_vaf))}}
    
    normal_vaf <- as.vector(myMergedData[i, paste0(base_normal, "AF_n")])
    if (normal_vaf < 0.50) {print(c("WBC:", normal_vaf))}
    # print(c("WBC:", normal_vaf))
    
    i <- i + 1
    
  } # end of while loop
  
  
  
  
  
  
  
}







