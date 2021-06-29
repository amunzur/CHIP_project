library(vcfR)
library(tidyverse)

path_to_filtered_mutect <- "/groups/wyattgrp/users/amunzur/chip_project/mutect_results_filtered/" # this helps get read depth
path_to_common <- "/groups/wyattgrp/users/amunzur/chip_project/bed_variants/" # this helps get the reads that pass germline filtering and found in both smaple groups
path_to_pileup <- "/groups/wyattgrp/users/amunzur/chip_project/variant_bed_files/"

type <- "ctDNA_prognosis" # ctDNA_prognosis or GU_finland_download
key <- "-cfDNA" # -cfDNA or -WBC

setwd(paste0(path_to_filtered_mutect))
for (file in list.files(key, pattern = '.bam_vcf_FILTERED_vcf$')){

    vcf <- read.vcfR(file, verbose = FALSE)

    if (dim(vcf@fix)[1] == 0) {
      print("No mutations called.") # skip if the file is empty
    } else {
      
      common_name <- paste(strsplit(file, key)[[1]][[1]], "tsv", sep = ".") # this has the filtered annovar results that we will use to subset the mutect results
      common <- as.data.frame(read.delim(paste0(path_to_common, common_name)), header = TRUE) %>% select(CHROM, POS)

      # GET VARIANT ALLELE READ DEPTH FROM VCF
      gt <- extract_gt_tidy(vcf) %>% select(gt_DP) # this has all the format information in separate cols
      variants <- vcfR2tidy(vcf)$fix %>% select(CHROM, POS) # chrom and pos info is taken from here
      vcf <- cbind(variants, gt) # redefine a vcf object with what we need only

      # subset the vcf based on annovar results 
      joined <- inner_join(vcf, common, by = c("CHROM", "POS")) # this df has the read depths of the variants found in both. yay! 
      if (dim(joined)[1] != dim(common)[1]) {noquote("Something is seriously wrong. Call Asli.")} # sanity check

      # GET TOTAL READ DEPTH FROM PILEUP
      pileup_name <- paste(strsplit(file, key)[[1]][[1]], "pileup", sep = ".") # this has the filtered annovar results that we will use to subset the mutect results
      pileup <- read.delim(paste0(path_to_pileup, type, "/RG_", pileup_name), header = FALSE)[, c(1, 2, 3, 4)]
      names(pileup) <- c("CHROM", "POS", "something", "TOTAL_DEPTH")
      pileup$something <- NULL # delete the col we don't need

      #second merge, merging with the total read depth at that position
      joined2 <- inner_join(joined, pileup, by = c("CHROM", "POS"))

      # calculate vaf 
      joined2 <- joined2 %>% 
        mutate(VAF = )







    
      }  

}