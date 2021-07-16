# CHIP_project

bed_variants: contains essentially the annovar results, but after germline filtering. In tsv format. Only contains the common variants found in both sample groups.

variant_bed_files: contains the called variants in the bed file format. This was used to run mpileup on the bam files based on these locations only. 
The bed files here are generated from the tsv files in the bed_variant folder. Its name might also be common_variants! These bed files here only contain 
the variants called in the targeted genes panel. 

Bed files here are generated using the script here: /groups/wyattgrp/users/amunzur/chip_project/scripts/mpileup/run_mpileup.sh
The script above both generates the bed files we have here, and then run samtools mpileup based on them. 

## Sample filtering and variant calling
Scripts located in `scripts/Mutect2` filter the bam files, run Mutect2 on them and remove the false positives using GATK's functions. The filtering consists of removing duplicates, fixing read mates, adding read groups and removing low quality reads. More detail abo

## Calculating coverage 
We choose which samples to work on based on coverage. `scripts/compute_coverage/compute_coverage.sh` script computes the coverage using `samtools depth` for either all the bam files in a given directory, or in the files written down in a text file. The files should contain the full path. Low quality reads and duplicates should be removed prior to running calculating coverage, because otherwise the metrics will be inflated. The bam file list (in the form of a text file) should be located in the same place, `scripts/compute_coverage`, and it should be named `bamList_coverage.txt`. By default the script outputs coverage information to `metrics/coverage_information`. User should specify a subdirectory here for output. 