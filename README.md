# CHIP_project


## Sample filtering and variant calling
Scripts located in `scripts/Mutect2` filter the bam files, run Mutect2 on them and remove the false positives using GATK's functions. The filtering consists of removing duplicates, fixing read mates, adding read groups and removing low quality reads. After running and filtering Mutect2, the nest section of R scripts clean up the results and compare the found variants to the results from other pipelines. 

## Filtering called variants and comparing to tnvstats 
Scripts in this section are found in `scripts/Rscripts`. The outputs of the scripts along with a summary of what they do is given below.

1. **`prepare_bamslist_for_tnvstats.R`**: This unnessarily complicated script prepares a list of bam files in a text file in order to run tnvstats. The important thing is to make sure that the tumor and wbc ids match in a given row. It outputs the text file to `finland_bams/bamslist/tnvstats_bamList.csv` and overwrites each time the script is run. 

2. **`find_common_variants.R`**: As the first step after running and filtering Mutect, this script identifies the common variants between the wbc and tumor samples, and saves the output to separate csv files for each samples in `common_variants/`, since the mutations that are likely to be CHIP must exist in both WBC and tumor samples.

3. **`get_vaf_from_tnvstats.R`**: For each sample, this script cross-checks the information from Mutect2 with tnvstats. Mainly, for each identified variant, it identifies the tumor vaf & read support, wbc vaf & read support. Based on the thresholds given in the script, it then filters the mutect results to retain the significant variants only.

	- `subsetted/sample_name.csv`: A csv file with information from both tnvstats and the Mutect2. It contains vaf, read support, ref and alt information about all variants identified in a sample. 
	- `subsetted/ALL_MERGED.csv`: All mutations identified in all samples by Mutect2 before any filtering. 
	- `subsetted/ALL_MERGED_FILTERED_0.1.csv`: All mutations identified in all samples by Mutect2 with vaf less than 10%. 
	- `subsetted/ALL_MERGED_FILTERED_0.2.csv`: All mutations identified in all samples by Mutect2 with vaf less than 20%. 

**`compute_vaf.R`**: The purpose of this script was to calculate the vaf variants based on samtools and Mutect2 calculations, but instead decided to use tnvstats entirely. Keeping it here in case we need to refer back to it. For now, **`get_vaf_from_tnvstats.R`** does the job instead. 

4. **`curate_mutect.R`**: This script curates the mutect outputs based on manual curation and insights we gain from tnvstats. So essentially this script introduces one more level of filtering after step 1. Before running this, we need to manually go through all the IGV snapshots for the mutect mutations and save the ones we think are real somewhere. Then this script will go through those png files and filter. 

	- `subsetted/curated_muts.csv`: Main output of the script, contains the curated list of mutations afrter scrutinous curation and filtering.
	- `tnvstats_mutect_compared/both.csv`: Mutations identified by both pipelines
	- `tnvstats_mutect_compared/mutect_only.csv`: Mutations identified by mutect only, may contain duplicate mutations that occurred in more than one patient.
	- `tnvstats_mutect_compared/mutect_only_filtered.csv`: Almost the same as above, but this one has unique mutations only, no duplicates. 
	- `tnvstats_mutect_compared/elie_only.csv`: Mutations identified by in-house pipeline only
	- `tnvstats_mutect_compared/combined.csv`: All muts are here, indicating whether they are found in only one pipeline or both. Good concise file for referencing later on. 

5. **`look_for_missing_variants.R`**: IN PROGRESS. This script goes back to tnvstats to understand why the in-house pipeline failed to detect some variants. It needs `tnvstats_mutect_compared/combined.csv` as input.

6. **`match_variants_with_annovar.R`**: This step can be skipped if we decide to not use annovar for annotation. This script will match the identified variants with annovar to annotate them. Variants are filtered based on two criteria here: 
	- They are removed if annovar identifies them as germ line mutations. 
	- They are removed if they occur outside of the 72 gene panel. To do so, we filter based on the bed file. 
	- `bed_variants`: Outputs for each sample are saved here as a separate tsv file. 

7. **`plotting.R`**: IN PROGRESS. This script plots the mutations taken from `tnvstats_mutect_compared/combined.csv` (step 5 above). It also indicates if a mutation was identified by both pipelines or by one of them only. 

## Calculating coverage 
We choose which samples to work on based on coverage. `scripts/compute_coverage/compute_coverage.sh` script computes the coverage using `samtools depth` for either all the bam files in a given directory, or in the files written down in a text file. The files should contain the full path. Low quality reads and duplicates should be removed prior to running calculating coverage, because otherwise the metrics will be inflated. The bam file list (in the form of a text file) should be located in the same place, `scripts/compute_coverage`, and it should be named `bamList_coverage.txt`. By default the script outputs coverage information to `metrics/coverage_information`. User should specify a subdirectory here for output. In the output directory, the files with the "ALL" string in the name contains the average depth computed with all the positions in the bam file, not subsetted to the 72 gene panel file.  

Besides calculating coverage, we also make plots to show the distribution of coverage across samples: 
8. **`make_coverage_plots.R`**: This script makes histograms showing the distribution of coverage in wbc and tumor samples. The outputted figures are located at `figures/coverage_plots/`.

## Information about bam files 
As of July 16, the analysis encompasses three groups of data: 
1. **ctDNA_prognosis**: Tumor samples from 72 gene targeted sequencing data 
2. **GU_finland_download**: WBC samples from 72 gene targeted sequencing data - matches to data in step 1 
3. **new_finland_download**: WBC and tumor samples are together in this group, containing bam files from OZM, GU, VPC, VCC, M1RP cohorts. Bam files starting with the `filtered_RG` have been filtered to remove duplicates and low quality reads. Filtered Mutect2 results for these files are located in `chip_project/mutect_results_filtered/new_finland_download`. 