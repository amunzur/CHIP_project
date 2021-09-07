# Looking for CHIP mutations in cfDNA samples

## Sample filtering and variant calling
Scripts located in `scripts/Mutect2` filter the bam files, run Mutect2 on them and remove the false positives using GATK's functions. The filtering consists of removing duplicates, fixing read mates, adding read groups and removing low quality reads. After running and filtering Mutect2, the nest section of R scripts clean up the results and compare the found variants to the results from other pipelines. 

## Filtering called variants and comparing to tnvstats 
Scripts in this section are found in `scripts/Rscripts`. The outputs of the scripts along with a summary of what they do is given below.

1. **`prepare_bamslist_for_tnvstats.R`**: This unnessarily complicated script prepares a list of bam files in a text file in order to run tnvstats. The important thing is to make sure that the tumor and wbc ids match in a given row. It outputs the text file to `finland_bams/bamslist/tnvstats_bamList.csv` and overwrites each time the script is run. 

2. **`find_common_variants.R`**: As the first step after running and filtering Mutect, this script identifies the common variants between the wbc and tumor samples, and saves the output to separate csv files for each samples in `common_variants/`, since the mutations that are likely to be CHIP must exist in both WBC and tumor samples. The script will only consider samples with both tumor and wbc vcf files.

3. **`get_vaf_from_tnvstats.R`**: For each sample, this script cross-checks the information from Mutect2 with tnvstats. Mainly, for each identified variant, it identifies the tumor vaf & read support, wbc vaf & read support. Based on the thresholds given in the script, it then filters the mutect results to retain the significant variants only. The script assumes that the common variants between WBC and tumor samples in mutect results have been identified. 

	- `subsetted/sample_name.csv`: A csv file with information from both tnvstats and the Mutect2. It contains vaf, read support, ref and alt information about all variants identified in a sample. 
	- `subsetted/ALL_MERGED.csv`: All mutations identified in all samples by Mutect2 before any filtering. 
	- `subsetted/ALL_MERGED_FILTERED_0.1.csv`: All mutations identified in all samples by Mutect2 with vaf less than 10%. 
	- `subsetted/ALL_MERGED_FILTERED_0.2.csv`: All mutations identified in all samples by Mutect2 with vaf less than 20%. 

Furthermore, the script also checks whether the tnvstats and the common variants csv files exist for each sample and outputs a metrics file thatis useful for troubleshooting later, located here: `chip_project/metrics/file_situation/new_finland_download`.

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

## Calculating metrics
Besides mutation calling, it's also possible to compute various metrics, such as sample coverage (used synonymously as sample depth) and the number of variants mutect identified

1. **COVERAGE:** We choose which samples to work on based on coverage. `scripts/compute_coverage/compute_coverage.sh` script computes the coverage using `samtools depth` for either all the bam files in a given directory, or in the files written down in a text file. The files should contain the full path. Low quality reads and duplicates should be removed prior to running calculating coverage, because otherwise the metrics will be inflated. The bam file list (in the form of a text file) should be located in the same place, `scripts/compute_coverage`, and it should be named `bamList_coverage.txt`. By default the script outputs coverage information to `metrics/coverage_information`. User should specify a subdirectory here for output. In the output directory, the files with the "ALL" string in the name contains the average depth computed with all the positions in the bam file, not subsetted to the 72 gene panel file.  

Besides calculating coverage, we also make plots to show the distribution of coverage across samples: 
**`make_coverage_plots.R`**: This script makes histograms showing the distribution of coverage in wbc and tumor samples. The outputted figures are located at `figures/coverage_plots/`.

It is a good idea to add the exact coverage information for the SNVs to the sheet we curated (`curated_muts.csv`). `add_coverage_mq_info.R` does exactly that and saves a new copy of the `curated_muts.csv` file to the same location with two extra columns added for coverage and mapping quality.

2. **NUMBER OF VARIANTS:** Shell script located at `chip_project/scripts/Mutect2/compute_variant_numbers.sh` calculates the number of variants Mutect2 called, after filtering for false positives. The output is a text file consisting of 2 columns where the first one is the vcf file, and the second one is the number of mutations mutect called. Usually WBC and tumor samples appear one after another on the file.

The script also calculates the number of common variants found across WBC and tumor samples. As an input it uses the csv files located at `chip_project/common_variants/new_finland_download`, and the metrics are saved at `chip_project/metrics/mutect_variant_numbers/new_finland_download_COMMON`. 

## About running tnvstats on the samples 
Tnvstats scripts aren't located in this in this project. They are run separately, but in each tnvstats directory for each batch of samples, bamlist used is found in the same place.

## Information about bam files 
As of July 16, the analysis encompasses three groups of data: 
1. **ctDNA_prognosis**: Tumor samples from 72 gene targeted sequencing data 
2. **GU_finland_download**: WBC samples from 72 gene targeted sequencing data - matches to data in step 1 
3. **new_finland_download**: WBC and tumor samples are together in this group, containing bam files from OZM, GU, VPC, VCC, M1RP cohorts. Bam files starting with the `filtered_RG` have been filtered to remove duplicates and low quality reads. Filtered Mutect2 results for these files are located in `chip_project/mutect_results_filtered/new_finland_download`. 

And these are the misc data, relating to the main 3 groups mentioned above:  
4. **new_finland_download_UNCOMP**: Bam files here have been uncompressed after filtering and adding the readgroups. The reason was to see if computing tnvstats on uncompressed bams would be faster.

## IGV snapshots
Taking IGV snapshots of the candidate CHIP mutatinos helps validate our findings. The snapshots scripts are located at `gene_panel_pipeline_modified/igv_scripts`, not in the `chip_project` directory. `make_igv_batch_file_general.py` is the main script that can run on various types of inputs. It creates a batch file that should be run on a linux or Windows system.

## Plotting 
There are multiple R scripts present that generate various plots that are saved to `chip_project/figures`, as follow: 
1. `make_mutect_plots.R`: This one generates 3 plots and outputs to `chip_project/figures/vaf_comparison/`. The plots it generates are: 
	- Classic vaf comparison plot1: A dot plot comparing the WBC and tumor VAFs of mutations. Dots are labelled to show chromosome and position. 
	- Classic vaf comparison plot2: A dot plot comparing the WBC and tumor VAFs of mutations. Dots are labelled to show samples. 
	- Classic vaf comparison plot3: A dot plot, x axis has the mutations and tumor & WBC vafs are colored differently. 

2. `make_coverage_plots.R`: This script plots the average coverage (depth) of samples in histograms and outputs the figures to `chip_project/figures/coverage_plots`. These plots help decide which cohorts to work with.

## Running the inhouse pipeline for SNV and indel calling
Besides running public software we also use a mutation calling software optimized for the types of data we work with. This pipeline consists of multiple scripts that need to be run in order, as explained below. All scripts are located at `chip_project/inhouse_pipeline/src`. 

- `pipeline_mutation_calling_all.bash`: A bash script to call `pipeline_mutation_calling_one.bash` on each tumor-normal samples, submits individual jobs to the cluster for each pair. Creates a folder for each tumour sample and stores all files in `chip_project/inhouse_pipeline/results_kidney_samples/${name_tumor}_bg_error_${bm}/variant_intermediate_files`. This output directory can be changed. 

For targeted data, use `genepanel_allchr_errorrates.txt` for background error rates. For whole exome sequencing, `hg38_wxs_allchr_18may18_errorrates.txt` should be used instead. 

- `pipeline_mutation_calling_one.bash`: This one calls snvs, and then applies low complexity, read end filers, and filters out potential gemline mutations. As an input, it uses the tnvstats and the tumor - WBC bam pairs.  

Process: 
#### SNV calling

1. Load BAM files and tnvstats
2. Run `call_snv.py` to query SNVs and annotate with ANNOVAR. Outputs to `${output}_both_somatic_snv_anno_flagged.tsv`. Some important parameters are: 
	- Upper VAF bound for a region to be classified as a germline heterzygous SNP &rightarrow; 0.8. 
	- Min read support from tumour &rightarrow; 10 (need to adjust base on cfDNA depth, the current one is assuming ~750X or more)
	- Min read support from WBC &rightarrow; 8 (need to adjust base on cfDNA depth, the current one is assuming ~500X or more. For kidney samples around 100X WBC, use 2 reads.)
	- Background error rate &rightarrow; 10
	- Min VAF &rightarrow; 1% (though some manual tweaking later is acceptable)
	Note that in the output file, it will be specified which samples are don't have enough read depth so we can go back and make our best judegement of including it or not. 

3. Apply low complexity filter on SNVs by calling `lowcomplex_filter_snv.py`. Outputs to `${output}_snv_lcf.tsv`.
4. Apply read end filter on SNVs by calling `read_end_filter_snv.py`. Outputs to `${output}_snv_lcf_re.tsv`.
5. Further filtering to   
	- Remove potential germline mutations based on database search (ExAC, and anything above 30% VAF).
	- Remove positions called in more than 2 samples.
	- Ensure VAF in WBC and cfDNA differ from each other by more than 5%.
	Outputs to `${output}_snv_lcf_re_pgf.tsv`.

#### Indel calling

6. Run `call_indel.py` to query INDELS and annotate with ANNOVAR. Outputs to `${output}_somatic_anno_indel.tsv`.
7. Apply low complexity filter on INDELS by calling `lowcomplex_filter_indel.py`. Outputs to: `${output}_indel_lcf.tsv`.
8. Further filtering to   
	- Remove potential germline mutations based on database search (ExAC, and anything above 30% VAF).
	- Remove positions called in more than 2 samples.
	- Ensure VAF in WBC and cfDNA differ from each other by more than 5%.
	Outputs to `${output}_indel_lcf_pgf.tsv`.

##### Abbreviations
- lcf: low complexity filter
- re: read end filter
- pgf: potential germline filter

#### Tabulating the output 
9. Loop through all samples and tabulate `${output}_snv_lcf_re_pgf.tsv` and `${output}_indel_lcf_pgf.tsv` to subset to genes of interest: ATM, AXSL1, TP53, BRCA2, BRCA1, RUNX1, GNAS, KMT2C (this is for the PCa panel) and restrict to coding variants that are more likely to be functional.

`Combine_CHIP_muts.py` for SNV  
`Combine_CHIP_indels.py` for INDELs

10. Take IGV snapshots and manuall curate variants. 
11. Produce visualizations:
	- WBC VAF and cfDNA VAF scatter plot
	- Per gene bar plot
	- `Combine_first_second_results.py`: WBC and cfDNA per mutation VAF scatter plot combined with bar plot showing read depth at the mut or combined with bar plot showing tumor content estimate.

### Result interpretation:
- WBC and tumor VAFs should follow the `y = x` line where they are rougly equal.
- Mutations should be in the coding region.
- Can cross reference the data on cBioportal. Two datasets of interest are "curated set of non-redundant studies" and "Cancer Therapy and Clonal Hematopoiesis (MSK, Nat Gent 2020)"" to see if the mutations called are recorded as potential oncogenetic in the literature.
- Other mutation functional predictor software includes polyfen and sift.

### TODOs
1. Run the bladder samples as Alex suggested
2. Call cfDNA muts of our PCa samples by running in house pipeline then run TC estimate script get a rough estimate of TC.
3. The read end filter is not working right now since we don't take into account of small indels' impact on read position. 



 


