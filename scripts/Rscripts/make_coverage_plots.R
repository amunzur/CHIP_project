library(tidyverse)
# This script makes coverage plots (histograms) for wbc and tumor samples.

# my very cool theme
cool_theme <- 
  
  theme(panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black", size = 1), 
      axis.ticks = element_line(colour = "black", size = 2),
      axis.text = element_text(size=10),
      axis.text.x = element_text(vjust=0.5, colour = "black", size=8),
      axis.text.y = element_text(vjust=0.5, colour = "black", size=6),
      axis.title = element_text(size=10,face="bold"), 
      legend.title = element_text(color = "black", size = 12),
      legend.text = element_text(color = "black", size = 12),
      axis.ticks.length=unit(0.15, "cm"), 
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 20)), 
      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

# INPUTS and OUTPUTS for the initial group of samples
# path_coverage_wbc <- "/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/GU_finland_download_72gene.csv"
# path_coverage_tumor <- "/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/ctDNA_prognosis_coverage_72gene.csv"
# path_wbc_plot <- "/groups/wyattgrp/users/amunzur/chip_project/figures/coverage_plots/first_batch/wbc_coverage.pdf"
# path_tumor_plot <- "/groups/wyattgrp/users/amunzur/chip_project/figures/coverage_plots/first_batch/tumor_coverage.pdf"

# INPUTS and OUTPUTS for the second group of samples
path_coverage_combined_FILTERED <- "/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/new_finland_download_FILTERED.csv"
path_wbc_plot_FILTERED <- "/groups/wyattgrp/users/amunzur/chip_project/figures/coverage_plots/second_batch/wbc_coverage_FILTERED.pdf"
path_tumor_plot_FILTERED <- "/groups/wyattgrp/users/amunzur/chip_project/figures/coverage_plots/second_batch/tumor_coverage_FILTERED.pdf"

path_coverage_combined_RAW <- "/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/new_finland_download_RAW.csv"
path_wbc_plot_RAW <- "/groups/wyattgrp/users/amunzur/chip_project/figures/coverage_plots/second_batch/wbc_coverage_RAW.pdf"
path_tumor_plot_RAW <- "/groups/wyattgrp/users/amunzur/chip_project/figures/coverage_plots/second_batch/tumor_coverage_RAW.pdf"

main <- function(path_coverage_wbc = NULL, 
				path_coverage_tumor = NULL,
				path_coverage_combined = NULL, 
				cool_theme, 
				SAVE, 
				what_are_you_plotting, 
				fill_color) {

	# what_are_you_plotting needs to be "first_batch" or "second_batch", given as strings.

	if (!is.null(path_coverage_wbc) & !is.null(path_coverage_tumor)) { # if they are given in separate files

		wbc <- as.data.frame(read.csv(path_coverage_wbc, header = FALSE))
		tumor <- as.data.frame(read.csv(path_coverage_tumor, header = FALSE))

		names(wbc) <- c("sample_path", "depth")
		names(tumor) <- c("sample_path", "depth")

		} else if (!is.null(path_coverage_combined)) {

			combined <- as.data.frame(read.csv(path_coverage_combined, header = FALSE))
			names(combined) <- c("sample_path", "depth")

			# separate into tumor and wbc samples 
			idx <- grepl("cfDNA|Baseline", combined$sample_path, ignore.case = TRUE)
			tumor <- combined[idx, ]

			idx <- grepl("WBC", combined$sample_path, ignore.case = TRUE)
			wbc <- combined[idx, ] } # end of if loop

	# now sort out a couple of misc stuff
	if (what_are_you_plotting == "initial_cohort"){

		plot_wbc <- ggplot(data = wbc, aes(x = depth)) +
			geom_hline(yintercept = seq(0, 20, 2), color = "black", linetype="dotted") + # horizontal grid lines
			geom_histogram(bins=25, fill="lightblue") + 
			geom_vline(aes(xintercept=mean(depth)), color="blue", linetype="dashed", size=1) +
			xlab("Average read depth (coverage)") + 
			ylab("Number of samples") +
			ggtitle("Histogram of average read depth in the first batch of CHIP samples - WBC", subtitle = "Dashed line indicates mean. 25 bins.") +
			scale_x_continuous(breaks = seq(0, 1000, 100)) + 
			scale_y_continuous(breaks = seq(0, 20, 2)) +
			theme_bw() + 
			cool_theme + 
			theme(axis.text.y = element_text(vjust=0.5, colour = "black", size=8))

		plot_tumor <- ggplot(data = tumor, aes(x = depth)) +
			geom_hline(yintercept = seq(0, 18, 2), color = "black", linetype="dotted") + # horizontal grid lines\ 
			geom_histogram(bins=25, fill="lightblue") + 
			geom_vline(aes(xintercept=mean(depth)), color="blue", linetype="dashed", size=1) +
			xlab("Average read depth (coverage)") + 
			ylab("Number of samples") +
			ggtitle("Histogram of average read depth in the first batch of CHIP samples - TUMOR", subtitle = "Dashed line indicates mean. 25 bins.") +
			scale_x_continuous(breaks = seq(0, 1600, 200)) + 
			scale_y_continuous(breaks = seq(0, 18, 2)) +
			theme_bw() + 
			cool_theme + 
			theme(axis.text.y = element_text(vjust=0.5, colour = "black", size=8)) 

	} else if (what_are_you_plotting == "second_batch") {

		plot_wbc <- ggplot(data = wbc, aes(x = depth)) +
			geom_hline(yintercept = seq(0, 80, 5), color = "black", linetype="dotted") + # horizontal grid lines\ 
			geom_histogram(bins=25, fill="lightblue") + 
			geom_vline(aes(xintercept=mean(depth)), color="blue", linetype="dashed", size=1) +
			xlab("Average read depth (coverage)") + 
			ylab("Number of samples") +
			ggtitle("Histogram of average read depth in the second batch of filtered CHIP samples - WBC", subtitle = "Dashed line indicates mean. 25 bins.") +
			scale_x_continuous(breaks = seq(0, 3000, 500)) + 
			scale_y_continuous(breaks = seq(0, 80, 5)) +
			theme_bw() + 
			cool_theme + 
			theme(axis.text.y = element_text(vjust=0.5, colour = "black", size=8))

		plot_tumor <- ggplot(data = tumor, aes(x = depth)) +
			geom_hline(yintercept = seq(0, 80, 5), color = "black", linetype="dotted") + # horizontal grid lines\ 
			geom_histogram(bins=25, fill="lightblue") + 
			geom_vline(aes(xintercept=mean(depth)), color="blue", linetype="dashed", size=1) +
			xlab("Average read depth (coverage)") + 
			ylab("Number of samples") +
			ggtitle("Histogram of average read depth in the second batch of filtered CHIP samples - TUMOR", subtitle = "Dashed line indicates mean. 25 bins.") +
			scale_x_continuous(breaks = seq(0, 4000, 500)) + 
			scale_y_continuous(breaks = seq(0, 80, 5)) +
			theme_bw() + 
			cool_theme + 
			theme(axis.text.y = element_text(vjust=0.5, colour = "black", size=8))} # end of plotting if loop

		return(c(list(plot_tumor), list(plot_wbc)))

	} # end of function

##########################################
# PLOTTING
##########################################
# to plot the data from 72 gene panel, the initial group of samples we have been working on 
# main(path_coverage_wbc = path_coverage_wbc, 
# 	 path_coverage_tumor = path_coverage_tumor,
# 	 path_coverage_combined = NULL, 
# 	 path_wbc_plot = path_wbc_plot, 
# 	 path_tumor_plot = path_tumor_plot, 
# 	 cool_theme =  cool_theme, 
# 	 SAVE = FALSE, 
# 	 what_are_you_plotting = "first_batch")

PLOT_TUMOR <- main(path_coverage_wbc = NULL, 
	 path_coverage_tumor = NULL,
	 path_coverage_combined = path_coverage_combined_FILTERED, 
	 cool_theme = cool_theme, 
	 SAVE = TRUE, 
	 what_are_you_plotting = "second_batch")[[1]]

PLOT_WBC <- main(path_coverage_wbc = NULL, 
	 path_coverage_tumor = NULL,
	 path_coverage_combined = path_coverage_combined_FILTERED, 
	 cool_theme = cool_theme, 
	 SAVE = TRUE, 
	 what_are_you_plotting = "second_batch")[[2]]

ggsave(path_wbc_plot_FILTERED, PLOT_WBC, width = 12, height = 12, units = "in")
ggsave(path_tumor_plot_FILTERED, PLOT_TUMOR, width = 12, height = 12, units = "in")

PLOT_TUMOR <- main(path_coverage_wbc = NULL, 
	 path_coverage_tumor = NULL,
	 path_coverage_combined = path_coverage_combined_RAW, 
	 cool_theme = cool_theme, 
	 SAVE = TRUE, 
	 what_are_you_plotting = "second_batch")[[1]]

PLOT_WBC <- main(path_coverage_wbc = NULL, 
	 path_coverage_tumor = NULL,
	 path_coverage_combined = path_coverage_combined_RAW, 
	 cool_theme = cool_theme, 
	 SAVE = TRUE, 
	 what_are_you_plotting = "second_batch")[[2]]

ggsave(path_wbc_plot_RAW, PLOT_WBC, width = 12, height = 12, units = "in")
ggsave(path_tumor_plot_RAW, PLOT_TUMOR, width = 12, height = 12, units = "in")