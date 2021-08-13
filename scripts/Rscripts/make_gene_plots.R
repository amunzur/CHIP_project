library(tidyverse)

	cool_theme <-

	  theme(panel.border = element_blank(),
	      panel.grid.major = element_blank(),
	      panel.grid.minor = element_blank(),
	      axis.line = element_line(colour = "black", size = 1),
	      axis.ticks = element_line(colour = "black", size = 2),
	      axis.text = element_text(size=10),
	      axis.text.x = element_text(vjust=0.5, colour = "black", size=12),
	      axis.text.y = element_text(vjust=0.5, colour = "black", size=12),
	      axis.title = element_text(size=14, face="bold"),
	      legend.title = element_text(color = "black", size = 12),
	      legend.text = element_text(color = "black", size = 12),
	      axis.ticks.length=unit(0.15, "cm"),
	      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 20)),
	      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

cohort <- "kidney_samples_second"

# Code below cleans up the df and saves it, makes it ready to plot the coverage plots. 
path_to_df <- file.path("/groups/wyattgrp/users/amunzur/chip_project/subsetted", cohort, "curated_muts_panel_annotated_germFiltered.csv")
df <- as.data.frame(read_csv(path_to_df))

df$sample <- unlist(lapply(as.list(df$SAMPLE_ID), function(some_name) strsplit(some_name, "-cfDNA")[[1]][[1]]))
df$label_long <- paste(df$CHROM, df$POSITION, df$sample, sep = "_")
df$label_long <- factor(df$label_long, levels = unique(df$label_long))

paths <- list.files(file.path("/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/individual_samples", "kidney_samples"), full.names = TRUE)
depth_list_cfDNA <- list()
depth_list_WBC <- list()

for (row in 1:nrow(df)) {
	message(df[row, ]$sample)

    # load the related coverage df
    coverage_df_cfDNA <- as.data.frame(read.delim(paths[which(grepl(df[row, ]$sample, paths) & grepl("cfDNA", paths))], header = FALSE))
    coverage_df_WBC <- as.data.frame(read.delim(paths[which(grepl(df[row, ]$sample, paths) & grepl("WBC", paths))], header = FALSE))

    names(coverage_df_cfDNA) <- names(coverage_df_WBC) <- c("CHROM", "POSITION", "DEPTH")

    # filter for the position and chrom
    depth_list_cfDNA <- append(depth_list_cfDNA, filter(coverage_df_cfDNA, CHROM == df[row, ]$CHROM, POSITION == df[row, ]$POSITION)[, 3])
    depth_list_WBC <- append(depth_list_WBC, filter(coverage_df_WBC, CHROM == df[row, ]$CHROM, POSITION == df[row, ]$POSITION)[, 3]) } # end of for loop

# add depth info to df
df$depth_cfDNA <- as.numeric(unlist(depth_list_cfDNA))
df$depth_WBC <- as.numeric(unlist(depth_list_WBC))

write_csv(df, "/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/misc_dfs/kidney_samples/make_gene_plot_df.csv")

##########################################################
# PLOTTING FUNCTION
##########################################################
make_plots <- function(cohort, depth_df_path, chosen_theme, plot_title_gene, plot_title_vaf, plot_title_depth, fig_path_main) {

	path_to_df <- file.path("/groups/wyattgrp/users/amunzur/chip_project/subsetted", cohort, "curated_muts_panel_annotated_germFiltered.csv")
	df <- as.data.frame(read_csv(path_to_df))

	# 1. GENE PLOT ##############################################################
	p_gene <- ggplot(data = df, aes(x = Gene.refGene)) + 
			geom_histogram(stat="count") +
			xlab("Gene") + 
			ylab("Counts") +
			ggtitle(plot_title_gene) +
			theme_bw() + 
			cool_theme + 
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

	# vaf plot
	df_vaf <- gather(df, "VAF_type", "VAF_value", c("VAF", "WBC_VAF")) # from wide to long
	df_vaf <- df_vaf[order(df_vaf$sample), ] # order by sample id
	df_vaf$VAF_type <- ifelse(df_vaf$VAF_type == "VAF", "Tumor", "WBC") # minor renaming 

	df_vaf$sample <- unlist(lapply(as.list(df_vaf$SAMPLE_ID), function(some_name) strsplit(some_name, "-cfDNA")[[1]][[1]]))
	df_vaf$label_long <- paste(df_vaf$CHROM, df_vaf$POSITION, df_vaf$sample, sep = "_")
	df_vaf$label_long <- factor(df_vaf$label_long, levels = unique(df_vaf$label_long))

	# 2. VAF PLOT ##############################################################
	p_vaf <- ggplot(data = df_vaf, aes(x = VAF_value, y = label_long, color = VAF_type)) + 
			geom_point(size = 3) + 
			scale_x_continuous(breaks = seq(0, 0.20, 0.02)) + 
			geom_hline(yintercept = seq(0, length(unique(df_vaf$label_long)), 1), linetype = "dotted", color = "gray") + 
			geom_vline(xintercept = seq(0, 0.20, 0.02), linetype = "dotted", color = "gray") + 
			ggtitle(plot_title_vaf) + 
			ylab("") +
			xlab("VAF") + 
			theme_bw() + 
			cool_theme

	p_vaf2 <- ggplot(data = df, aes(x = VAF, y = WBC_VAF)) + 
			geom_point(size = 1) + 
			geom_smooth(method = "lm", se = FALSE) + 
			scale_x_continuous(breaks = seq(0, 0.20, 0.02)) + 
			scale_y_continuous(breaks = seq(0, 0.20, 0.02)) + 
			geom_hline(yintercept = seq(0, 0.20, 0.02), linetype = "dotted", color = "gray") + 
			geom_vline(xintercept = seq(0, 0.20, 0.02), linetype = "dotted", color = "gray") + 
			ggtitle(plot_title_vaf) + 
			ylab("WBC VAF") +
			xlab("Tumor VAF") + 
			theme_bw() + 
			cool_theme


	# 3. COVERAGE PLOT ##############################################################
	df <- as.data.frame(read_csv(depth_df_path))
	df <- gather(df, "sample_type", "depth", c("depth_cfDNA", "depth_WBC")) # from wide to long
	df$sample_type <- ifelse(df$sample_type == "depth_cfDNA", "Tumor", "WBC")

	# plotting
	p_coverage <- ggplot(data = df, aes(x = label_long, y = depth, fill = sample_type)) + 
				geom_bar(position="dodge", stat="identity") + 
				scale_y_continuous(breaks = seq(0, 6000, 1000), expand = c(0, 0)) + 
				geom_hline(yintercept = seq(0, 6000, 1000), linetype = "dotted", color = "gray") + 
				# geom_vline(xintercept = seq(0, 0.20, 0.02), linetype = "dotted", color = "gray") + 
				ggtitle(plot_title_depth) + 
				ylab("Depth") +
				xlab("") + 
				theme_bw() + 
				cool_theme + 
				coord_flip()

	# 4. SAVE THE PLOTS ##############################################################
	cohort <- "kidney_samples"
	output_gene <- file.path(fig_path_main, cohort, "genes.pdf")
	output_vaf <- file.path(fig_path_main, cohort, "vaf.pdf")
	output_vaf2 <- file.path(fig_path_main, cohort, "vaf2.pdf")
	output_depth <- file.path(fig_path_main, cohort, "depth.pdf")

	ggsave(output_gene, p_gene, width = 10, height = 10, units = "cm")
	ggsave(output_vaf, p_vaf, width = 10, height = 10, units = "cm")
	ggsave(output_vaf2, p_vaf, width = 10, height = 10, units = "cm")
	ggsave(output_depth, p_coverage, width = 10, height = 10, units = "cm")

} # end of function



make_plots <- function(cohort = "kidney_samples_second", 
						depth_df_path = "/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/misc_dfs/kidney_samples/make_gene_plot_df.csv", 
						chosen_theme = cool_theme, 
						plot_title_gene = "Genes where mutations appear in kidney samples", 
						plot_title_vaf = "Tumor and WBC VAF in kidney samples", 
						plot_title_depth = "Depth at SNVs", 
						fig_path_main = "/groups/wyattgrp/users/amunzur/chip_project/figures/main_figures")


curated_df <- gather(curated_df, "Tumor_or_normal", "Sample_name", c("SAMPLE_ID", "sample_n")) # from wide to long










