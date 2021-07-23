library(tidyverse)
library(ggrepel)

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

curated_df_path <- "/groups/wyattgrp/users/amunzur/chip_project/subsetted/new_finland_download/curated_muts.csv"
curated_df <- as.data.frame(read_csv(curated_df_path))

# adding new cols 
curated_df$sample <- unlist(lapply(as.list(curated_df$sample_n), function(some_name) strsplit(some_name, "-WBC")[[1]][[1]]))
curated_df$label <- paste(curated_df$CHROM, curated_df$POSITION, sep = "_")
names(curated_df)[c(4, 5)] <- c("Tumor", "WBC")
curated_df$label_long <- paste(curated_df$CHROM, curated_df$POSITION, curated_df$sample, sep = "_")

# label by chrom and position 
p_chrompos <- ggplot(data = curated_df, aes(x = Tumor, y = WBC, label = label)) + 
	geom_point() + 
	geom_smooth(method = "lm", se = FALSE) + 
	geom_text_repel() + 
	scale_x_continuous(breaks = seq(0, 0.18, by = 0.03)) + 
	scale_y_continuous(breaks = seq(0, 0.18, by = 0.03)) + 
	ggtitle("Manually curated list of mutations - Mutect2, second batch of samples") + 
	xlab("Tumor VAF") + 
	ylab("WBC VAF") + 
	theme_bw() + 
	cool_theme

# label by sample id
p_samples <- ggplot(data = curated_df, aes(x = Tumor, y = WBC, label = sample)) + 
	geom_point() + 
	geom_smooth(method = "lm", se = FALSE) + 
	geom_text_repel() + 
	scale_x_continuous(breaks = seq(0, 0.18, by = 0.03)) + 
	scale_y_continuous(breaks = seq(0, 0.18, by = 0.03)) + 
	ggtitle("Manually curated list of mutations - Mutect2, second batch of samples") + 
	xlab("Tumor VAF") + 
	ylab("WBC VAF") + 
	theme_bw() + 
	cool_theme

# classic plot
curated_df <- gather(curated_df, "VAF_type", "VAF_value", c("Tumor", "WBC")) # from wide to long
curated_df <- curated_df[order(curated_df$sample), ] # order by sample id
curated_df$label_long <- factor(curated_df$label_long, levels = unique(curated_df$label_long))

p_classic <- ggplot(data = curated_df, aes(x = VAF_value, y = label_long, color = VAF_type)) + 
	geom_point(size = 3) + 
	scale_x_continuous(breaks = seq(0, 0.18, 0.03)) + 
	geom_hline(yintercept = seq(0, length(unique(curated_df$label_long)), 1), linetype = "dotted", color = "gray") + 
	ggtitle("Manually curated list of mutations - Mutect2, second batch of samples") + 
	xlab("VAF") + 
	theme_bw() + 
	cool_theme

######################################################
# SAVE
######################################################
p_chrompos_path <- 
p_samples_path <- 
p_classic <- 