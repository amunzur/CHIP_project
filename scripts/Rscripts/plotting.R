library(tidyverse)
library(ggrepel)
library(grid)
library(plotly)

GeomLabelRepel$draw_key <- function (data, params, size) { draw_key_rect(data) }

# to upload
path_to_curated_combined <-"/groups/wyattgrp/users/amunzur/chip_project/tnvstats_mutect_compared/combined.csv" # all muts are here, inbdicating whether they are found in only one pipeline or both

# to write 
path_to_figure <- "/groups/wyattgrp/users/amunzur/chip_project/figures/combined_fig.pdf"

df <- as.data.frame(read_csv(path_to_curated_combined)) # read the main data we will be plotting 
df <- transform(df, labelling = paste(CHROM, POSITION, sep=", ")) # combine two cols, they will be used to label the dots 

# non interactive plot
p1 <- ggplot(df, aes(x = VAF, y = WBC_VAF, color = STATUS)) + 
	geom_point(size = 2) + 
	geom_label_repel(aes(label = labelling),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
	theme_bw()

# interactive plot
p2 <- plot_ly(data = df,
        x = ~VAF, y = ~WBC_VAF,
        opacity = 1,
        color = ~STATUS,
        type = "scatter",
        mode = "markers",
        marker = list(size = 5), 
        text = labelling) %>% 
  layout(legend= list(itemsizing='constant'))


ggsave(path_to_figure, pl)

