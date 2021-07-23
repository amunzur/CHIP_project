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
