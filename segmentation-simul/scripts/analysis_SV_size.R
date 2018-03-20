library(ggplot2)
source("scripts/analysis_common.R")

D <- read_files(snakemake@input)
p <- ggplot(D) + 
  aes(SV_found/SV_total, bp_true/bp_total, shape = Method, col = bp_per_Mb) + 
  ylab("Precision") +
  xlab("True positive rate (sensitivity)") +
  geom_line(alpha = 0.5, size = 0.5) +
  geom_point(alpha = 0.8, size = 1.5) +
  facet_grid(SV_size ~ SV_class) +
  scale_x_continuous(labels = scales::percent, breaks = c(0,0.5,1)) + 
  scale_y_continuous(labels = scales::percent, breaks = c(0,0.5,1)) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  guides(shape = guide_legend(override.aes = list(size=3))) +
  scale_color_gradient2(low = "dodgerblue3", 
                        mid = "darkorange", 
                        high = "red", 
                        midpoint = mean(unique(D$bp_per_Mb)),
                        name = "Breakpoints per Mb") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(p, filename = snakemake@output[[1]], width=10, height=6)
