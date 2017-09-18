
d = fread("simulation.input.txt")
ggsave(ggplot(d) + 
           aes(V3-V2) + 
           geom_histogram(binwidth = 20000) + 
           scale_x_continuous(breaks=pretty_breaks(10), label=comma) + 
           xlab("SV size"), 
       filename = "simulation.input.SVsize.pdf",
       width=6,
       height=4)

ggsave(ggplot(d) + 
           aes(V1) + 
           geom_bar() + 
           coord_flip() + 
           xlab(NULL) + 
           scale_y_continuous(breaks = pretty_breaks(8)),
       filename = "simulation.input.chromosomes.pdf",
       width=4,
       height=6)

