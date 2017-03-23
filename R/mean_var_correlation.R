library(data.table)
library(ggplot2)
library(hexbin)


message("Trio samples")

d = fread("data/counts/trio.100k.q10.sed.txt")
d <- d[grep('^chr[1-9X][0-9]?$', chrom),] # no chromosome Y
d[, c("chrom", "cell", "sample") := .( 
    factor(chrom),
    factor(substr(sample, nchar(sample)-2, nchar(sample))),
    factor(substr(sample, 1, nchar(sample)-4)) )  ]
setkey(d, sample, cell, chrom)

e <- d[, .(mean = mean(w+c), var = var(w+c), median = median(w+c)), by = .(sample, cell)]
e_len <- d[, .(len = length(unique(cell))), by = sample]
plt <- ggplot( e ) +
    aes(mean, var) + 
    geom_point() +
    geom_smooth(method="lm") +
    ggtitle("Trio samples") +
    facet_wrap(~sample, ncol = 4) + 
    scale_x_continuous( limits=c(0, quantile(e$mean,0.98))) +
    scale_y_continuous( limits=c(0, quantile(e$var,0.98))) +
    geom_label(data=e_len, x=0, y=Inf, hjust=0, vjust=1, aes(label = paste(len, "cells")))
ggsave(plt, filename = "data/plots/trio.100k.q10.mean_var_corr.pdf", width = 10, height=14)




message("Chromothripsis samples")

d = fread("data/counts/chromothripsis.500k.q10.sed.txt")
d <- d[grep('^[1-9X][0-9]?$', chrom),] # no chromosome Y
d[, c("chrom", "cell", "sample") := .( 
    factor(chrom),
    factor(substr(sample, nchar(sample)-2, nchar(sample))),
    factor(substr(sample, 1, nchar(sample)-4)) )  ]
setkey(d, sample, cell, chrom)

e <- d[, .(mean = mean(w+c), var = var(w+c), median = median(w+c)), by = .(sample, cell)]
e_len <- d[, .(len = length(unique(cell))), by = sample]
plt <- ggplot( e ) +
    aes(mean, var) + 
    geom_point() +
    geom_smooth(method="lm") +
    ggtitle("Chromothripsis samples") +
    facet_wrap(~sample, ncol = 4) + 
    scale_x_continuous( limits=c(0, quantile(e$mean,0.98))) +
    scale_y_continuous( limits=c(0, quantile(e$var,0.98))) +
    geom_label(data=e_len, x=0, y=Inf, hjust=0, vjust=1, aes(label = paste(len, "cells")))
ggsave(plt, filename = "data/plots/chromothripsis.500k.q10.mean_var_corr.pdf", width = 10, height=20)


