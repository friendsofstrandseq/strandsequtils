library(data.table)
library(ggplot2)
library(hexbin)

plotit <- function(file, out_file, min_reads = 5, limit = 70)
{
    d = fread(file)
    d <- d[grep('^(chr)?[1-9XY][0-9]?$', chrom),]
    d[, c("chrom", "cell", "sample") := .( 
        factor(chrom),
        factor(substr(sample, nchar(sample)-2, nchar(sample))),
        factor(substr(sample, 1, nchar(sample)-4)) )  ]
    setkey(d, sample, cell, chrom)
    
    cairo_pdf(out_file, width=10, height=10, onefile=T)
    for (s in unique(d$sample)) {
        message(s)
        e <- d[sample == s & w >= min_reads & c >= min_reads, ]
        # Print max 12 cells:
        e <- e[cell %in% unique(e$cell)[1:12], ]
        e$cell = factor(as.character(e$cell))
        corr <- e[, .(r = cor(w,c)), by = cell]
        plt <- ggplot(e) + aes(w,c) + 
            #geom_point(alpha=0.1, size=0.5) + # comment out if files get too big
            stat_density2d() +
            scale_x_continuous(limits = c(0,limit)) +
            scale_y_continuous(limits = c(0,limit)) +
            #coord_cartesian(xlim=c(0,limit), ylim=c(0,limit)) +
            facet_wrap(~cell) +
            ggtitle(paste0(s,", w>=", min_reads, " & c>=", min_reads)) + 
            geom_label(data = corr, x=Inf, y=Inf, hjust=1,vjust=1, aes(label = paste0("r^{2}==", round(r,3))), parse=T )
        print(plt)
    }
    dev.off()
}


plotit("data/counts/trio.100k.q10.sed.txt",
       "data/plots/trio.100k.q10.wc_corr.pdf")

plotit("data/counts/chromothripsis.500k.q10.sed.txt",
       "data/plots/chromothripsis.500k.q10.wc_corr.pdf",
       min_reads=5, limit=100)
