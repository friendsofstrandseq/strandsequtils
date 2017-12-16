library(tilingArray)
library(data.table)
library(tidyverse)
library(assertthat)


segmentation <- function(d, windowsize = 200e3, bp_per_Mb = 0.2) 
{
    assert_that(is.data.table(d))
    assert_that(all(c("chrom","start","end","sample","cell","c","w") %in% colnames(d)))
    n_col = nrow(unique(d[, .(sample,cell)]))
    
    BKP = list()
    for (chrom_ in unique(d$chrom))
    {
        d_ = d[chrom == chrom_,]
        e_ = melt(d_[,.(start,sample,cell,c,w)], measure.vars = c("c","w"), variable.name = "strand", value.name = "count")
        e_ = dcast.data.table(e_, start ~ sample + cell + strand, value.var = "count", sep = "___")
        Y = as.matrix(select(e_,-start))
        
        # Calculating max_k and max_segment
        chrom_size  = max(d_$end)
        max_length  = nrow(e_)
        max_number  = max(10,as.integer(floor(chrom_size/1e6 * bp_per_Mb))) # max number of segments
        message("Segmenting chrom ", chrom_, " across ", n_col, " cells. Max. ", max_number, " breakpoints.")
        
        # Run tiling array
        S = segment(Y, max_number, max_length)
        bkps = data.table(k   = rep(0:(max_number-1), 0:(max_number-1)), 
                          bkp = unlist(S@breakpoints) -1)
        BKP[[chrom_]] = bkps
    }
    return (BKP)
}



args = commandArgs(trailingOnly = T)

f_in = args[1]
windowsize = as.integer(args[2])
f_out = args[3]

d = fread(paste("zcat", f_in))
BKP = segmentation(d)
#plot_segmentation(d, BKP, file = "out.pdf", width = 16, height = 9)


BKP2 = NULL
for (chrom_ in names(BKP))
    BKP2 = rbind(BKP2, cbind(BKP[[chrom_]], chrom = chrom_))


message("Writing file ", f_out)
write.table(BKP2, file = f_out, quote=F, row.names=F, col.names=F, sep="\t")
