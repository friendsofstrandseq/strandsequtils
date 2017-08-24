library(data.table)
library(tidyverse)
library(assertthat)
library(jointseg)

log(1.5, base=10)

segmentation <- function(d, bp_per_Mb = 0.2, normalize = T)
{
    assert_that(is.data.table(d))
    assert_that(all(c("chrom","start","end","sample","cell","c","w") %in% colnames(d)))
    n_col = nrow(unique(d[, .(sample,cell)]))
    
    BKP = list()
    for (chrom_ in unique(d$chrom))
    {
        d_ = d[chrom == chrom_,]
        
        # Prepare Y matrix
        e_c = dcast(d_, chrom + start + end ~ sample + cell, value.var = "c", sep = "___")
        e_w = dcast(d_, chrom + start + end ~ sample + cell, value.var = "w", sep = "___")
        
        #  Generate matrix of total coverage
        Y_total = e_c[,4:ncol(e_c),with=F] + e_w[,4:ncol(e_w),with=F]
        
        colsd   = sapply(Y_total, sd)
        
        #  Y consists of Crick columns + Watson columns
        Y = cbind(e_c[,4:ncol(e_c),with=F], e_w[,4:ncol(e_c),with=F])
        colnames(Y) = c(paste0("c_", colnames(e_c)[4:ncol(e_c)]),
                        paste0("w_", colnames(e_w)[4:ncol(e_w)]) )
        
        # Scale by standard deviation of total coverage per sample/cell
        if (normalize) {
            Y = as.data.table(as.matrix(Y) / as.numeric(colsd))
        } else {
            Y = as.data.table(as.matrix(Y))
        }
        
        # Calculating max_k and max_segment
        chrom_size  = max(d_$end)
        max_length  = nrow(e_c)
        max_number  = max(10,as.integer(floor(chrom_size/1e6 * bp_per_Mb))) # max number of segments
        
        message("Segmenting chrom ", chrom_, " across ", n_col, " cells. Max. ", max_number, " breakpoints.")
        
        # Run jointSeg
        S = jointseg::jointSeg(Y, K = max_number, method = "GFLars", verbose = F)
        bkps = data.table(k   = rep(0:(max_number), 0:(max_number)),
                          bkp = unlist(unlist(S$dpBkpList)) -1)
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
