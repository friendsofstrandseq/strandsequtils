suppressWarnings(library(data.table))
library(ggplot2)
library(assertthat)

args = commandArgs(trailingOnly = T)
if (length(args) < 2 | length(args) > 3 || !file.exists(args[1]) || !grepl('.pdf$', args[2]) ) {
    message("Usage: Rscript cluster_segments.R mosaicatcher.out.txt out.pdf [out.clusterfile]")
    message("")
    message("Read a (non-gzipped) mosaicatcher output file and cluster cells by their strand states.")
    message("The order will be different for each chromosome. Creates a plot and can optionally save")
    message("this order in a file.")
    quit()
}
f_in = args[1]
f_pdf = args[2]

# Read all counts - remove "None" bins
e = fread(f_in, header=T)
e = e[class != "None",]

# Variables for output
clustering <- list()
clustering[[ "sample.cell.names" ]] = unique(paste(e$sample, e$cell, sep="_"))
plots <- list()

message("Finished reading ", length(clustering$sample.cell.names), " cells.")

for (chrom_ in unique(e$chrom)) 
{
    message("Clustering ", chrom_)
    
    # Subselect a chromosome
    e_ = e[chrom == chrom_,]
    N = length(unique(paste(e_$sample,e_$cell))) # number of cells
    
    # Get into wide format 
    e_cast <- dcast(e_, chrom+start+end ~ sample+cell, value.var = "class")
    assert_that(all(colnames(e_cast)[4:ncol(e_cast)] == clustering$sample.cell.names))
    
    # Calculate pairwise distance
    M = nrow(e_cast)
    distance = matrix(nrow = N, ncol = N)
    for (i in 1:N)
        for (j in 1:N)
            distance[i,j] = (M - sum(e_cast[,i+3,with=F] == e_cast[,j+3,with=F]))/M
    
    # Re-order cells
    order = hclust(as.dist(distance))$order
    clustering[[ chrom_ ]] = order
    e_$cell_name = factor(paste(e_$sample,e_$cell,sep="_"), levels = clustering$sample.cell.names[order], ordered = T)

    # Plot chromosome
    p <- ggplot(e_) +
        aes(xmin = start, xmax = end, ymin = as.numeric(cell_name), ymax = as.numeric(cell_name)+1, fill = class) +
        geom_rect() +
        scale_fill_brewer(type="qual", palette=3) +
        ggtitle(paste(chrom_, "across", N, "cells")) +
        theme_classic()
    plots <- list(plots, list(p))
}   

if (length(args)>2) {
    message("Saving cluster order to ", args[3])
    save(clustering, file = args[3])
}

message("Writing plot ...")
pdf(f_pdf, width=29.7, height = 21)
invisible(lapply(plots, print))
dev.off()

