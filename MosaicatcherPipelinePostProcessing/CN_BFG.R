## Load required libraries
library(data.table)
library(GenomicAlignments)

bin.size <- 100000
maximumCN <- 200
confidence.interval = 0.95
sample.id <- "C7_data"
#regions <- "/media/porubsky/Elements/StrandSeqNation/C7/BFB_CN_estimates/BFBregion_C7_CNV_segments.txt"
dir <- "/home/maryam/research/Mosaicatcher/BFB_analysis/C7-7Jun"
outputfolder <- dir
#bam.folder <- "/media/porubsky/Elements/StrandSeqNation/C7/"

# get the file names and the directory paths of input files (Mosaicatcher outputs)
countsFile <- file.path(dir, paste0(format(bin.size, scientific = F),"_fixed.txt.gz"))
infoFile <- file.path(dir, paste0(format(bin.size, scientific = F),"_fixed.info"))
stateFile <- file.path(dir, "final.txt")

# read the input files
counts <- fread(paste("zcat", countsFile))
#segs <- getCountsPerRegions(bam.folder=bam.folder, outputfolder=outputfolder, regions=regions, sample.id=sample.id)
#segs <- as(segs, 'data.table')
info <- fread(infoFile)
strand <- fread(stateFile)
segs <- fread("/home/maryam/research/Mosaicatcher/BFB_analysis/C7_data.RegionCounts.txt")

# source mosaiclassifier
file.sources = list.files(file.path(dir, "mosaiClassifier/"), pattern="*.R", full.names = TRUE)
for (file in file.sources) {
  source(file)
}

get_CN_NB_prob <- function(counts, segs, info, strand, manual.segs=TRUE) {
  
  # set the None classes to WW in the counts data table
  counts[class=="None", class:="WW"]
  # prepare the data to be used in SV classification
  probs <- mosaiClassifierPrepare(counts, info, strand, segs, manual.segs = manual.segs)
  probs[, nb_r:=expected*nb_p/(1-nb_p)]
  
  probs <- probs[, cbind(.SD, CN=1:maximumCN), 
                 by=.(sample, cell, chrom, start, end)]
  
  
  CNprobs <- probs[, .(cov=sum(W+C), 
                       expected=sum(expected), 
                       nb_p=nb_p[1], 
                       nb_r=sum(nb_r)), 
                   by=.(sample, chrom, start, end)]
  # first observation
  (CNprobs$cov / CNprobs$expected) * 2
  
  CNprobs <- CNprobs[, cbind(.SD, CN=1:maximumCN), 
                     by=.(sample, chrom, start, end)]
  
  
  # computing NB probs and calling most probable CN in probs table per segment per cell
  probs[, CN_ll:=dnbinom(x=C+W, size=CN*(nb_r/2), prob=nb_p)]
  cellCNcalls <- probs[, .SD[which.max(CN_ll)], 
                       by=.(sample, cell, chrom, start, end)]
  
  # computing NB probs and calling CN in CNprobs table per segment
  CNprobs[, CN_ll:=dnbinom(x=cov, size=CN*(nb_r/2), prob=nb_p)]
  CNcalls <- CNprobs[, .SD[which.max(CN_ll)], 
                     by=.(sample, chrom, start, end)]
  
  
  # removing extra columns for cleaning the data for plotting phylogenetic tree
  cellCNcalls <- cellCNcalls[start!=0]
  cellCNcalls[, `:=`(nb_p=NULL, class=NULL, C=NULL, W=NULL,expected=NULL, scalar=NULL, nb_r=NULL, CN_ll=NULL)]
  
  # normalize probs (all CN probs for a single-cell should sum up to 1)
  probs[, CN_ll:=CN_ll/sum(CN_ll), by=.(sample, chrom, start, end, cell)]
  
  # sort probs by CN_ll in each cell
  setorder(probs, sample, chrom, start, end, cell, -CN_ll)
  
  # add an extra column for the sum of CN probabilities
  probs[, cumsum_CN_ll:=cumsum(CN_ll),
        by = .(sample, chrom, start, end, cell)]
  
  # shrink the table to have only the parts of the table that fall into the confidence interval
  conf.probs <- conf.probs[, num_conf_CNs:=min(length(which(cumsum_CN_ll < confidence.interval))+1, .N),
                           by = .(sample, chrom, start, end, cell)]
  
  conf.probs <- conf.probs[, head(.SD, num_conf_CNs[1]),
                           by = .(sample, chrom, start, end, cell)]
  
  # define the table that includes the most likely and the confidence interval range for each segment and cell
  conf.CNs <- conf.probs[, cbind(head(.SD,1), min_CN=min(CN), max_CN=max(CN)),
                         by = .(sample, chrom, start, end, cell)]
  
  conf.CNs[, `:=`(CN_ll=NULL, cumsum_CN_ll=NULL, num_conf_CNs=NULL)]
  
  write.table(conf.CNs, 
              file.path(dir, paste0(format(bin.size, scientific = F),"_BFB_cell_confident_CN.table")),
              sep = "\t",
              quote = F)
  
  write.table(cellCNcalls, 
              file.path(dir, paste0(format(bin.size, scientific = F),"_BFB_cell_CNs.table")),
              sep = "\t",
              quote = F)
  
  write.table(probs, 
              file.path(dir, paste0(format(bin.size, scientific = F),"_BFB_cell_CN_probs.table")),
              sep = "\t",
              quote = F)
  
}