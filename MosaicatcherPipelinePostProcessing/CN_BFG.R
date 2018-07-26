bin.size <- 25000
maximumCN <- 200
dir <- "/home/maryam/research/Mosaicatcher/BFB_analysis/C7-7Jun"

# manually define the breakpoints here
brs <- c(33998805, 34185186, 34211372, 35497754, 37452507, 39578189)
# round the breakpoints to be the multiples of bin sizes
breakpoints <- data.table(k = length(unique(round(brs/bin.size))), chrom = "chr10", bps = unique(round(brs/bin.size)))
# write down the breakpoints in a file
brFile <- "/home/maryam/research/Mosaicatcher/BFB_analysis/manual-regions.txt"
write.table(breakpoints, file = brFile, quote = F, sep = "\t", row.names = F)

# get the file names and the directory paths of input files
countsFile <- file.path(dir, paste0(format(bin.size, scientific = F),"_fixed.txt.gz"))
#brFile <- file.path(dir, "segmentation2/C7/manual-regions.txt")
infoFile <- file.path(dir, paste0(format(bin.size, scientific = F),"_fixed.info"))
stateFile <- file.path(dir, "final.txt")

# read the input files
# counts[chrom=="chr10" & start > 33998805 & end < 34185186]
counts <- fread(paste("zcat", countsFile))
segs <- fread(brFile)
info <- fread(infoFile)
strand <- fread(stateFile)


# manually define the segs with read counts in each cell TODO: take the bed file and call getCountsPerRegions
segs <- fread("/home/maryam/research/Mosaicatcher/BFB_analysis/C7_data.RegionCounts.txt")

get_CN_NB_prob <- function(counts, segs, info, strand, manual.segs=TRUE) {

  
  # set the None classes to WW in the counts data table
  counts[class=="None", class:="WW"]
  # prepare the data to be used in SV classification
  probs <- mosaiClassifierPrepare(counts, info, strand, segs, manual.segs = T)
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
  
  write.table(cellCNcalls, 
              file.path(dir, paste0(format(bin.size, scientific = F),"_BFB_cell_CNs.table")),
              sep = "\t",
              quote = F)

}
