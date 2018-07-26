bin.size <- 100000
maximumCN <- 200
sample.id <- "C7_data"
regions <- "/media/porubsky/Elements/StrandSeqNation/C7/BFB_CN_estimates/BFBregion_C7_CNV_segments.txt"
dir <- "/media/porubsky/Elements/StrandSeqNation/C7/BFB_CN_estimates/"
outputfolder <- dir
bam.folder <- "/media/porubsky/Elements/StrandSeqNation/C7/"

# get the file names and the directory paths of input files (Mosaicatcher outputs)
countsFile <- file.path(dir, paste0(format(bin.size, scientific = F),"_fixed.txt.gz"))
infoFile <- file.path(dir, paste0(format(bin.size, scientific = F),"_fixed.info"))
stateFile <- file.path(dir, "final.txt")

# read the input files
counts <- fread(paste("zcat", countsFile))
segs <- getCountsPerRegions(bam.folder=bam.folder, outputfolder=outputfolder, regions=regions, sample.id=sample.id)
info <- fread(infoFile)
strand <- fread(stateFile)

# run get_CN_NB_prob function
get_CN_NB_prob(counts=counts, segs=segs, info=info, strand=strand, manual.segs=TRUE)

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