## To run the function below
library(GenomicAlignments)
sample.id <- "C7_data"
bam.folder <- "/media/porubsky/Elements/StrandSeqNation/C7/"
outputfolder <- "/media/porubsky/Elements/StrandSeqNation/C7/BFB_CN_estimates/"
regions <- "/media/porubsky/Elements/StrandSeqNation/C7/BFB_CN_estimates/BFBregion_C7_CNV_segments.txt"

getCountsPerRegions(bam.folder=bam.folder, outputfolder=outputfolder, regions=regions, sample.id=sample.id) -> count.table

####################################################################################################################

#' This funcion will go through all BAM files present in a 'bam.folder' and will
#' counts directional reads in specified 'regions'.
#' 
#' NOTE: Works only for paired-end reads!!!
#' 
#' @param bam.folder A folder where are stored BAMs.
#' @param outputfolder A folderr where the final counts will be stored.
#' @param regions A bed file with regions to count reads in.
#' @param sample.id A sample identifier.
#' 
#' @author David Porubsky

getCountsPerRegions <- function(bam.folder=NULL, outputfolder=NULL, regions=NULL, sample.id=NULL) {
 
  ## Load required packages
  if ("GenomicAlignments" %in% rownames(installed.packages())) {
    require(GenomicAlignments)
  } else {
    stop("Package 'GenomicAlignments' cannot be loaded!!!")
  } 
  
  BFB.regions <- read.table(regions, header=FALSE, stringsAsFactors = FALSE)
  BFB.regions.gr <- GenomicRanges::GRanges(seqnames=BFB.regions$V1, ranges=IRanges(start=BFB.regions$V2, end=BFB.regions$V3))
  
  bams <- list.files(bam.folder, pattern = "\\.bam$", full.names = TRUE)
  
  file.header <- Rsamtools::scanBamHeader(bams[1])[[1]]
  chrom.lengths <- file.header$targets
  seqlengths(BFB.regions.gr) <- chrom.lengths[seqlevels(BFB.regions.gr)]
  
  #allReads <- GRangesList()
  binned.count.grl <- GRangesList()
  for (bam in bams) {
    filename <- basename(bam)
    cell <- unlist(strsplit(filename, "\\."))[1]
    message("Processing ", filename, " ...")
    ## Load data (paired-end only)
    suppressWarnings( data.raw <- GenomicAlignments::readGAlignmentPairs(file = bam, param=Rsamtools::ScanBamParam(which=range(BFB.regions.gr), what='mapq', flag=Rsamtools::scanBamFlag(isDuplicate=FALSE))) )
    ## Use only first mate of the pair
    data <- as(GenomicAlignments::first(data.raw), "GRanges")
    ## filter by mapping quality 10
    data <- data[mcols(data)$mapq >= 10]
    ## Count reads and creat bins object
    bins <- BFB.regions.gr
    bins$sample <- sample.id
    bins$cell <- cell
    bins$C <- countOverlaps(bins, data[strand(data) == "+"])
    bins$W <- countOverlaps(bins, data[strand(data) == "-"])
    bins$class <- "None"
    binned.count.grl[[bam]] <- bins
  } 
  
  binned.count.gr <- unlist(binned.count.grl, use.names = FALSE)
  binned.count.gr <- sort(binned.count.gr)
  count.table.df <- as.data.frame(binned.count.gr)
  colnames(count.table.df)[1] <- 'chrom'
  remove.cols <- which(colnames(count.table.df) %in% c('width','strand'))
  count.table.df <- count.table.df[,-remove.cols]
  
  destination <- file.path(outputfolder, paste0(sample.id, ".RegionCounts.txt"))
  write.table(count.table.df, file=destination, quote = FALSE, row.names = FALSE, sep = "\t")
  
  return(count.table.df)
}