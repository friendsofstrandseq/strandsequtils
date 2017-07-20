#' Count the number of Watson and Crick reads
#' @param cells.alignmemts A list of \code{\link{GRanges}} objects with Strand-specific read data
#' @param segments A list \code{\link{GRanges}} objects (genomic intervals)
#' @author Maryam Ghareghani
#' @export
#' 

WCreadCounts <- function(segments, cells.alignments)
{
  numCells = length(cells.alignments)
  
  df = data.frame(chromosome = seqnames(segments), start = start(segments), end = end(segments))
  
  for (i in 1:numCells)
  {
    reads <- cells.alignments[[i]]
    wdata <- reads[strand(reads) == "-"]
    cdata <- reads[strand(reads) == "+"]
    
    df[[paste0("W",i)]] = GenomicRanges::countOverlaps(segments, wdata)
    df[[paste0("C",i)]] = GenomicRanges::countOverlaps(segments, cdata)
  }
  df
}
