#' Count the number of Watson and Crick reads and return a dataframe of read counts
#' 
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


#' splits a read count data frame based on chromosomes and returns a list containing the splited dataframes
#' 
#' @param counts dataframe containing W and C read counts in some segments
#' @author Maryam Ghareghani
#' @export
#' 

split.chromosomes <- function(counts)
{
  chr = unique(counts$chromosome)
  list.counts = list()
  for (i in 1:length(chr))
  {
    list.counts[[i]] = counts[counts$chromosome == as.character(chr)[i],]
  }
  list.counts
}

# considers only autosomes