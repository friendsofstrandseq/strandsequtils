#' Read bam files from a directory and output them as a \code{\link{GRanges}} object. It also writes the name of the bam files in a file
#' @param directory directory containing the bam files
#' @param bamFilenames name of a file to write the names of the bam files in
#' @param unq Only the unique alignments are given as output, if unq = TRUE
#' @author Maryam Ghareghani
#' @export
#' 

read.bams <- function(directory, bamFilenames = "", unq = TRUE)
{
  setwd(directory)
  files = list.files(pattern = "\\.bam$")
  if (bamFilenames != "")
  {
    write(noquote(files), file = bamFilenames, append = FALSE, sep = "\n")
  }
  
  grlist = list()
  for (i in 1:length(files))
  {
    print(paste("reading", files[i]))
    start.time = Sys.time()
    suppressWarnings(data.raw <- GenomicAlignments::readGAlignmentPairs(file = files[i],
                                                                    param=Rsamtools::ScanBamParam(what='mapq', flag=scanBamFlag(isDuplicate=F))))
                                                                   # , which = GRanges(paste0("chr",chrNum), IRanges(1,chrLen)))))
    print(paste("time for reading bam file:",Sys.time()-start.time))
    data.first <- as(GenomicAlignments::first(data.raw), 'GRanges')
    print(paste("time for subsetting the first mates:",Sys.time()-start.time))
    # quality filter--- getting unique reads
    if (unq)
    {
      data.first <- data.first[mcols(data.first)$mapq > 10]
    }
    print(paste("time for quality filtering:",Sys.time()-start.time))
    grlist[[i]] <- data.first
  }
  grlist
}

# need to be changed... We shouldn't count first pairs only
