#' set the cell/chromosome types
#' @param alignmemts A list of \code{\link{GRanges}} objects with Strand-specific read data
#' @param K number of chromosomes
#' @author Maryam Ghareghani
#' @export
#' 

getTypes <- function(alignments, K = 22)
{
  types = matrix(, nrow = K, ncol = length(alignments))
  for (j in 1:K)
  {
    for (i in 1:length(alignments))
    {
      aln = alignments[[i]][seqnames(alignments[[i]]) == paste0("chr",j)]
      Nw = length(which(strand(aln) == "-"))
      Nc = length(which(strand(aln) == "+"))
     
      if (Nc == 0)
      {
	types[j,i] = "ww"
	break()
      }
      WCfrac = Nw/Nc
      
      if (WCfrac < 0.1)
      {
        types[j,i] = "cc"
      }
      else if (WCfrac > (2/3) && WCfrac < 1.5)
      {
        types[j,i] = "wc"
      }
      else if (WCfrac > 10)
      {
        types[j,i] = "ww"
      }
      else
      {
        types[j,i] = "?"
      }
    }
  }
  return(types)
}
