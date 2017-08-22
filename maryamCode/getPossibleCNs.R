#' returns all CNs with maximum probability
#' @param segCounts W and C counts of one segment (a count dataframe with one row)
#' @param p the p parameter of the NB distribution
#' @param chrCellsDispPars the dispersion parameters of cells for the chromosome where the segment resides
#' @binLength the length of the bins used for NB parameter estimation
#' @alpha the coefficient for the number of background reads
#' @author Maryam Ghareghani
#' @export
#' 


getPossibleCNs = function(segCounts, p, chrCellsDispPars, binLength, alpha)
{
  CN = NULL
  totalR = sum(chrCellsDispPars)
  totalCount = sum(as.integer(segCounts[,4:ncol(segCounts)]))
  segLen = as.integer(segCounts[,3]) - as.integer(segCounts[,2]) + 1
  
  CNp = getCNprob(totalCount, p, totalR, segLen)
  
  if (sum(CNp) > 0) # CN <= maxCN = 5
  {
    CN = which(CNp == max(CNp))-1
  }
  
  CN
}
