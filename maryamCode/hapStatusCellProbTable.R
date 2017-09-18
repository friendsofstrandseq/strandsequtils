#' @param hapStatus a vector of decoded halotype status (strings)
#' @param counts a one row dataframe containing the read counts of the segment
#' @param chrCellTypes cell types in the chromosome where the segment lies
#' @param chrCellsDispPars dispersion parameters of cells in the chromosome where the segment lies
#' @author Maryam Ghareghani
#' @export
#' 


getCellStatProbabilities = function(hapStatus, counts, chrCellTypes, p, chrCellsDispPars, binLength, alpha)
{
  #source('./getSegType.R')
  #source('./dispersionParGivenType.R')
  numCells = length(chrCellTypes)
  pr = matrix(, nrow = length(hapStatus), ncol = numCells)
  rownames(pr) = hapStatus
  colnames(pr) = 1:numCells
  segLen = as.integer(counts[,3]) - as.integer(counts[,2])
  
  for (i in 1:length(hapStatus))
  {
    for (j in 1:numCells)
    {
      segType = getSegType(chrCellTypes[j], hapStatus[i])
      disp = dispersionPar(segType, chrCellsDispPars[j], segLen, binLength, alpha)
      pr[i,j] = dnbinom(as.integer(counts[,2*j+2]), size = disp[1], prob = p)*dnbinom(as.integer(counts[,2*j+3]), size = disp[2], prob = p)
    }
  }
  
  for (j in 1:numCells)
  {
    if (sum(pr[,j]) > 0)
      pr[,j] = pr[,j]/sum(pr[,j])
  }
  
  pr
}