# not tested- not used
#' calls status for each single cell in a segment given a probability table
#' @param probTable a probability table
#' @author Maryam Ghareghani
#' @export
#' 

singleCellStatCall = function(probTable)
{
  status = NULL
  for (i in 1:ncol(probTable))
  {
    maxIdx = which(probTable[,i]==max(probTable[,i]))
    if (length(maxIdx) == 1)
    {
      status = c(status, rownames(probTable)[maxIdx])
    } else {
      status = c(status, NA)
    }
  }
  status
}