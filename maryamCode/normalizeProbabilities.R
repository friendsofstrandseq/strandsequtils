# tested
#' normalize a probability table to a table in which each column sums up to 1
#' @param probTable a probability table
#' @author Maryam Ghareghani
#' @export
#' 

normalizeProbTable = function(probTable)
{
  colsum = colSums(probTable)
  newProbTable = probTable
  for (j in 1:ncol(probTable))
  {
    if (colsum[j] != 0)
    {
      newProbTable[,j] = newProbTable[,j]/colsum[j]
    }
  }
  
  newProbTable
}