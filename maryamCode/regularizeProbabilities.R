#' We assume that the posteriori probability for different status is uniform with probability regFactor and a probability computed based on the NB model with probability (1-regFactor)
#' regularize the probability table according to the mentioned assumption
#' @param probTable A non regularized prob table in which every row corresponds to a status and each column corresponds to a single cell
#' @param regFactor the regularization factor

regularizeProbTable = function(probTable, regFactor = 0.001)
{
  newProbTable = probTable
  for (j in 1:ncol(probTable))
  {
    newProbTable[,j] = regFactor/nrow(probTable) + (1-regFactor)*newProbTable[,j]
    newProbTable[,j] = newProbTable[,j]/sum(newProbTable[,j])
  }
  newProbTable
}