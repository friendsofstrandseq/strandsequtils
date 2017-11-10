#' compute the list of all possible haplotypes for each CN between 0 and maxCN
#' 
#' @param maxCN maximum number of copy numbers
#' @author Maryam Ghareghani
#' @export
#' 

getHapStatesForCN = function(maxCN = 5)
{
  hapStates = list()
  for (i in 0:maxCN)
  {
    hapStates[[i+1]] = allStatus(3,i)
  }
  
  hapStates
}