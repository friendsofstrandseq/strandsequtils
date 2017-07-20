#tested
#' generates all binary strings with n 0s and m 1s
#' @param n the number of 0s
#' @param m the number of 1s
#' @author Maryam Ghareghani
#' @export
#' 


allStatus = function(n,m)
{
  allStat = NULL
  status = initialState(n,m)
  
  while(status != FALSE)
  {
    allStat = c(allStat, status)
    status = getNextState(status)
  }
  
  allStat
}