# tested
#' computes the next binary status of a given status with equal number of 0s and 1s
#' @param currentState the current binary status
#' @author Maryam Ghareghani
#' @export
#' 

getNextState = function(currentState)#, n, m)
{
  source('./initializeState.R')
  nextState = currentState
  pos = str_locate_all(currentState, "01")[[1]][,1]
  if (length(pos) == 0)
  {
    return(FALSE) # can't be incremented
  }
  else
  {
    pos = pos[length(pos)] # last occurence of the pattern
    nextState = paste0(substr(nextState,1,pos-1),"10",substr(nextState,pos+2,nchar(nextState)))
    nextState = paste(nextState, collapse = "")
    
    if (pos + 1 < nchar(nextState))
    {
      c1 = str_count(substr(nextState, pos+2, nchar(nextState)),"1")
      c0 = nchar(nextState) - pos - 1 - c1
      nextState = paste0(substr(nextState,1,pos+1), initialState(c0,c1))
    }
  }
  nextState
}