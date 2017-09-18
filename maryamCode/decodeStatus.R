#tested
#' Converts a binary status to an integer vector status
#' @param binaryStatus a binary string in which the number of 1s between every two consecutive zeros (or before the first zero or after the last zero) indicates a copy number
#' @author Maryam Ghareghani
#' @export
#' 

decodeStatus = function(binaryStatus)
{
  status = NULL
  pos = str_locate_all(binaryStatus, "0")[[1]][,1]
  pos = c(0, pos, nchar(binaryStatus)+1)
  
  for (i in 1:length(pos)-1)
  {
    status = c(status, pos[i+1]-pos[i]-1)
  }
  status
}