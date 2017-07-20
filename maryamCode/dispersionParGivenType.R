#' Compute the dispersion paramters, given a segment type
#' @param segType a string containing a number of w characters followed by a number of c characters
#' @param r the dispersion parameter for the total (W+C) number of read counts (for CN = 2)
#' @param alpha fraction of the number of reads in minor and major strands
#' @author Maryam Ghareghani
#' @export
#' 

dispersionPar = function(segType, r, segLength, binLength = 100000, alpha = 0.05)
{
  CNw = str_count(segType, "w")
  CNc = str_count(segType, "c")
  
  disp = rep(r/2,2)*(segLength/binLength)*c(CNw, CNc)
  
  for (i in which(disp == 0))
  {
    disp[i] = r*alpha
  }
  disp
}
