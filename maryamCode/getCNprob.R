#' generate CN probability vector
#' @param readCount total number of read counts (W+c) in all cells in this segment
#' @param p the p parameter in the NB distribution
#' @param r dispersion parameter for the total number of read counts (W+C in all cells) in a bin with copy number 2
#' @param segLen the length of the segment
#' @param binLen the length of the bins used for estimating the NB parameters
#' @param maxCN the maximum number of CNs
#' @param alpha the coefficient of NB dispersion parameter for the number of background reads
#' @author Maryam Ghareghani
#' @export
#' 


getCNprob = function(readCount, p, r, segLen, binLen = 100000, maxCN = 5, alpha = 0.05)
{
  CNprob = NULL
  disp = (r/2)*(segLen/binLen) # disp = dispersion par for copy number 1 for this segLen
  
  for (i in 0:maxCN)
  {
    CNprob = c(CNprob, dnbinom(readCount, size = disp*(max(i,alpha)), prob = p))
  }
  
  if (sum(CNprob) != 0)
  {
    CNprob = CNprob/sum(CNprob) # normalization
  }
  else
  {
    expRC = (1-p)*disp/p
    CN = round(readCount/expRC)
    CNprob = c(rep(0,CN),1)
  }
  
  CNprob
}