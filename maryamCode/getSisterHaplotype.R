#' takes as input a haplotype state and returns the opposite haplotype state which has the same genotype
#' @param hapState a decoded haplotype state
#' @author Maryam Ghareghani

sisterHaplotype = function(hapState)
{
  sisterHap = ""
  for (i in 1:(nchar(hapState)/4))
  {
    sisterHap = paste0(sisterHap, substr(hapState,i+2,i+3), substr(hapState,i,i+1))
  }
  sisterHap
}