# tested
#' computes the genotype status corresponding to the input haplotype status
#' @param haplotyeStatus a decoded haplotype status
#' @author Maryam Ghareghani
#' @export
#' 

getGenotypeStatus = function(haplotypeStatus)
{
  hapStatus = as.integer(strsplit(haplotypeStatus,"")[[1]])
  genotypeStatus = NULL
  for (i in 1:(length(hapStatus)/4))
  {
    genotypeStatus = c(genotypeStatus, hapStatus[4*i-3]+hapStatus[4*i-1])
    genotypeStatus = c(genotypeStatus, hapStatus[4*i-2]+hapStatus[4*i])
  }
  
  paste(genotypeStatus, collapse = "")
}