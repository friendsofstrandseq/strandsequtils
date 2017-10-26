#' estimates r parameters
#' @param counts dataframe containing W and C read counts in bins
#' @param p estimated p parameter
#' @author Maryam Ghareghani
#' @export
#' 

estimateR <- function(counts, p)
{
  numCells = (ncol(counts[[1]])-3)/2
  avg = NULL
  s = NULL
  dispersion = matrix(, nrow = length(counts), ncol = numCells)
  for (i in 1:length(counts))
  {
    for (j in 1:numCells)
    {
      avg = mean(counts[[i]][,2*j+2]+counts[[i]][,2*j+3])
      dispersion[i,j] = avg*p/(1-p)
    }
  }
  
  dispersion
}