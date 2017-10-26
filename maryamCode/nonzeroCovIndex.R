#' balcklist segments with no coverage
#' @param counts dataframe containing W and C read counts in some segments
#' @author Maryam Ghareghani
#' @export
#' 

nonzero.cov.bins <- function(counts)
{
  idx = list()
  for (i in 1:length(counts))
  {
    idx[[i]] = which(rowSums(counts[[i]][,4:ncol(counts[[i]])]) != 0)
  }
  idx
}

# need to be changed... We shouldn't count first pairs only
