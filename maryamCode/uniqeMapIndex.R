#' exclude segments with no coverage
#' @param counts dataframe containing W and C read counts in segments
#' @param uniqe.counts dataframe containing W and C unique read counts in segments
#' @author Maryam Ghareghani
#' @export
#' 

unique.mappable.bins <- function(counts, unique.counts, alpha = 0.05)
{
  idx = list()
  for (i in 1:length(counts))
  {
    idx[[i]] = which(rowSums(unique.counts[[i]][,4:ncol(unique.counts[[i]])])/rowSums(counts[[i]][,4:ncol(counts[[i]])]) > (1-alpha))
  }
  idx
}

# need to be changed... We shouldn't count first pairs only
