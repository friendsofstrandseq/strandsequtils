#' estimates p parameter 
#' @param counts dataframe containing W and C read counts in bins
#' @author Maryam Ghareghani
#' @export
#' 

estimateP <- function(counts, directory = "")
{
  numCells = (ncol(counts[[1]])-3)/2
  avg = NULL
  s = NULL
  for (i in 1:length(counts))
  {
    for (j in 1:numCells)
    {
      readcount = counts[[i]][,2*j+2]+counts[[i]][,2*j+3] # chr i, cell j
      avg = c(avg, mean(readcount))
      s = c(s, var(readcount))
    }
  }
 
  p = sum(avg*avg)/sum(s*avg)

  if (directory != "")
  { 
    pdf(paste0(directory, "mean-var-plot.pdf"))
    plot(avg, s, xlab = "mean", ylab = "var")

    lines(avg, avg/p, col = "red")
  
    dev.off()
  }

  p
}

