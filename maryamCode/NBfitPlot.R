#' create seperate PDF files for different chromosomes and in each file, create plots for W, C, and total read coverage in all cells with the fitted NB distribution
#' @param directory the directory in which we put the plots of NB fitting
#' @param counts list of dataframes containing W and C read counts
#' @param p the estimated p parameter
#' @param dispersion the estimated dispersion parameters
#' @author Maryam Ghareghani
#' @export
#' 


NBfitplots <- function(directory, counts, cellTypes, p, dispersion, bin.size)
{
  #setwd(directory)
  numCells = (ncol(counts[[1]])-3)/2
  for (i in 1:length(counts))
  {
    pdf(paste0(directory,"chr", i, "NBfit.pdf"))
    par(mfrow=c(2,3))
    for (j in 1:numCells)
    {
      Wcount = counts[[i]][,2*j+2]
      Ccount = counts[[i]][,2*j+3]
      
      r = dispersionPar(cellTypes[i,j], dispersion[i,j], binLength = bin.size)
      r = c(r, dispersion[i,j]) # adding the dispersion par for the total number of read counts
      
      allCounts = list(Wcount, Ccount, Wcount + Ccount)
      name = c("Watson", "Crick", "Total")
      
      for (q in 1:3)
      {
        ymax = max(table(allCounts[[q]]))/sum(table(allCounts[[q]]))
        hist(allCounts[[q]], breaks = (-1:max(allCounts[[q]]))+0.5, freq = FALSE, xlab = "", ylab = "", 
             main = paste("cell", j, name[q], ", t =", cellTypes[i,j]), ylim = c(0, ymax))
        x = 0:max(allCounts[[q]])
        v = dnbinom(x, size = r[q], prob = p)
        lines(x,v,col="red")
      }
    }
    dev.off()
  }
}
