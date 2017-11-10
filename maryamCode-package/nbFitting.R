#' estimates p parameter and saves the empirical mean-var plot in the directory
#' 
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


#' estimates r parameters
#' 
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


#' create seperate PDF files for different chromosomes and in each file, create plots for W, C, and total read coverage in all cells with the fitted NB distribution
#' 
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
