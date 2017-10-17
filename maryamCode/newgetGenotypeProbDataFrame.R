#' takes as input a haplotype probTable dataFrame and outout the genotype probTable
#' @param hapProbDF a dataFrame containing haplotype prob table
#' @author Maryam Ghareghani

getGenotypeProbDataFrame = function(hapProbDF)
{
  GTs = list()
  haplotypeNames = colnames(hapProbDF[14:ncol(hapProbDF)])
  GTnames = NULL
  
  for (i in 1:length(haplotypeNames))
  {
    haplotypeNames[i] = substr(haplotypeNames[i], 2, nchar(haplotypeNames[i]))
  }
  
  for (i in 1:length(haplotypeNames))
  {
    sisterHap = sisterHaplotype(haplotypeNames[i])
    sisPos = match(sisterHap, haplotypeNames)
    if (sisPos > i-1)
    {
      GTs[[length(GTs)+1]] = unique(c(sisPos,i))
      GTnames = c(GTnames, haplotypeNames[i])
    }
  }
  
  GTprobDF = hapProbDF[,1:13]
  for (i in 1:length(GTs))
  {
    if (length(GTs[[i]]) == 1)
    {
      GTprobDF = cbind(GTprobDF, hapProbDF[,13+GTs[[i]]])
    } else
    {
      GTprobDF = cbind(GTprobDF, hapProbDF[,13+GTs[[i]][1]]+hapProbDF[,13+GTs[[i]][2]])
    }
  }
  
  colnames(GTprobDF)[14:ncol(GTprobDF)] = GTnames
  
  GTprobDF
}