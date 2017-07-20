# tested
#' compute a genotype probability table given the haplotype probability table
#' @param haplotypeProbTable a haplotype probability table with the set of decoded (non-binary) haplotypes for rows and cells for columns
#' @author Maryam Ghareghani
#' @export
#' 

getGenotypeProbTable = function(haplotypeProbTable)
{
  haplotypes = rownames(haplotypeProbTable)
  genotypes = NULL
  hapIdx = list()
  
  for (i in 1:length(haplotypes))
  {
    #genotypeStatus = getGenotypeStatus(decodeStatus(haplotypes[i]))
    genotypeStatus = getGenotypeStatus(haplotypes[i])
    idx = which(genotypes == genotypeStatus)
    if (length(idx) == 0)
    {
      genotypes = c(genotypes, genotypeStatus)
      hapIdx[[length(hapIdx)+1]] = i
    }
    else
    { # length(idx) should be 1
      hapIdx[[idx]] = c(hapIdx[[idx]], i)
    }
  }
  
  genotypeProbTable = matrix(, nrow = length(genotypes), ncol = ncol(haplotypeProbTable))
  rownames(genotypeProbTable) = genotypes
  
  for (i in 1:length(genotypes))
  {
    if (length(hapIdx[[i]]) == 1)
    {
      genotypeProbTable[i,] = haplotypeProbTable[hapIdx[[i]],]
    } else {
      genotypeProbTable[i,] = colSums(haplotypeProbTable[hapIdx[[i]],])#/length(hapIdx[[i]]) # I think we shouldn't devide it by the number of haplotypes corresponding to the genotype
    }
  }
  
  genotypeProbTable
}