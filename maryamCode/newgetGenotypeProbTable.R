#' compute the genotype probability table given the haplotype probability table
#' @param haplotypeProbTable a haplotype probability table with the set of decoded (non-binary) haplotypes for rows and cells for columns
#' @author Maryam Ghareghani
#' @export
#' 

newgetGenotypeProbTable = function(haplotypeProbTable)
{
  haplotypes = rownames(haplotypeProbTable)
  genotypes = NULL
  hapIdx = list()
  
  #for (i in 1:length(haplotypes))
  #{
  #  haplotypes[i] = substr(haplotypes[i], 2, nchar(haplotypes[i]))
  #}
  
  for (i in 1:length(haplotypes))
  {
    sisterHap = sisterHaplotype(haplotypes[i])
    sisPos = match(sisterHap, haplotypes)
    if (sisPos > i-1)
    {
      hapIdx[[length(hapIdx)+1]] = unique(c(sisPos,i))
      genotypes = c(genotypes, haplotypes[i])
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
      genotypeProbTable[i,] = colSums(haplotypeProbTable[hapIdx[[i]],])
    }
  }
  
  genotypeProbTable
}