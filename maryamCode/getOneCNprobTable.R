

oneCNprobTable = function(genotypeProbTables, CN, AggRowIdx, colNum = 14)
{
  startCol = colNum + (CN*(CN+1)/2)
  
  newTable = apply(as.matrix(genotypeProbTables[,startCol:(startCol+CN)]), 1, as.numeric)
  for (i in 1:nrow(newTable))
  {
    newTable[i,aggRowIdx] = prod(newTable[i,setdiff(1:ncol(newTable), AggRowIdx)])
  }
  newTable = cbind(genotypeProbTables[,1:7], t(normalizeProbTable(newTable)))
  colnames(newTable) = c(colnames(genotypeProbTables)[1:7], colnames(genotypeProbTables)[startCol:(startCol+CN)])
  newTable
}