

cellsStatusProbTable = function(genotypeProbTables, segmentsSet)
{
  mat = genotypeProbTables[[segmentsSet[1]]]
  for (i in 2:length(segmentsSet))
  {
    mat = rbind(mat, genotypeProbTables[[segmentsSet[i]]])
  }
  t(mat)
}