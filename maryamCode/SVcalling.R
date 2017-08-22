

SVcalling = function(aggProbTable, maximumCN = 5)
{
  CN = NULL
  SV = NULL
  names = colnames(aggProbTable)
  
  for (i in 1:nrow(aggProbTable))
  {
    mCN = which(as.numeric(aggProbTable[i,8:(maximumCN+8)]) == max(as.numeric(aggProbTable[i,8:(maximumCN+8)])))
    mStat = which(as.numeric(aggProbTable[i,(maximumCN+9):ncol(aggProbTable)]) == max(as.numeric(aggProbTable[i,(maximumCN+9):ncol(aggProbTable)])))
    if (length(mCN) == 1)
    {
      CN = c(CN, names[7+mCN])
    } else
    {
      CN = c(CN, NA)
    }
    
    if (length(mStat) == 1)
    {
      SV = c(SV, names[maximumCN+8+mStat])
    } else
    {
      SV = c(SV, NA)
    }
  }
  
  cbind(aggProbTable[,1:7], data.frame(CNstat = CN, invCNstat = SV))
}