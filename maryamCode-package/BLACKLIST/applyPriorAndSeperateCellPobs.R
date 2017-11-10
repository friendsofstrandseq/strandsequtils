setwd("/local/home/mgharegh/research/codes/serverCode-Oct17/maryamCode")
source('./WCcount.R')
source('./readBams.R')
source('./splitChr.R')
source('./nonzeroCovIndex.R')#buggy-- the server version is true
source('./uniqeMapIndex.R')#buggy-- the server version is true
source('./filter.R')#buggy-- the server version is true
source('./estimateP.R')
source('./estimateR.R')
source('./NBfitPlot.R')# check and test it later
source('./setCellTypes.R')#buggy-- the server version is true
source('./convertHaplotypeToGenotypeProbTable.R')
source('./convertHaplotypeToGenotypeStatus.R')
source('./decodeStatus.R')
source('./dispersionParGivenType.R')
source('./getSegType.R')
source('./hapStatusCellProbTable.R')
source('./initializeState.R')
source('./nextStatus.R')
source('./singleCellStatatusCalling.R')# not tested-not used
source('./generateAllStates.R')
source('./getCNprob.R')
source('./getPossibleCNs.R')
source('./hapStateForCN.R')
source('./chrNum.R')
source('./SVcalling.R')# add this and the next functions
source('./getCellStatProbTable.R')
source('./normalizeProbabilities.R')
source('./getOneCNprobTable.R')
source('./Plot_ReadDistributions_Simple.R')
source('./regularizeProbabilities.R')
source('./getSisterHaplotype.R')
source('./newgetGenotypeProbDataFrame.R')
source('./newSVcalling.R')

args = commandArgs(trailingOnly=TRUE)
priorFile = args[1]
directory = args[2]
subdir = args[3]

bin.size = 50000
setwd(paste0(directory, "newFormat/", subdir))

getCellNum <- function(Str)
{
  as.numeric(substr(Str, 6, nchar(Str)))+1
}

getAggProbTable = function(GTprobs, prior, numCells)
{
  aggP = data.frame()
  
  for (st in unique(GTprobs$start))
  {
    GTP = GTprobs[GTprobs$start == st,]
    logGTP = log(as.matrix(GTP[,14:ncol(GTP)]))
    agg = c(unlist(GTP[1,1:3]), numCells+1, "all", sum(GTP$Wcount), sum(GTP$Ccount), colSums(logGTP) + log(prior))
    aggP = rbind(aggP, agg, make.row.names = F, stringsAsFactors =F)
  }
  
  colnames(aggP) = colnames(GTprobs)[c(1:7,14:ncol(GTprobs))]
  
  aggP
}


variants = read.table(paste0(directory, "variants.", subdir,".txt"), header = T, stringsAsFactors = F)
prior = read.table(priorFile)[,3]
hapProbs = read.table("allSegCellProbs.table", stringsAsFactors = F, header = T)
GTprobs = getGenotypeProbDataFrame(hapProbs)

#applyPriorAndSeperateCellsPop = function(GTprobs, variants, prior, bin.size) # variants is a dataframe containing SVs and the cell IDs that have the SVs
#{
  GTnames = colnames(GTprobs[,14:ncol(GTprobs)])
  
  SVrows = NULL # rows in the GTprobs that correspond to the cells that have SV
  
  
  variants$start = round(variants$start/bin.size)*bin.size
  starts = intersect(unique(variants$start), unique(GTprobs$start))
  
  for (st in starts)
  {
    GTProwIdx = which(GTprobs$start == st)
    GTP = GTprobs[GTProwIdx,]
    var = variants[variants$start == st,]
    
    SVcellsPop = sapply(var$cell, getCellNum)
    
    subrows = setdiff(match(SVcellsPop,GTP$cells), NA)
    SVrows = c(SVrows, GTProwIdx[subrows])
  }
  
  write.table(GTprobs[SVrows,], file = "SVcellsGTprobs.data", row.names = F, quote = F)
  write.table(GTprobs[setdiff(1:nrow(GTprobs), SVrows),], file = "nonSVcellsGTprobs.data", row.names = F, quote = F)
  
#}




GTprobs = read.table("SVcellsGTprobs.data", stringsAsFactors = F, header = T)
GTprobsNoSVs = read.table("nonSVcellsGTprobs.data", stringsAsFactors = F, header = T)
prior = prior/sum(prior)
variants$start = round(variants$start/bin.size)*bin.size
m = match(unique(variants$start), variants$start)
var = variants[m,]
var = var[which(sapply(var$chrom, chrNumber) < 23),]

agg = getAggProbTable(GTprobs, prior, 100)
aggNoSV = getAggProbTable(GTprobsNoSVs, prior, 100)

SVcellsSV = newSVcalling(agg)
noSVcellsSV = newSVcalling(aggNoSV)

ord = order(var$start)
var = var[ord,]
SVcellsSV = SVcellsSV[ord,]

SVcellsSV = cbind(SVcellsSV, trueSV = var$SV_type)

#SVcodes = c(hom_ref = "X1010", het_del = "X0010", hom_del = "X0000", het_dup = "X1020", hom_dup = "X2020", het_inv = "X0110", hom_inv = "X0101", inv_dup = "X1011")
SVdecodes = c(X1010 = "hom_ref", X0010 = "het_del", X0000 = "hom_del", X1020 = "het_dup", X2020 = "hom_dup", X0110 = "het_inv", X0101 = "hom_inv", X1011 = "inv_dup")

confusionMat = matrix(0L, nrow = 9, ncol = 9)
colnames(confusionMat) = c(SVdecodes, "false_del")
rownames(confusionMat) = c(SVdecodes, "otherGTs")

decodeSV = function(code, SVdecodes)
{
  decode = ""
  if (is.na(match(code, names(SVdecodes))))
  {
    decode = "otherGTs"
  } else {
    decode = SVdecodes[[code]]
  }
  
  decode
}

for (i in 1:nrow(SVcellsSV))
{
  confusionMat[decodeSV(as.character(SVcellsSV$invCNstat)[i], SVdecodes), as.character(SVcellsSV$trueSV)[i]] = confusionMat[decodeSV(as.character(SVcellsSV$invCNstat)[i], SVdecodes), as.character(SVcellsSV$trueSV)[i]]+1
}
t = table(noSVcellsSV$invCNstat)
confusionMat[sapply(names(t), decodeSV, SVdecodes = SVdecodes), "hom_ref"] = confusionMat[sapply(names(t), decodeSV, SVdecodes = SVdecodes), "hom_ref"] + as.vector(t)
write.table(confusionMat, quote = F, file = "confusion.ods", sep = "\t")
