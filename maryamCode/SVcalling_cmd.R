library(ggplot2)
library(data.table)
library(hexbin)
library(gridExtra)
library(GenomicAlignments)
library(GenomicRanges)
library(MASS)
library(IRanges)
library(stringr)

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
#source('./normalizeProbTable.R')
source('./normalizeProbabilities.R')
source('./getOneCNprobTable.R')
source('./Plot_ReadDistributions_Simple.R')
source('./regularizeProbabilities.R')

args = commandArgs(trailingOnly=TRUE)
directory = args[1]
segCountFile = args[2]

# set parameters
bin.size = 100000
K = 22 # number of chromosomes
maximumCN = 5

# new ...
# reading files... (alternative to the previous lines)
#directory = "/home/mgharegh/research/data/strand-seq/allCells-server/clean/"
p = read.table(paste0(directory, "p.data"))[1,1]
r = read.table(paste0(directory, "r.data"), stringsAsFactors = FALSE)
colnames(r) = r[1,]
r = r[2:nrow(r),]
segmentsCounts = read.table(paste0(directory, "readCounts_", segCountFile), stringsAsFactors = FALSE)
colnames(segmentsCounts) = segmentsCounts[1,]
segmentsCounts = segmentsCounts[2:nrow(segmentsCounts),]
cellTypes = read.table(paste0(directory, "cellTypes.data"), stringsAsFactors = FALSE)
colnames(cellTypes) = colnames(r)
numCells = ncol(cellTypes)

CNhaplotypes = getHapStatesForCN()

hapProbTables = list()
genotypeProbTables = list()

for (i in 1:nrow(segmentsCounts))
{
  print(i)
  hapStates = NULL
  segCounts = segmentsCounts[i,]
  chr = chrNumber(segCounts[1,1])
  
  if (chr > K)
    break()
  
  CN = getPossibleCNs(segCounts, p, as.numeric(r[chr,]), bin.size)
  if (length(CN) > 0 && CN[1] < maximumCN)
  {
    # for (j in 1:length(CN))
    # {
    #   hapStates = c(hapStates, allStatus(3,CN[j]))
    # }
    for (j in 0:maximumCN)
    {
      hapStates = c(hapStates, allStatus(3,j))
    }
    for (j in 1:length(hapStates))
    {
      hapStates[j] = paste(decodeStatus(hapStates[j]), collapse = '')
    }
    hapProbTables[[i]] = getCellStatProbabilities(hapStates, segCounts, as.character(cellTypes[chr,]), p, as.numeric(r[chr,]), binLength = bin.size, alpha = 0.05)
    genotypeProbTables[[i]] = getGenotypeProbTable(hapProbTables[[i]])
  }
}

filterSeg = which(as.numeric(segmentsCounts$end) - as.numeric(segmentsCounts$start) > 100)


# sort cells based on type in each chr
chrOrder = list()
for (i in 1:nrow(cellTypes))
{
  chrOrder[[i]] = c(which(cellTypes[i,] == "wc"), which(cellTypes[i,] == "ww"), which(cellTypes[i,] == "cc"))
}
### new
nonReghapProbTables = hapProbTables
nonRegGenotypeProbTables = genotypeProbTables
for (i in 1:length(hapProbTables))
{
  if (! is.null(hapProbTables[[i]]))
  {
    hapProbTables[[i]] = regularizeProbTable(hapProbTables[[i]])
    genotypeProbTables[[i]] = regularizeProbTable(genotypeProbTables[[i]])
  }
}
### new

for (i in filterSeg)
{
  print(i)
  if (i > length(hapProbTables))
  {
    break()
  }
  
  if (is.null(hapProbTables[[i]])) # CN >= maxCN
  {
    next()
  }
  chr = chrNumber(segmentsCounts[i,1])
  df = data.frame(cell = chrOrder[[chr]], type = as.character(cellTypes[chr, chrOrder[[chr]]]), stringsAsFactors = F)
  
  Wcount = Ccount = NULL
  for (j in chrOrder[[chr]])
  {
    Wcount = c(Wcount, as.integer(segmentsCounts[i,2*j+2]))
    Ccount = c(Ccount, as.integer(segmentsCounts[i,2*j+3]))
  }
  
  df = cbind(df, Wcount)
  df = cbind(df, Ccount)
  
  cellsCNprob = matrix(, nrow = maximumCN+1, ncol = ncol(cellTypes))
  cellsCNprob[1,] = genotypeProbTables[[i]][1,]
  j = 2
  for (l in 1:maximumCN)
  {
    cellsCNprob[l+1,] = colSums(genotypeProbTables[[i]][j:(j+l),])
    j = j+l+1
  }
  
  for (j in 1:(maximumCN+1))
  {
    df = cbind(df, cellsCNprob[j,chrOrder[[chr]]])#rownames(genotypeProbTables[[i]])[j] = )
  }
  
  for (j in 1:nrow(genotypeProbTables[[i]]))
  {
    df = cbind(df, genotypeProbTables[[i]][j,chrOrder[[chr]]])#rownames(genotypeProbTables[[i]])[j] = )
  }
  
  # aggregate cells statistics
  agg = c(ncol(cellTypes) + 1, "all", sum(df[,3]), sum(df[,4]))
  
  for (j in 5:ncol(df))
  {
    agg = c(agg, prod(df[,j]))
  }
  # normalizing
  if (sum(as.numeric(agg[5:10]))>0)
  {
    agg[5:10] = as.numeric(agg[5:10])/sum(as.numeric(agg[5:10]))
  }
  if (sum(as.numeric(agg[11:length(agg)]))>0)
  {
    agg[11:length(agg)] = as.numeric(agg[11:length(agg)])/sum(as.numeric(agg[11:length(agg)]))
  }
  df = rbind(df, agg)
  
  segNames = segmentsCounts[i,1:3][rep(seq_len(nrow(segmentsCounts[i,1:3])), each=nrow(df)),]
  rownames(segNames) = NULL
  df = cbind(segNames, df)
  
  names(df) = c("chr", "start", "end", "cells", "types", "Wcount", "Ccount", paste0(rep("CN",maximumCN+1), as.character(0:maximumCN)), rownames(genotypeProbTables[[i]]))
  
  #dataFrames[[i]] = df
  if (i == 1)
  {
    probTables = df
  } else
  {
    probTables = rbind(probTables, df)
  }
}

write.table(probTables, file = paste0(directory,"allSegCellProbs",segCountFile), sep = "\t", quote = FALSE, row.names = FALSE)
aggProbTable = probTables[probTables$cells == numCells+1,]
SVs = SVcalling(aggProbTable)
SVs = cbind(SVs[,1:3], width = as.numeric(SVs[,3])-as.numeric(SVs[,2]), SVs[,4:ncol(SVs)])
write.table(aggProbTable, file = paste0(directory,"allSegAggregateProbs",segCountFile), sep = "\t", quote = FALSE, row.names = FALSE)

# AshleySeg = read.table("/local/home/mgharegh/research/data/strand-seq/allCells-server/AshleyInversions/segmentsWithStatus.bed", stringsAsFactors = FALSE)
# AshleySeg = AshleySeg[match(SVs[,2], AshleySeg[,2]),]
# avgCount = apply(as.matrix(r), 2, as.numeric)*(1-p)/(p*2)
# avgTotCount = rowSums(avgCount)
# SVs = cbind(SVs[,1:8], data.frame(observedToAverageReadCountFraction = 
#                                     ((as.numeric(SVs$Wcount)+as.numeric(SVs$Ccount))*bin.size)/(avgTotCount[sapply(SVs$chr,chrNumber)]*SVs$width))
#             , SVs[,9:ncol(SVs)], AshleyCall = AshleySeg[,4])
# 
write.table(SVs, file = paste0(directory, "allSegSV", segCountFile), sep = "\t", quote = FALSE, row.names = FALSE)

