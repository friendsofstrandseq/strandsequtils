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
source('./newSVcalling.R')
source('./newhapStatusCellProbTable.R')
source('./getSisterHaplotype.R')

args = commandArgs(trailingOnly=TRUE)
directory = args[1]

# set parameters
bin.size = 50000
K = 22 # number of chromosomes
maximumCN = 5
hapltypeMode = FALSE


binRC = read.table(paste0(directory,"binReadCounts.data"), stringsAsFactors = F)
cellTypes = read.table(paste0(directory,"cellTypes.data"), stringsAsFactors = F)
segmentsCounts = read.table(paste0(directory, "segmentCounts.data"), header = T, stringsAsFactors = F)

numCells = (ncol(binRC)-3)/2

colnames(binRC) = c("chromosome", "start", "end", paste0(rep(c("W","C"),numCells), ceiling(1:(2*numCells)/2)))
colnames(segmentsCounts) = colnames(binRC)

unique.counts = binRC
unique.counts = split.chromosomes(unique.counts)

# blacklisting
nonzeroIndex = nonzero.cov.bins(unique.counts)
unique.counts = filt(unique.counts, nonzeroIndex)

# estimate NB parameters
p = estimateP(unique.counts, directory)
disp = estimateR(unique.counts, p)
rownames(disp) = paste0(rep("chr",K), 1:K)
colnames(disp) = paste0(rep("cell", numCells), 1:numCells)

r = disp

start.time = Sys.time()
print("NB fitting...")
NBfitplots(paste0(directory, "NBfitPlots/"), unique.counts, cellTypes, p, disp, bin.size)
print(Sys.time()-start.time)

# output the NB parameters
write(p, file = paste0(directory, "p.data"))
write.matrix(disp, paste0(directory, "r.data"))

CNhaplotypes = getHapStatesForCN()

hapProbTables = list()
genotypeProbTables = list()

# store the subset of non SCE cells in each chromosome
nonSCEcells = list()
for (i in 1:K)
{
  nonSCEcells[[i]] = which(cellTypes[i,] != "?")
}

hapStates = NULL
for (j in 0:maximumCN)
{
  hapStates = c(hapStates, allStatus(3,j))
}
for (j in 1:length(hapStates))
{
  hapStates[j] = paste(decodeStatus(hapStates[j]), collapse = '')
}

aggProbTable = matrix(, nrow = nrow(segmentsCounts), ncol = length(hapStates))
colnames(aggProbTable) = hapStates

print("computing haplotype probabilities...")
for (i in 1:nrow(segmentsCounts))
{
  print(i)
  segCounts = segmentsCounts[i,]
  chr = chrNumber(segCounts[,1])
  CN = getPossibleCNs(segCounts, p, as.numeric(r[chr,]), bin.size)
  if (length(CN) > 0 && CN[1] < maximumCN)
  {
    ### changed
    hapProbTables[[i]] = newgetCellStatProbabilities(hapStates, segCounts, as.character(cellTypes[chr,]), p, as.numeric(r[chr,]), binLength = bin.size, alpha = 0.05)
    # regularization
    hapProbTables[[i]] = regularizeProbTable(hapProbTables[[i]])
    # computing aggregate Probabilities (based on non SCE cells only)
    for (h in 1:length(hapStates))
    {
      print(h)
      aggProbTable[i,h] = sum(log(hapProbTables[[i]][h,nonSCEcells[[chr]]]))
    }
    #normalizing hapProbTable
    hapProbTables[[i]] = normalizeProbTable(hapProbTables[[i]])
    
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

print("constructing the dataframe...")
aggProbDF = data.frame()

for (i in filterSeg)
{
  print(i)
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
  
  for (h in 1:nrow(hapProbTables[[i]]))
  {
    df = cbind(df, hapProbTables[[i]][h,chrOrder[[chr]]])#rownames(genotypeProbTables[[i]])[j] = )
  }
  
  aggProbVec = as.data.frame(c(segmentsCounts[i,1:3], ncol(cellTypes) + 1, "all", sum(df[,3]), sum(df[,4]), aggProbTable[i,]))
  names(aggProbVec) = c("chr", "start", "end", "cells", "types", "Wcount", "Ccount",hapStates)
  aggProbDF = rbind(aggProbDF, aggProbVec)
  print(paste("nrow =", nrow(aggProbDF)))
  
  segNames = segmentsCounts[i,1:3][rep(seq_len(nrow(segmentsCounts[i,1:3])), each=nrow(df)),]
  rownames(segNames) = NULL
  df = cbind(segNames, df)
  
  names(df) = c("chr", "start", "end", "cells", "types", "Wcount", "Ccount", paste0(rep("CN",maximumCN+1), as.character(0:maximumCN)), rownames(hapProbTables[[i]]))
  
  if (i == 1)
  {
    probTables = df
  } else
  {
    probTables = rbind(probTables, df)
  }
}

print("output tables...")
write.table(probTables, file = paste0(directory,"allSegCellProbs.table"), sep = "\t", quote = FALSE, row.names = FALSE)
SVs = newSVcalling(aggProbDF)
write.table(SVs, file = paste0(directory, "allSegSV.table"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(aggProbDF, file = paste0(directory,"allSegAggProbs.table"), sep = "\t", quote = FALSE, row.names = FALSE)
