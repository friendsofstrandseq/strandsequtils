#' Wrapper function to run SV calling after having the NB parmeters and segment read counts
#'
#' @param bin.size the size of the bins
#' @param K the number of chromosomes (autosomes)
#' @param maximumCN maximum CN in the segments for SV calling
#' @param outputDir the outputDir containing the input and output files
#' @param segCountFile the name of the segments read counts file
#' @author Maryam Ghareghani

# 
# library(ggplot2)
# library(data.table)
# library(hexbin)
# library(gridExtra)
# library(GenomicAlignments)
# library(GenomicRanges)
# library(MASS)
# library(IRanges)
# library(stringr)
# 
# setwd("/local/home/mgharegh/research/codes/serverCode-Oct17/maryamCode")
# source('./WCcount.R')
# source('./readBams.R')
# source('./splitChr.R')
# source('./nonzeroCovIndex.R')#buggy-- the server version is true
# source('./uniqeMapIndex.R')#buggy-- the server version is true
# source('./filter.R')#buggy-- the server version is true
# source('./estimateP.R')
# source('./estimateR.R')
# source('./NBfitPlot.R')# check and test it later
# source('./setCellTypes.R')#buggy-- the server version is true
# source('./convertHaplotypeToGenotypeProbTable.R')
# source('./convertHaplotypeToGenotypeStatus.R')
# source('./decodeStatus.R')
# source('./dispersionParGivenType.R')
# source('./getSegType.R')
# source('./hapStatusCellProbTable.R')
# source('./initializeState.R')
# source('./nextStatus.R')
# source('./singleCellStatatusCalling.R')# not tested-not used
# source('./generateAllStates.R')
# source('./getCNprob.R')
# source('./getPossibleCNs.R')
# source('./hapStateForCN.R')
# source('./chrNum.R')
# source('./SVcalling.R')# add this and the next functions
# source('./getCellStatProbTable.R')
# #source('./normalizeProbTable.R')
# source('./normalizeProbabilities.R')
# source('./getOneCNprobTable.R')
# source('./Plot_ReadDistributions_Simple.R')
# source('./regularizeProbabilities.R')
# source('./newSVcalling.R')
# source('./newhapStatusCellProbTable.R')
# source('./getSisterHaplotype.R')
# source('./newgetGenotypeProbTable.R')
# source('./newgetGenotypeProbDataFrame.R')

#args = commandArgs(trailingOnly=TRUE)
#outputDir = args[1]
#segCountFile = args[2]

#non cmd
# outputDir = "/home/mgharegh/research/data/strand-seq/allCells-server/clean/HGSVC/HG00512/test/"
# segCountFile = "readCounts_5simple_inversions_HGSVC.txt"
# 
# # set parameters
# bin.size = 100000
# K = 22 # number of chromosomes
# maximumCN = 5

SVcalling.wrapper.func = function(bin.size, K, maximumCN, segmentsCounts, r, p, cellTypes, outputDir, hapMode)
{
  # blacklisting chr X and Y
  #segmentsCounts = segmentsCounts[which(sapply(segmentsCounts$chromosome, chrNumber) < 23),]
  
  colnames(cellTypes) = colnames(r)
  numCells = ncol(cellTypes)
  
  #the following line has bugs
  #CNhaplotypes = getHapStatesForCN()
  
  hapProbTables = list()
  normInvCNstatus = list()
  
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
  
  for (i in 1:nrow(segmentsCounts))
  {
    print(i)
    segCounts = segmentsCounts[i,]
    chr = chrNumber(segCounts[1,1])
    CN = getPossibleCNs(segCounts, p, as.numeric(r[chr,]), bin.size)
    if (length(CN) > 0 && CN[1] < maximumCN)
    {
      # computing haplotype probabilities
      hapProbTables[[i]] = newgetCellStatProbabilities(hapStates, segCounts, as.character(cellTypes[chr,]), p, as.numeric(r[chr,]), binLength = bin.size, alpha = 0.05, haplotypeMode = hapMode)
      # regularization
      hapProbTables[[i]] = regularizeProbTable(hapProbTables[[i]])
      # computing aggregate Probabilities (based on non SCE cells only), I think it doesn't matter if we do it after normalization as well
      for (h in 1:length(hapStates))
      {
        aggProbTable[i,h] = sum(log(hapProbTables[[i]][h,nonSCEcells[[chr]]]))
      }
      #normalizing hapProbTable
      hapProbTables[[i]] = normalizeProbTable(hapProbTables[[i]])
      
      #normalInvCNstatus
      normInvCNstatus[[i]] = getGenotypeProbTable(hapProbTables[[i]])
    }
  }
  
  # filter out small segments
  filterSeg = which(as.numeric(segmentsCounts$end) - as.numeric(segmentsCounts$start) > 100)
  
  
  # sort cells based on type in each chr
  chrOrder = list()
  for (i in 1:nrow(cellTypes))
  {
    chrOrder[[i]] = c(which(cellTypes[i,] == "wc"), which(cellTypes[i,] == "ww"), which(cellTypes[i,] == "cc"))
  }
  
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
    cellsCNprob[1,] = normInvCNstatus[[i]][1,]
    j = 2
    for (l in 1:maximumCN)
    {
      cellsCNprob[l+1,] = colSums(normInvCNstatus[[i]][j:(j+l),])
      j = j+l+1
    }
    
    for (j in 1:(maximumCN+1))
    {
      df = cbind(df, cellsCNprob[j,chrOrder[[chr]]])
    }
    
    for (h in 1:nrow(hapProbTables[[i]]))
    {
      df = cbind(df, hapProbTables[[i]][h,chrOrder[[chr]]])
    }
    
    aggProbVec = as.data.frame(c(segmentsCounts[i,1:3], ncol(cellTypes) + 1, "all", sum(df[,3]), sum(df[,4]), aggProbTable[i,]))
    names(aggProbVec) = c("chr", "start", "end", "cells", "types", "Wcount", "Ccount", hapStates)
    aggProbDF = rbind(aggProbDF, aggProbVec)
    #print(paste("nrow =", nrow(aggProbDF)))
    
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
  
  GTprobDF = getGenotypeProbDataFrame(probTables)

  
  write.table(probTables, file = paste0(outputDir,"allSegCellProbs.table"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(GTprobDF, file = paste0(outputDir,"allSegCellGTprobs.table"), sep = "\t", quote = FALSE, row.names = FALSE)
  SVs = newSVcalling(aggProbDF)
  write.table(SVs, file = paste0(outputDir, "allSegSV.table"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(aggProbDF, file = paste0(outputDir,"allSegAggProbs.table"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# calling the function:
#SVcalling.wrapper.func(bin.size = 100000, K = 22, maximumCN = 5, outputDir = "/home/mgharegh/research/data/strand-seq/allCells-server/clean/HGSVC/HG00512/test/", segCountFile = "readCounts_5simple_inversions_HGSVC.txt")
#get segmentCounts as input instead of the segCountsFile