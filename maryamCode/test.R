#! /TL/opt/bin/R

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
source('./nonzeroCovIndex.R')
source('./uniqeMapIndex.R')
source('./filter.R')
source('./estimateP.R')
source('./estimateR.R')
source('./NBfitPlot.R')
source('./setCellTypes.R')
source('./dispersionParGivenType.R')
source('./convertHaplotypeToGenotypeProbTable.R')
source('./convertHaplotypeToGenotypeStatus.R')
source('./decodeStatus.R')
source('./dispersionParGivenType.R')
source('./getSegType.R')
source('./initializeState.R')
source('./nextStatus.R')
source('./singleCellStatatusCalling.R')
source('./generateAllStates.R')
source('./getCNprob.R')
source('./getPossibleCNs.R')
source('./hapStateForCN.R')
source('./chrNum.R')
source('./hapStatusCellProbTable.R')

args = commandArgs(trailingOnly=TRUE)
temp = args[1]
bamDir = args[2]

# set parameters
bin.size = 100000
K = 22 # number of chromosomes
maximumCN = 5

# get length of chromosomes
chrLens = read.table(paste0(temp,"chrLens.data"))[,1]

# read bam files
start.time = Sys.time()
cell.alignments = read.bams(bamDir, paste0(temp,"bamFilenames.data"), FALSE)
print(Sys.time() - start.time)
numCells = length(cell.alignments)
#cell.unique.alignments = read.bams("/local/home/mgharegh/research/data/strand-seq/bams/")
cell.unique.alignments = list()
for (i in 1:length(cell.alignments))
{
  cell.unique.alignments[[i]] = cell.alignments[[i]][mcols(cell.alignments[[i]])$mapq > 10]
}
print(Sys.time() - start.time)

# Determining cell types
cellTypes = getTypes(cell.unique.alignments)
write.matrix(cellTypes, file = paste0(temp, "cellTypes.data"))

# make consecutive intervals of length bin.size
bins = data.frame()
for (i in 1:K)
{
  numbins = floor(chrLens[i]/bin.size)
  df = as.data.frame(successiveIRanges(rep(bin.size, numbins)))
  df = cbind(chromosome = rep(paste0("chr",i), numbins), df)
  bins = rbind(bins, df)
}

# constructing a granges object from bins
segments <- GRanges(seqnames=bins$chromosome, ranges=IRanges(start=bins$start, end=bins$end))

# count the number of W and C reads in the bins
counts = WCreadCounts(segments, cell.alignments)
unique.counts = WCreadCounts(segments, cell.unique.alignments)

# output bin read counts to a file
write.table(counts, file = paste0(temp,"binReadCounts.data"), quote = FALSE, row.names = FALSE)
write.table(unique.counts, file = paste0(temp,"binUniqeReadCounts.data"), quote = FALSE, row.names = FALSE)

counts = split.chromosomes(counts)
unique.counts = split.chromosomes(unique.counts)

# blacklisting
nonzeroIndex = nonzero.cov.bins(counts)
counts = filt(counts, nonzeroIndex)
unique.counts = filt(unique.counts, nonzeroIndex)

mappable.bins.index = unique.mappable.bins(counts, unique.counts)
counts = filt(counts, mappable.bins.index)
unique.counts = filt(unique.counts, mappable.bins.index)

# estimate NB parameters
p = estimateP(unique.counts, temp)
disp = estimateR(unique.counts, p)
rownames(disp) = paste0(rep("chr",K), 1:K)
colnames(disp) = paste0(rep("cell", numCells), 1:numCells)

start.time = Sys.time()
NBfitplots(paste0(temp, "NBfitPlots/"), unique.counts, cellTypes, p, disp)
print(Sys.time()-start.time)

# output the NB parameters
write(p, file = paste0(temp, "p.data"))
write.matrix(disp, paste0(temp, "r.data"))

# SV Calling...
# counting W and C reads in segments:
df.segm = read.table(paste0(temp,"segments.bed"), header=F)
colnames(df.segm) = c("chromosome", "start", "end")
segments <- GRanges(seqnames=df.segm$chromosome, ranges=IRanges(start=df.segm$start, end=df.segm$end))
seg.unique.counts = WCreadCounts(segments, cell.unique.alignments)
write.table(seg.unique.counts, file = paste0(temp,"segUniqeReadCounts.data"), quote = FALSE, row.names = FALSE)

# seg.unique.counts = split.chromosomes(seg.unique.counts)

# reading files... (alternative to the previous lines)
p = read.table(paste0(temp, "p.data"))[1,1]
r = read.table(paste0(temp, "r.data"), stringsAsFactors = FALSE)
colnames(r) = r[1,]
r = r[2:nrow(r),]
segmentsCounts = read.table(paste0(temp, "segUniqeReadCounts.data"), stringsAsFactors = FALSE)
colnames(segmentsCounts) = segmentsCounts[1,]
segmentsCounts = segmentsCounts[2:nrow(segmentsCounts),]
cellTypes = read.table(paste0(temp, "cellTypes.data"), stringsAsFactors = FALSE)
colnames(cellTypes) = colnames(r)

CNhaplotypes = getHapStatesForCN()

hapProbTables = list()
genotypeProbTables = list()

for (i in 1:nrow(segmentsCounts))
{
  print(i)
  hapStates = NULL
  segCounts = segmentsCounts[i,]
  chr = chrNumber(segCounts[1,1])
  CN = getPossibleCNs(segCounts, p, as.numeric(r[chr,]))
  if (length(CN) > 0 && CN[1] < maximumCN)
  {
    for (j in 1:length(CN))
    {
      hapStates = c(hapStates, allStatus(3,CN[j]))
    }
    for (j in 1:length(hapStates))
    {
      hapStates[j] = paste(decodeStatus(hapStates[j]), collapse = '')
    }
    hapProbTables[[i]] = getCellStatProbabilities(hapStates, segCounts, as.character(cellTypes[chr,]), p, as.numeric(r[chr,]), binLength = 100000, alpha = 0.05)
    genotypeProbTables[[i]] = getGenotypeProbTable(hapProbTables[[i]])
  }
}

# sort cells based on type in each chr
chrOrder = list()
for (i in 1:nrow(cellTypes))
{
  chrOrder[[i]] = c(which(cellTypes[i,] == "wc"), which(cellTypes[i,] == "ww"), which(cellTypes[i,] == "cc"))
}

dataFrames = list()
for (i in 1:nrow(segmentsCounts))
{
  if (is.null(hapProbTables[[i]])) # CN >= maxCN
  {
    next()
  }
  chr = chrNumber(segmentsCounts[i,1])
  df = data.frame(cell = chrOrder[[chr]], type = as.character(cellTypes[chr, chrOrder[[chr]]]))
  for (j in 1:nrow(genotypeProbTables[[i]]))
  {
    df = cbind(df, genotypeProbTables[[i]][j,chrOrder[[chr]]])#rownames(genotypeProbTables[[i]])[j] = )
  }
  names(df) = c("cells", "types", rownames(genotypeProbTables[[i]]))
  
  dataFrames[[i]] = df
}

# dataFrames = list()
# for (i in 1:nrow(segmentsCounts))
# {
#   if (is.null(hapProbTables[[i]])) # CN >= maxCN
#   {
#     next()
#   }
#   chr = chrNumber(segmentsCounts[i,1])
#   df = data.frame(cell = chrOrder[[chr]], type = as.character(cellTypes[chr, chrOrder[[chr]]]))
#   
#   Wcount = Ccount = NULL
#   for (j in chrOrder[[chr]])
#   {
#     Wcount = c(Wcount, as.integer(segmentsCounts[i,2*j+2]))
#     Ccount = c(Ccount, as.integer(segmentsCounts[i,2*j+3]))
#   }
#   
#   df = cbind(df, Wcount)
#   df = cbind(df, Ccount)
#   
#   cellsCNprob = matrix(, nrow = maximumCN+1, ncol = ncol(cellTypes))
#   cellsCNprob[1,] = genotypeProbTables[[i]][1,]
#   j = 2
#   for (l in 1:maximumCN)
#   {
#     cellsCNprob[l+1,] = colSums(genotypeProbTables[[i]][j:(j+l),])
#     j = j+l+1
#   }
#   
#   for (j in 1:(maximumCN+1))
#   {
#     df = cbind(df, cellsCNprob[j,chrOrder[[chr]]])#rownames(genotypeProbTables[[i]])[j] = )
#   }
#   
#   for (j in 1:nrow(genotypeProbTables[[i]]))
#   {
#     df = cbind(df, genotypeProbTables[[i]][j,chrOrder[[chr]]])#rownames(genotypeProbTables[[i]])[j] = )
#   }
#   
#   names(df) = c("cells", "types", "Wcount", "Ccount", paste0(rep("CN",maximumCN+1), as.character(0:maximumCN)), rownames(genotypeProbTables[[i]]))
#   
#   dataFrames[[i]] = df
# }
