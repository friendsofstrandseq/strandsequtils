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
source('./SVcalling.R')
source('./getCellStatProbTable.R')
source('./normalizeProbabilities.R')
source('./getOneCNprobTable.R')
source('./Plot_ReadDistributions_Simple.R')

args = commandArgs(trailingOnly=TRUE)
temp = args[1]
bamDir = args[2]
outputDir = args[3]

# set parameters
bin.size = 100000
K = 22 # number of chromosomes
maximumCN = 5

# get length of chromosomes
chrLens = read.table(paste0(temp,"chrLens.data"))[,1]

# read bam files
start.time = Sys.time()
cell.alignments = read.bams(bamDir, paste0(outputDir,"bamFilenames.data"), FALSE)
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
write.matrix(cellTypes, file = paste0(outputDir, "cellTypes.data"))

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
write.table(counts, file = paste0(outputDir,"binReadCounts.data"), quote = FALSE, row.names = FALSE)
write.table(unique.counts, file = paste0(outputDir,"binUniqeReadCounts.data"), quote = FALSE, row.names = FALSE)

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
p = estimateP(unique.counts, outputDir)
disp = estimateR(unique.counts, p)
rownames(disp) = paste0(rep("chr",K), 1:K)
colnames(disp) = paste0(rep("cell", numCells), 1:numCells)

start.time = Sys.time()
NBfitplots(paste0(outputDir, "NBfitPlots/"), unique.counts, cellTypes, p, disp)
print(Sys.time()-start.time)

# output the NB parameters
write(p, file = paste0(outputDir, "p.data"))
write.matrix(disp, paste0(outputDir, "r.data"))

# counting W and C reads in segments:
for (s in c("simple_inversions_HGSVC.txt", "complex_inversions_HGSVC.txt"))
{
  df.segm = read.table(paste0(temp,s), header=F)
  colnames(df.segm) = c("chromosome", "start", "end")
  segments <- GRanges(seqnames=df.segm$chromosome, ranges=IRanges(start=df.segm$start, end=df.segm$end))
  seg.unique.counts = WCreadCounts(segments, cell.unique.alignments)
  write.table(seg.unique.counts, file = paste0(outputDir,"readCounts_",s), quote = FALSE, row.names = FALSE)
}
