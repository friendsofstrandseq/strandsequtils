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
source('./SVcalling.R')
source('./getCellStatProbTable.R')
source('./normalizeProbabilities.R')
source('./getOneCNprobTable.R')
source('./Plot_ReadDistributions_Simple.R')
source('./regularizeProbabilities.R')

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
  chr = chrNumber(segCounts[,1])
  CN = getPossibleCNs(segCounts, p, as.numeric(r[chr,]))
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
    hapProbTables[[i]] = getCellStatProbabilities(hapStates, segCounts, as.character(cellTypes[chr,]), p, as.numeric(r[chr,]), binLength = 100000, alpha = 0.05)
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

# regularization
nonReghapProbTables = hapProbTables
nonRegGenotypeProbTables = genotypeProbTables
for (i in 1:nrow(segmentsCounts))
{
  if (! is.null(hapProbTables[[i]]))
  {
    hapProbTables[[i]] = regularizeProbTable(hapProbTables[[i]])
    genotypeProbTables[[i]] = regularizeProbTable(genotypeProbTables[[i]])
  }
}

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

subDir = "SVcallingDataRegularizedProb/"
write.table(probTables, file = paste0(temp, subDir,"allSegCellProbs.table"), sep = "\t", quote = FALSE, row.names = FALSE)
aggProbTable = probTables[probTables$cells == 116,]
SVs = SVcalling(aggProbTable)
SVs = cbind(SVs[,1:3], width = as.numeric(SVs[,3])-as.numeric(SVs[,2]), SVs[,4:ncol(SVs)])
write.table(aggProbTable, file = paste0(temp, subDir, "allSegAggregateProbs.table"), sep = "\t", quote = FALSE, row.names = FALSE)

AshleySeg = read.table("/local/data/maryam/data/strand-seq/segmentsWithStatus.bed", stringsAsFactors = FALSE)
AshleySeg = AshleySeg[match(SVs[,2], AshleySeg[,2]),]
avgCount = apply(as.matrix(r), 2, as.numeric)*(1-p)/(p*2)
avgTotCount = rowSums(avgCount)
SVs = cbind(SVs[,1:8], data.frame(observedToAverageReadCountFraction =
                                    ((as.numeric(SVs$Wcount)+as.numeric(SVs$Ccount))*bin.size)/(avgTotCount[sapply(SVs$chr,chrNumber)]*SVs$width))
            , SVs[,9:ncol(SVs)], AshleyCall = AshleySeg[,4])

write.table(SVs, file = paste0(temp, subDir, "allSegSV.table"), sep = "\t", quote = FALSE, row.names = FALSE)

# dirty part
# clustering cells and ordering them based on clusters
#filtSeg = which(as.numeric(segmentsCounts$end) - as.numeric(segmentsCounts$start) > 10000)
gg = cellsStatusProbTable(genotypeProbTables, 1:ncol(cellTypes))# filtSeg)
cellDist = dist(gg, method = "euclidean")
cellsClust = hclust(cellDist, method = "ward.D")
plot(cellsClust)
ord = order.dendrogram(as.dendrogram(cellsClust))

# reading the prob table
gp = read.table(paste0(temp, subDir, "allSegCellProbs.table"), sep = "\t", stringsAsFactors = FALSE)
names = gp[1,]
gp = gp[2:nrow(gp),]
names[13+c(1,2,4,7,11,16)] = c("00","01","02","03","04","05")
colnames(gp) = names
aggP = read.table(paste0(temp, subDir, "allSegAggregateProbs.table"), stringsAsFactors = FALSE)
colnames(aggP) = names
aggP = aggP[2:nrow(aggP),]

# dirty part
# clustering cells and ordering them based on clusters
#filtSeg = which(as.numeric(segmentsCounts$end) - as.numeric(segmentsCounts$start) > 10000)
gg = cellsStatusProbTable(genotypeProbTables, 1:ncol(cellTypes))# filtSeg)
cellDist = dist(gg, method = "euclidean")
cellsClust = hclust(cellDist, method = "ward.D")
plot(cellsClust)
ord = order.dendrogram(as.dendrogram(cellsClust))

# reading the prob table
gp = read.table(paste0(temp, "allSegCellProbs.table"), sep = "\t", stringsAsFactors = FALSE)
names = gp[1,]
gp = gp[2:nrow(gp),]
names[13+c(1,2,4,7,11,16)] = c("00","01","02","03","04","05")
colnames(gp) = names
aggP = read.table(paste0(temp, "allSegAggregateProbs.table"), stringsAsFactors = FALSE)
colnames(aggP) = names
aggP = aggP[2:nrow(aggP),]

# confusion matrix

# output the starting point of the interesting segments
starts = c("19663449","89264184","60881137","62456779")#c("72976532","38966689","50347139","22158427") #c("16226720","50589805","207645514","10284","68334","49336926","10002","131563803","12602621","85741915","103206144","4311483","90024202","76660","86862156","63769526","43323835","12150227")
SVnames = c("inv-dup","het-inv","hom-inv","normal")
for (i in 1:length(starts))#(st in starts)
{
  test = gp[gp$start == starts[i],]
  write.table(test, paste0(temp,"correctlyCalledSegments/",SVnames[i],".table"), quote = F, row.names = F)
}
# test...
st = "72976532"#"62456779"#"179633980"#"49603109"
test = gp[gp$start == st,]
dim(test)
nmat = apply(as.matrix(test[,8:ncol(test)]), 2, as.numeric)
#nmat2 = getBinarySVvector(nmat)
heatmap(nmat, Rowv = as.dendrogram(hclust(dist(nmat, method = "euclidean"), method = "ward.D")), Colv = NA, labCol = colnames(nmat))
test[which(test$CN2 == 0),]
bamNames[as.numeric(test$cells[which(test$`20` == 0)])]
aggRowIdx = which(test$types == "all")
testCN3 = oneCNprobTable(test, 1, aggRowIdx)
nmatCN3 = apply(as.matrix(testCN3[,8:ncol(testCN3)]), 2, as.numeric)
heatmap(nmatCN3, Rowv = as.dendrogram(hclust(dist(nmatCN3, method = "euclidean"), method = "ward.D")), Colv = NA, labCol = colnames(nmatCN3))
testCN4 = oneCNprobTable(test, 2, aggRowIdx)
nmatCN4 = apply(as.matrix(testCN4[,8:ncol(testCN4)]), 2, as.numeric)
heatmap(nmatCN4, Rowv = as.dendrogram(hclust(dist(nmatCN4, method = "euclidean"), method = "ward.D")), Colv = NA, labCol = colnames(nmatCN4))
ord = order.dendrogram(as.dendrogram(hclust(dist(nmat, method = "euclidean"), method = "ward.D")))
test[ord,4:8]
getCNprob((as.numeric(aggP$Wcount)+as.numeric(aggP$Ccount))[which(aggP$start == st)], p, sum(as.numeric(r[7,])), (as.numeric(aggP$end)-as.numeric(aggP$start))[which(aggP$start == st)])
write.table(test, paste0(temp,"SVcalling/test.table"), quote = F, row.names = F)
# global clustering of cells
test2 = t(genotypeProbTables[[which(segmentsCounts$start == st)]])
heatmap(test2, Rowv = as.dendrogram(cellsClust), Colv = NA, labCol = names[14:length(names)])

# bar plot
barplot(as.numeric(aggP[1,][14:ncol(aggP)]), names.arg = names[14:ncol(aggP)])

# some statistics
genotypes = as.character(names[,14:ncol(names)])
confMat = confusionMatrix(SVs, genotypes)
testSeg = intersect(which(SVs$invCNstat == "00"), which(SVs$AshleyCall != "normal"))
SVs[testSeg,]

# WC plots:

