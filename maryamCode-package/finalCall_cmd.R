#' Counts the number of Watson and Crick reads of single cells in bins and segments and fitting NB distribution
#'
#' @param bin.size the size of the bins
#' @param K the number of chromosomes (autosomes)
#' @param temDir the directory containing the chromosomes length file
#' @param bamDir the directory containing all of the bam files
#' @param directory the directory containing the input and output files
#' @param segmentsFile the name of the segments bed file
#' @author Maryam Ghareghani

library(ggplot2)
library(data.table)
library(hexbin)
library(gridExtra)
library(GenomicAlignments)
library(GenomicRanges)
library(MASS)
library(IRanges)
library(stringr)


args = commandArgs(trailingOnly=TRUE)
print(args)
#args = strsplit(args, " ")[[1]]
#print(args)
print(class(args))
args = as.data.frame(strsplit(args, split = "="), stringsAsFactors = F)
print(args)

Rdirectory = args[2,match("Rdirectory", as.character(args[1,]))]
print(Rdirectory)

prevDir = getwd()
setwd(Rdirectory)

source('chrNum.R')
source('filter.R')
source('generateHaploStates.R')
source('getDispParAndSegType.R')
source('getPossibleCNs.R')
source('hapAndGTprobsComputing.R')
source('nbFitting.R')
source('Plot_ReadDistributions_Simple.R')
source('readBams.R')
source('SVcalling.R')
source('WCcount.R')
source('wrapper_formatChange.R')
source('wrapper_readBamsWCcountAndNBfit.R')
source('wrapper_SVcalling.R')

setwd(prevDir)

binRCfile = args[2,match("binRCfile", as.character(args[1,]))]
BRfile = args[2,match("BRfile", as.character(args[1,]))]
infoFile = args[2,match("infoFile", as.character(args[1,]))]
stateFile = args[2,match("stateFile", as.character(args[1,]))]
outputDir = args[2,match("outputDir", as.character(args[1,]))]
bin.size = as.numeric(args[2,match("bin.size", as.character(args[1,]))])
K = as.numeric(args[2,match("K", as.character(args[1,]))])
maximumCN = as.numeric(args[2,match("maximumCN", as.character(args[1,]))])
haplotypInfo=F
if (any(as.character(args[1,])=="haplotypeInfo")){haplotypInfo = T}

print(paste("Rdirectory =", Rdirectory))
print(paste("binRCfile =", binRCfile))
print(paste("BRfile =", BRfile))
print(paste("infoFile =", infoFile))
print(paste("stateFile =", stateFile))
print(paste("outputDir =", outputDir))
print(paste("bin.size =", bin.size))
print(paste("K =", K))
print(paste("maximumCN =", maximumCN))

binRC = split.chromosomes(changeRCformat(binRCfile, outputDir))
cellTypes = changeCellTypesFormat(stateFile)
NBparams = changeNBparamsFormat(infoFile, K)
p = NBparams[[1]]
r = NBparams[[2]]
segmentsCounts = getSegReadCounts(binRC, BRfile, K, bin.size)

SVcalling.wrapper.func(bin.size, K, maximumCN, segmentsCounts, r, p, cellTypes, outputDir, hapMode = haplotypInfo)

# Rscript finalCall_cmd.R Rdirectory=/local/home/mgharegh/research/codes/serverCode-Oct17/maryamCode/ binRCfile=/local/home/mgharegh/research/HDhackathon/simulatedData/counts.cov5.vaf0.5.small.p0.3.txt.gz BRfile=/local/home/mgharegh/research/HDhackathon/simulatedData/breakpoints.cov5.vaf0.5.small.p0.3.txt infoFile=/local/home/mgharegh/research/HDhackathon/simulatedData/D2Rfb.50kb_fixed.info stateFile=/local/home/mgharegh/research/HDhackathon/simulatedData/sces.cov5.vaf0.5.small.p0.3.txt outputDir=/local/home/mgharegh/research/HDhackathon/simulatedData/ bin.size=50000 K=22 maxCN=5
# binRCfile="/local/home/mgharegh/research/HDhackathon/simulatedData/counts.cov5.vaf0.5.small.p0.3.txt.gz"
# stateFile="/local/home/mgharegh/research/HDhackathon/simulatedData/sces.cov5.vaf0.5.small.p0.3.txt" # "D2Rfb.strand_states.txt"
# outputDir="/local/home/mgharegh/research/HDhackathon/simulatedData/"
# bamNamesFile="bamNames.txt"
# infoFile="/local/home/mgharegh/research/HDhackathon/simulatedData/D2Rfb.50kb_fixed.info"
# breakpointsFile="/local/home/mgharegh/research/HDhackathon/simulatedData/breakpoints.cov5.vaf0.5.small.p0.3.txt"
# Rdirectory="/local/home/mgharegh/research/codes/serverCode-Oct17/maryamCode"

Rdirectory = "/home/mgharegh/research/codes/strandsequtils-shrinkingAnadMergingFuncs/maryamCode/"
binRCfile = "/local/home/mgharegh/research/HDhackathon/data/skin/D2Rfb.100000_fixed.txt.gz"
BRfile = "/local/home/mgharegh/research/HDhackathon/data/skin/D2Rfb.100000_fixed.many.txt"
infoFile = "/local/home/mgharegh/research/HDhackathon/data/skin/D2Rfb.100000_fixed.info"
stateFile = "/local/home/mgharegh/research/HDhackathon/data/skin/D2Rfb.final.txt"
outputDir = "/local/home/mgharegh/research/HDhackathon/data/skin/"
K = 22
maximumCN = 5
bin.size = 100000
