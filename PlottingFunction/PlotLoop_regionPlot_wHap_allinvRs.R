################# ################# ################# ################# ################# 
################# ################# ################# ################# ################# 
#################
### PLOTTING intInv list
### still in progress April 20 2017
library(ggplot2)
library(data.table)
library(hexbin)
source('./Plot_ReadDistributions_segD_feb20.R')
source("./Plot_regionPlot_wHapData_v2.R") # changes inv.p to include 

## prepare segdup track for plots:
sd<-(read.table('./hg38_SegDupsTrack.txt', header=T))
segD<-GRanges(sd[,2:7])
segD$pct <- sd$fracMatch

## prepare input data:
dataset<- "HG00733"
## READS to plot > GRanges object with read data
load(paste0("./", dataset, "_StrandseqData_Example.Rdata")) # example.data
#
## ROIs to plot > GRanges list of regions
load(paste0('./', dataset,'_integratedInversionList.Rdata')) # intInv 
regions<- intInv # inverted region to plot >> May 7

# region_gr<-GRanges(seqnames(intInv), IRanges(start=intInv$innerBP_start, end=intInv$innerBP_end)) # innerBP regions >> May 7

## plotting loop:
plots <- list()
for(i in 1:length(regions)){
  roi<- regions[i]

  ##adjust bin to size of variant
  if(width(roi) < 3000) {bin <- 20}else if(width(roi) > 300000){ bin<-20000}else{bin<- 2000}

  #if(length(roi$geno) == 0) { roi$geno<- paste0("roi_", i)} # in case roi file isn't genotyped
  #lab<-paste(as.character(roi$IDs), as.character(roi$GTs))
  #plt<- regionPlot_wHapData(composite.data, roi, ID=lab, bin, segD, regionData=region_gr, col_a="grey70", col_b<- "darkorchid4", col_c <- "grey7")
  plt<- regionPlot_wHapData(example.data, roi, ID=paste0("roi_", i), bin, segD, regionData=regions, col_a="grey70", col_b<- "darkorchid4", col_c <- "grey7")
  
  plots[[1+length(plots)]] <- plt
  
}

#######
message("saving plots")
### generate a single pdf with all ROIs in it
pdf(paste0("./PLOTs_", dataset,"compositeFile_invIntList.pdf"), height=8, width=11, onefile = TRUE) # initiate a pdf
bquiet = lapply(plots, print) # takes a long time!
dev.off()


