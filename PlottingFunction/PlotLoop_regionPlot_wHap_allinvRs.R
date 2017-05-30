################# ################# ################# ################# ################# 
################# ################# ################# ################# ################# 
#################
### PLOTTING intInv list
### still in progress April 20 2017
library(ggplot2)
library(data.table)
library(hexbin)
source('/Volumes/DISCO/Analysis_Code/Plot_ReadDistributions_segD_feb20.R')
source("/Volumes/DISCO/Analysis_Code/Plot_regionPlot_wHapData_v2.R") # changes inv.p to include 

## prepare segdup track for plots:
sd<-(read.table('/Volumes/DISCO/trio_project/hg38_SegDupsTrack.txt', header=T))
segD<-GRanges(sd[,2:7])
segD$pct <- sd$fracMatch

## read.data to Plot
#for(dataset in c("HG00512", "HG00513", "HG00514", "HG00731", "HG00732", "HG00733", "NA19238", "NA19239", "NA19240")){
# message("processing dataset: ", dataset)

dataset<- "HG00733"
## READS to plot
>  load(paste0("/Volumes/DISCO/trio_project/COMPOSITE_files/", dataset, "_CompositeFile_haplotag_FULL.Rdata"))

## ROIs to plot
> load(paste0('/Volumes/DISCO/inversion_comparisons/Results_April2017_Integration/', dataset,'_integratedInversionList.Rdata')) #intInv
## see Analysis_April2017_finatlIntegration_may1.R for labeling info (lines )
regions<- intInv # inverted region to plot >> May 7

region_gr<-GRanges(seqnames(intInv), IRanges(start=intInv$innerBP_start, end=intInv$innerBP_end)) # innerBP regions >> May 7


plots <- list()
for(i in 1:length(regions)){
  roi<- regions[i]
  
  #roi<- GRanges(as.data.frame(cbind(seqnames="chr3", start=129380821, end=129487939)))
  ##adjust bin to size of variant
  if(width(roi) < 3000) {bin <- 20}else if(width(roi) > 300000){ bin<-20000}else{bin<- 2000}

  #if(length(roi$geno) == 0) { roi$geno<- paste0("roi_", i)}
  lab<-paste(as.character(roi$IDs), as.character(roi$GTs))
  plt<- regionPlot_wHapData(composite.data, roi, ID=lab, bin, segD, regionData=region_gr, col_a="grey70", col_b<- "darkorchid4", col_c <- "grey7")
  
  #pdf(paste0(plotLoc, "genotypingOtherCallsets_PLOTS_example_smallestPB.pdf"), height=8, width=11, onefile = TRUE) # initiate a pdf
  #plt
  #dev.off()
  plots[[1+length(plots)]] <- plt
  
}

#######
message("saving plots")
### generate a single pdf with all ROIs in it
#pdf(paste0("./ROI_plots/", dataset,"_inversion_PLOTS.pdf"), height=8, width=11, onefile = TRUE) # initiate a pdf
plotLoc<-('/Volumes/DISCO/inversion_comparisons/Results_April2017_Integration/')
pdf(paste0(plotLoc, "PLOTs_", dataset,"compositeFile_invIntList.pdf"), height=8, width=11, onefile = TRUE) # initiate a pdf
# '/Volumes/korbel2/StrandSeq/Trio_project/PUR/HG00733/HG00733_CompositeFile_haplotag.Rdata'
bquiet = lapply(plots, print) # take a long time!
dev.off()
