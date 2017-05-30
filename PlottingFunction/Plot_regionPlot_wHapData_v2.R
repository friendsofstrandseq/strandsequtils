
#' Plot reads and phase data for a region in strand-aware fashion
#' @param read.data A \code{\link{GRanges}} object with haplotagged read data.
#' @param roi genomic coordinate that will be plotted, as \code{\link{GRanges}} object (seqnames, start, end)
#' @param ID identifying label for the region (str)
#' @param bin integer (e.g. 2000) used to bin data for plotting
#' @param segD A \code{\link{GRanges}} object with segmental duplication data for the region
#' @param regionData optional \code{\link{GRanges}} containing genomic coordinates, e.g. if more than one ROI is mapped to the area
#' @param col_a assigns mreads (i.e. Watson) col (e.g "grey", or rgb(1,2,3, max=255))
#' @param col_b assigns preads (i.e. Crick) col (e.g "red", or rgb(1,2,3, max=255))
#' @param col_c assigns colour of roi plot
#' @author Ashley Sanders, David Porubsky
#' @export
#' 
regionPlot_wHapData <- function(read.data, roi, ID="roi", bin, segD, regionData=0, col_a=rgb(243,165,97, max=255), col_b=rgb(103,139,139, max=255), col_c="grey7" ){
  library(ggplot2)
  library(data.table)
  library(hexbin)
  source('/Volumes/DISCO/Analysis_Code/Plot_ReadDistributions_segD_feb20.R')
  
  options(warn=-1)
  
  chr<-as.character(seqnames(roi))
  padding<-width(roi)*1.5
  startLoc<- start(roi)-padding
  if (startLoc < 1){startLoc<-1}
  endLoc  <- end(roi) +padding
  file.name<-paste0(chr,":", start(roi), "-",end(roi),".pdf")
  message(paste0("generating region plot for ", file.name))

  input.reads<-read.data[seqnames(read.data) == chr & start(read.data) > startLoc & end(read.data) < endLoc , ]
  segD.data<- segD[which(seqnames(segD)==chr & start(segD) > startLoc & end(segD) < endLoc ), ]
  
  #############################  
  ###  roi PLOT:
  #  roi location plot (simple)
  inv.data<- data.frame(start=start(roi), end=end(roi), genoT=ID) 
  p.inv <- ggplot(inv.data) + geom_rect(aes(ymin=5, ymax=10, xmin=start, xmax=end,fill=col_c), alpha=0.25) +
    theme(legend.position="none", axis.text.x = element_blank())+
    scale_fill_identity() +
    coord_cartesian(xlim=c(min(start(input.reads)), max(end(input.reads))))+
    scale_y_discrete(breaks=NULL, labels(NULL))+ xlab(NULL)+
    geom_text(aes(x=start+((end-start)/2), y=9, label=genoT), col="black")+
    ggtitle(paste0(chr, ":", start(roi), "-", end(roi), " (", prettyNum(width(roi), big.mark=","), " bp)")) + theme_light()

  # add ROI locations to inv plot (if more than one ROI is mapped to the area)
 # if(length(regionData) > 1){
  region_roi <-  as.data.table(regionData[seqnames(regionData) == seqnames(roi) & start(regionData) >= startLoc & end(regionData) <= endLoc,])
  p.inv = p.inv + geom_rect(data = region_roi, inherit.aes = F,
                    aes(xmin = start, ymin=4, xmax = end, ymax = 7), 
                    fill = "black", alpha=0.5)+
                    theme(legend.position="none", axis.text.x = element_blank())
  #}
  #############################  
  ###  read data PLOTS:  
  # increase binsize for thicker bars!
  if(length(input.reads) != 0){
  out<- read.plotR(input.reads, segD.data,  rm.outliers=T, bin.size=bin, col_a=col_a, col_b=col_b, wlim=0.75, clim=0.25)
    p.reads<-out[1][[1]] # barplot of plus/minus reads
    p.ratios<- out[2][[1]] # dot plot of wcRatios
    p.segD<-out[3][[1]] # density plot of segdups
  }else{
    p.reads <- ggplot()+ geom_rect(aes(ymin=0, ymax=10, xmin=startLoc, xmax=endLoc), fill="white", alpha=0)+ geom_text(aes(x=startLoc+((endLoc-startLoc)/2), y=5, label="no reads"), col="grey9")+ scale_fill_identity() + scale_x_continuous(limits=c( startLoc, endLoc )) + ylab("strand.idx")+ theme_light()+xlab(NULL)
    p.ratios<- p.reads +theme(axis.text.x = element_blank(), axis.title.x=element_blank(), legend.position="none") + ylab("ratio") 
    p.segD <- ggplot( as.data.frame(segD.data)) + geom_rect(aes(xmin=startLoc, xmax=endLoc, ymin=0, ymax=1), alpha=0) +  theme_bw() + theme (panel.border = element_blank())+  ylab("SD") + labs(x=NULL) +  scale_x_continuous(limits=c(startLoc,endLoc))+theme(axis.text.x = element_blank(), axis.text.y=element_blank())+ theme(legend.position="none", axis.ticks.y=element_blank())+guides(fill = guide_legend(reverse = TRUE))
  }
  #############################    
  ##### phasing PLOT:
  hapData<-read.data[subjectHits(findOverlaps(GRanges(chr, IRanges(start=startLoc, end=endLoc)), read.data)),] # pull out reads in roi
  hapData<-hapData[hapData$HP!="NA"]
  
  if(length(hapData) != 0){
    hapData$line<-0
    hapData$HP <- as.factor(hapData$HP)
    
    if(table(strand(hapData))["+"] !=0){ hapData[strand(hapData)=="+"]$line <- 1.8} # assign line based on strand (ref on top)
    if(table(strand(hapData))["-"] !=0){ hapData[strand(hapData)=="-"]$line <- -1.8} # assign line based on strand (non-ref on btm)
    hap.plot<- as.data.table(hapData)
    #chr.df <- data.frame(start=min(start(roi)), end=max(end(roi)) )# determine start and end of the selected inverted region to plot the middle line
    
    p<- ggplot(hap.plot) + 
      geom_linerange(aes(x=start, ymin=0, ymax=line, color=HP), alpha=0.3) +
      geom_point(aes(x=start, y=line, color=HP, shape=HP), size=3, alpha=0.4) + scale_shape_manual(values = c(1,2))+
      #geom_jitter(aes(x=start, y=line, color=HP, shape=HP), size=3, alpha=0.4) + scale_shape_manual(values = c(1,2))+
      geom_rect(aes(xmin=startLoc, xmax=endLoc, ymin=-0.3, ymax=0.3), inherit.aes=F)+
      xlab(NULL)+ylab("Ph") +theme_bw()+ylim(2.2,-2.2)+
      theme(legend.position="none", panel.border = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank())+
      coord_cartesian(xlim=c(min(start(input.reads)), max(end(input.reads))))  + ## ONLY if aligning to compData 
      scale_color_brewer(palette="Set1")
  }else{
    p<- ggplot() + 
      geom_rect(aes(xmin=startLoc, xmax=endLoc, ymin=-0.3, ymax=0.3), inherit.aes=F)+
      xlab(NULL)+ylab("Ph") +theme_bw()+ylim(2.2,-2.2)+
      theme(legend.position="none", panel.border = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
      coord_cartesian(xlim=c(min(start(input.reads)), max(end(input.reads))))  + ## ONLY if aligning to compData 
      scale_color_brewer(palette="Set1")
  }
  
  # add ROI location to phase plot
  p = p + geom_rect(data = inv.data, inherit.aes = F,
                    aes(xmin = start, ymin=-0.25, xmax = end, ymax = 0.25), 
                    fill = "grey", alpha=0.5)
  
  plt<- cowplot::plot_grid(p.inv, p.reads, p.ratios, p, p.segD,  ncol=1, align="v", rel_heights = c(0.5, 2, 1.25, 0.75, 0.25)) 
  return(plt)
}



