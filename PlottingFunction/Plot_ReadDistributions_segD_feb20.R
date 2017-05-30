
#' Plot reads based on directionality
#' @param input.reads A \code{\link{GRanges}} object with Strand-specific read data
#' @param segD.data A \code{\link{GRanges}} object with segmental duplication data for the region
#' @param bin_size a integer (e.g. 2000) used to bin data for plotting
#' @param re.outliers if T will limit the axis to remove outlier
#' @param wlim is the upper wRatio to be called 'ww' - everything b/w wlim & clim is 'wc'
#' @param clim is the lower wRatio to be called 'cc' - everything b/w wlim & clim is 'wc'
#' @param col_a assigns mreads (i.e. Watson) col (e.g "grey", or rgb(1,2,3, max=255))
#' @param col_b assigns preads (i.e. Crick) col (e.g "red", or rgb(1,2,3, max=255))
#' @author Ashley Sanders
#' @export
#' 
#' 
read.plotR <- function(input.reads, segD.data, bin.size, rm.outliers=F, wlim=0.7, clim=0.3, col_a=rgb(243,165,97, max=255), col_b=rgb(103,139,139, max=255))
{  
    
  library(GenomicRanges)
  library(ggplot2)
 

    #bin the data for plotting
    chr<- as.character(seqnames(input.reads)[1])
    max_pos <- max(start(input.reads))
    numbins <- floor(max_pos/bin.size) #calculate number of bins
    ir <- successiveIRanges(rep(bin.size,numbins)) #create continuous ranges based on number of bins and the binsize
    lastBinWidth <- max_pos - (bin.size*numbins)
    lastBinStart <- (max_pos - lastBinWidth) + 1
    lastRange <- IRanges(start = lastBinStart, width = lastBinWidth) #calculate last range
    ir <- c(ir,lastRange)
    #rows <- rep(0,length(ir))
    gr <- GRanges(seqnames=chr, ranges=ir) #initialize GRanges object to store bin read counts
    
    gr$midpoint <- start(gr) + (width(gr)/2)  #get the midposition of each bin
    #counts overlaps between bins and reads in input.reads
    gr$mreads <- GenomicRanges::countOverlaps(gr, input.reads[strand(input.reads)=='-']) 
    gr$preads <- GenomicRanges::countOverlaps(gr, input.reads[strand(input.reads)=='+'])
    gr$readNo <- gr$mreads + gr$preads
    gr$ratio <- gr$mreads/gr$readNo # calculate wcRatio
      gr$strandCall<-'NA' 
      if(length(which(gr$ratio>wlim))> 0){ gr[which(gr$ratio>wlim),]$strandCall <-'ww'}
      if(length(which(gr$ratio<clim))> 0){ gr[which(gr$ratio<clim),]$strandCall <-'cc'}
      if(length(which(gr$ratio>=clim&gr$ratio<=wlim))> 0){gr[which(gr$ratio>=clim & gr$ratio<=wlim),]$strandCall <-'wc'}
   
      ## remove outliers
    if(rm.outliers==T){
    ## identify outliers in pread and mreads   (similar to breakpointR)
      trim=5 # upper/lower 5% of reads will be trimmed to consider sd
      peakTh=0.33 # used to calculate zscores
      zlim=3.291 # outliers have zscore higher than zlim
      maxReads = 500 # will trim any read values higher than this
      
      if (length(which(gr$preads != 0)) > 0){
      trim.p <- gr[gr$preads > 0]$preads # pull out non-0 values
      sds.p <- sd(gr[quantile(trim.p, (trim/100) ) < gr$preads & gr$preads < quantile(trim.p, ((100-trim)/100) )]$preads)
        th.p <- max(gr$preads)*peakTh
        zscores.p <- (gr$preads - th.p) / sds.p
        outlier.idx.p <- which(zscores.p > zlim | gr$preads > maxReads) # find outliers based on zscore
      }else{outlier.idx.p <-0}
      
      if (length(which(gr$mreads != 0)) > 0){
      trim.m <- gr[gr$mreads > 0]$mreads # pull out non-0 values
      sds.m <- sd(gr[quantile(trim.m, (trim/100) ) < gr$mreads & gr$mreads < quantile(trim.m, ((100-trim)/100) )]$mreads)
        th.m<- max(gr$mreads)*peakTh
        zscores.m <- (gr$mreads - th.m) / sds.m
        outlier.idx.m <-which(zscores.m > zlim | gr$mreads > maxReads )# find outliers
      }else{outlier.idx.m <-0}
      
    }else{
      outlier.idx.p <-0
      outlier.idx.m<- 0
    }
      
      # replace outliers with max non-outlier value:
      gr$preads[outlier.idx.p] <- max(gr[-outlier.idx.p]$preads)
      gr$mreads[outlier.idx.m] <- max(gr[-outlier.idx.m]$mreads)
  
      
    ## prepare plots data  
    dfplot.reads <- as.data.frame(gr)
    dfplot.reads$preads <- - dfplot.reads$preads	# negative non-ref (preads) counts
    #select bin which exceeds set th > will be labeled as red in plot
    dfplot.points.minus <- dfplot.reads[outlier.idx.m,] 
    dfplot.points.plus <- dfplot.reads[outlier.idx.p,]
    ##
    
    # plotting function 
    ### p1 <- barplot of plus/minus reads
      p1 <- ggplot(dfplot.reads) + geom_rect(aes(ymin=0, ymax=mreads, xmin=start, xmax=end), fill=col_a, size=0.5)+ # minus (Watson) reads on top
                                  scale_fill_identity() +
                                  scale_x_continuous(limits=c( min(start(input.reads)), max(end(input.reads)) )) + #, breaks=NULL, labels(NULL))    
                                  #scale_y_continuous(limits=c( min(start(input.reads)), max(end(input.reads)) )) + #, breaks=NULL, labels(NULL))    
                                  geom_rect(aes(ymin=0, ymax=preads, xmin=start, xmax=end), fill=col_b, size=0.5) + ylab("depth") # + labs(x=chr) # plus (Crick) reads on btm
            #plot dots at the end of the outlier bins
            if(nrow(dfplot.points.minus) !=0){
              p1 <- p1 + geom_rect(data=dfplot.points.minus, aes(ymin=mreads-(mreads*0.01), ymax=mreads, xmin=start, xmax=end), fill='red')}
            if(nrow(dfplot.points.plus) !=0){
              p1 <- p1 + geom_rect(data=dfplot.points.plus, aes(ymin=preads-(preads*0.01), ymax=preads, xmin=start, xmax=end), fill='red')}
              p1 <- p1 + theme_light()
              
      ### p2 <- dot plot of read ratios
      p2<- ggplot(dfplot.reads) + ylim(c(0,1)) + theme_light()+
        geom_count(aes(x=midpoint, y=ratio, size=readNo, col=strandCall), alpha=0.4) +     
        ylab("ratio")  + #labs(x=chr) +
        #scale_color_brewer(palette="Dark2", direction = -1, drop=FALSE)+
        scale_color_manual(values=c('cc'="#8856a7", 'wc'="#fd8d3c",'ww'="#2ca25f"), drop=F) +
        scale_x_continuous(limits=c(min(start(input.reads)), max(end(input.reads)))) + 
       # theme(axis.text.x = element_blank(), axis.title.x=element_blank())+theme(legend.position="bottom", legend.box="horizontal")+guides(color = guide_legend(order=1, override.aes = list(size=3))) # places strandCall first and increases size of dot
        theme(axis.text.x = element_blank(), axis.title.x=element_blank(), legend.position="none")  # removes legend
        #

      
      ### p3 <- bar plot of segdups
      dfplot.segD <- as.data.frame(segD.data)
      if(nrow(dfplot.segD) >= 1){
        p3 <- ggplot(dfplot.segD) + 
          geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=pct), alpha=0.4) +
          theme_bw() + theme (panel.border = element_blank())+
          ylab("SD")  + labs(x=NULL) + # labs(x=chr) + 
          scale_x_continuous(limits=c(min(start(input.reads)), max(end(input.reads))))+
          theme(axis.text.x = element_blank(), axis.text.y=element_blank()) + #removes axis labels
          #theme(legend.position="bottom", legend.box="horizontal")+
          theme(axis.ticks.y=element_blank() , legend.position="none")+
          scale_fill_gradient(limits=c(0.9, 1), low="skyblue", high="steelblue4")
          #scale_fill_discrete(breaks = c("0.9", "0.95", "1.0"))# doesn't work!
      }else{
        p3 <- ggplot(dfplot.segD) + 
          geom_rect(aes(xmin=start(input.reads), xmax=end(input.reads), ymin=0, ymax=1), alpha=0) +
          theme_bw() + theme (panel.border = element_blank())+
          ylab("SD")  + labs(x=NULL) + # labs(x=chr) + 
          scale_x_continuous(limits=c(min(start(input.reads)), max(end(input.reads))))+
          theme(axis.text.x = element_blank(), axis.text.y=element_blank()) +
          theme(legend.position="none", axis.ticks.y=element_blank())+
          guides(fill = guide_legend(reverse = TRUE))
      }
         


    out<-list(p1, p2, p3)
      
      
    return(out)
}
  
