#' Plot reads based on directionality
#' @param input.reads A \code{\link{GRanges}} object with Strand-specific read data
#' @param bin_size a integer (e.g. 2000) used to bin data for plotting
#' @param col_a assigns mreads (i.e. Watson) col (e.g "grey", or rgb(1,2,3, max=255))
#' @param col_b assigns preads (i.e. Crick) col (e.g "red", or rgb(1,2,3, max=255))
#' @author Ashley Sanders
#' @export
#' 
#' 
read.plotR <- function(input.reads, bin.size=2000, col_a=rgb(243,165,97, max=255), col_b=rgb(103,139,139, max=255))
{  
  
  # library(GenomicRanges)
  # library(ggplot2)
  
  # bin the data for plotting as histogram
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
  
  
  ## prepare plots data  
  dfplot.reads <- as.data.frame(gr)
  dfplot.reads$preads <- - dfplot.reads$preads	# negative non-ref (preads) counts
  
  
  # plotting function 
  ### p1 <- barplot of plus/minus reads
  p1 <- ggplot(dfplot.reads) + geom_rect(aes(ymin=0, ymax=mreads, xmin=start, xmax=end), fill=col_a, size=0.5)+ # minus (Watson) reads on top
    scale_fill_identity() +
    scale_x_continuous(limits=c( min(start(input.reads)), max(end(input.reads)) )) + #, breaks=NULL, labels(NULL))    
    #scale_y_continuous(limits=c( min(start(input.reads)), max(end(input.reads)) )) + #, breaks=NULL, labels(NULL))    
    geom_rect(aes(ymin=0, ymax=preads, xmin=start, xmax=end), fill=col_b, size=0.5) ### + ylab("depth") # + labs(x=chr) # plus (Crick) reads on btm
  p1 <- p1  #+ theme_light()
  
  
  out<-list(p1)
  
  return(out)
}
