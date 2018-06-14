## To run the function load required packages below
library(ggplot2)
library(reshape)
library(cowplot)
library(zoo)

## Read in data
tab <- read.table("/media/porubsky/Elements/StrandSeqNation/C7/BFB_cell_CNs.table", header=T, stringsAsFactors = FALSE) 
plotCNheatmap(tab) -> hm.plot



#' Plot clustred heatmap
#'
#' This function plot heatmap aligned with phylogenetic tree from CN data.
#' 
#' @param data.tab A data.frame with following columns: <sample><cell><chrom><start><end><CN>
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @author David Porubsky

plotCNheatmap <- function(data.tab=NULL) {

  ## Cluster data
  CN.mat <- split(data.tab$CN, data.tab$cell) 
  CN.mat <- do.call(rbind, CN.mat)
  euc.dist <- dist(CN.mat, method = "euclidean")
  hc <- hclust(euc.dist, method = "ward.D2")
  
  ## Prepare data for plottting
  CN.mat.df <- as.data.frame(CN.mat)
  colnames(CN.mat.df) <- paste("BFBregion", 1:ncol(CN.mat.df))
  CN.mat.df$ID <- rownames(CN.mat.df)
  rownames(CN.mat.df) <- NULL
  CN.mat.df.ord <- CN.mat.df[hc$order,] #sort data based on clustering
  suppressMessages( plt.data <- melt(CN.mat.df.ord) ) #get long data format
  plt.data$ID <- factor(plt.data$ID, levels=unique(plt.data$ID))
  plt.data$y <- rep(0.5:length(levels(plt.data$ID)), length(unique(plt.data$variable)))
  
  ## Get color categories for CN
  quantile.range <- quantile(plt.data$value, probs = seq(0, 1, 0.2)) #split data into 5 categories
  color.palette <- colorRampPalette(colors = c("white", "black"))(length(quantile.range) - 1)
  label.text <- rollapply(round(quantile.range), width = 2, by = 1, FUN = function(i) paste(i, collapse = " : "))
  color.categs <- findInterval(plt.data$value, quantile.range, all.inside = TRUE)
  plt.data$ColCategs <- color.categs
  
  ## Plot the data
  #Plot Heatmap
  pltlist <- list()
  #plt <- ggplot(plt.data) + geom_tile(aes(x=variable, y=y, fill=value), col='black') + theme_bw() + ylab("") + scale_y_continuous(breaks = plt.data$y, labels = plt.data$ID, expand = c(0,0))
  plt <- ggplot(plt.data) + geom_tile(aes(x=variable, y=y, fill=factor(ColCategs)), col='black') + theme_bw() + ylab("") + scale_y_continuous(breaks = plt.data$y, labels = plt.data$ID, expand = c(0,0))
  plt <- plt + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(face="bold", size=10), axis.line=element_blank(), axis.title.x=element_blank(), plot.margin=margin(l=-0.8,unit="cm"))
  #plt <- plt + scale_fill_gradient2("CN", high = "black", low = "white")
  plt <- plt + scale_fill_manual(values = color.palette, name = "CN", labels = label.text)
  pltlist[['heatmap']] <- plt
  
  #Plot Dendrogram
  dhc <- stats::as.dendrogram(hc)
  ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
  ggdndr <- ggplot(ddata$segments) + geom_segment(aes(x=x, xend=xend, y=y, yend=yend)) + scale_y_continuous(breaks = plt.data$y, labels = plt.data$ID) + scale_x_continuous(expand=c(0,0.5))
  suppressMessages( ggdndr <- ggdndr + scale_y_reverse() )
  ggdndr <- ggdndr + coord_flip()
  ggdndr <- ggdndr + theme(panel.background=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.title=element_blank())
  pltlist[['dendro']] <- ggdndr
  
  ## plot both dendrogram and heatmap together
  complete.plt <- cowplot::plot_grid(plotlist=rev(pltlist), align='h', ncol=length(pltlist), rel_widths = c(2,4))
  return(complete.plt)
}


## old dendrograms
#plot(hc)
#ggtree(as.phylo(hc))
