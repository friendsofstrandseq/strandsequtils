## To run the function load required packages below
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(cowplot)
library(zoo)

## Read in data
tab1 <- read.table("/media/porubsky/Elements/StrandSeqNation/C7/BFB_cell_RPKM_CNs.table", header=T, stringsAsFactors = FALSE) 
tab2 <- read.table("/media/porubsky/Elements/StrandSeqNation/C7/BFB_CN_estimates/100000_BFB_cell_CNs.table", header=T, stringsAsFactors = FALSE) 

plotCNheatmap(tab1) -> hm.plot1 #RPKM
plotCNheatmap(tab2) -> hm.plot2 #CN estim

#' Plot clustred heatmap
#'
#' This function plot heatmap aligned with phylogenetic tree from CN data.
#' 
#' @param data.tab A data.frame with following columns: <sample><cell><chrom><start><end><CN>
#' @param continuous If to use contunous color scheme
#' @param max.CN Set outliers to max.CN
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @author David Porubsky

plotCNheatmap <- function(data.tab=NULL, continuous=TRUE, max.CN=0) {

  ## Set CN more than 'max.CN' to max.CN
  if (max.CN > 0) {
    mask <- data.tab$CN >= max.CN
    data.tab$CN[mask] <- max.CN
  }  
  
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
  if (continuous) {
    plt <- ggplot(plt.data) + geom_tile(aes(x=variable, y=y, fill=value)) + theme_bw() + ylab("") + scale_y_continuous(breaks = plt.data$y, labels = plt.data$ID, expand = c(0,0))
    plt <- plt + scale_fill_gradientn(colors = brewer.pal(name = 'BuPu', n = 9), name="CN")
  } else {
    plt <- ggplot(plt.data) + geom_tile(aes(x=variable, y=y, fill=factor(ColCategs)), col='black') + theme_bw() + ylab("") + scale_y_continuous(breaks = plt.data$y, labels = plt.data$ID, expand = c(0,0))
    plt <- plt + scale_fill_manual(values = color.palette, name = "CN", labels = label.text)
  }
  plt <- plt + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(face="bold", size=10), axis.line=element_blank(), axis.title.x=element_blank(), plot.margin=margin(l=-0.8,unit="cm"))
  #plt <- plt + scale_fill_manual(values = color.palette, name = "CN", labels = label.text)
  #plt <- plt + scale_fill_gradient2("CN", high = "black", low = "white")
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

############################################################################################################
# TEST code

tab <- read.table("/media/porubsky/Elements/StrandSeqNation/C7/BFB_CN_estimates/100000_BFB_cell_CNs.table", header=T, stringsAsFactors = FALSE) 
data.tab <- tab

#CN.mat <- split(data.tab$CN, data.tab$cell) 
#CN.mat <- do.call(rbind, CN.mat)
#euc.dist <- dist(CN.mat, method = "euclidean")
#hc <- hclust(euc.dist, method = "ward.D2")
#cell.ord <- hc$order

data.tab.bfb <- split(data.tab, data.tab$start)

plots <- list()
cell.orders <- list() 
for (i in 1:length(data.tab.bfb)) {
  plt.df <- data.tab.bfb[[i]]
  cell.orders[[i]] <- order(plt.df$CN)
} 

plots <- list()
for (i in (length(cell.orders)-1):1) {
  ord <- cell.orders[[i+1]]
  CN.ord <-data.tab.bfb[[i+1]]$CN[ord]
  
  prev.seg <- data.tab.bfb[[i]]
  prev.seg <- prev.seg[ord,]
  prev.seg$group <- CN.ord
  plt <- ggplot(prev.seg[order(prev.seg$CN),]) + geom_point(aes(x=1:nrow(prev.seg), y=CN, color=factor(group))) + scale_color_discrete(guide = FALSE) + xlab("")
  plots[[i]] <- plt
  #prev.seg <- split(data.tab.bfb[[i-1]], CN.ord)
  #prev.seg.srt <- lapply(prev.seg, function(x) x[order(x[7]),])
  #prev.seg.srt <- do.call(rbind, prev.seg.srt)
}
ord <- cell.orders[[length(plots)+1]]
CN.ord <-data.tab.bfb[[length(plots)+1]]$CN[ord]
last.seg <- data.tab.bfb[[length(plots)+1]]
last.seg <- last.seg[ord,]
last.seg$group <- CN.ord
plots[[length(plots)+1]] <- ggplot(last.seg[order(last.seg$CN),]) + geom_point(aes(x=1:nrow(prev.seg), y=CN, color=factor(group))) + scale_color_discrete(guide = FALSE) + xlab("")

plot_grid(plotlist = plots, ncol = 1)


ord <- cell.orders[[6]]
plots <- list()
for (i in length(cell.orders):1) {
  CN.ord <-data.tab.bfb[[i]]$CN[ord]
  
  seg <- data.tab.bfb[[i]]
  seg <- seg[ord,]
  seg$group <- CN.ord
  plt <- ggplot(seg) + geom_point(aes(x=1:nrow(prev.seg), y=CN, color=factor(group))) + scale_color_discrete(guide = FALSE) + xlab("")
  plots[[i]] <- plt
}
plot_grid(plotlist = plots, ncol = 1)


ord <- cell.orders[[6]]
plt.df <- data.tab.bfb
for (i in 1:length(plt.df)) {
  df <- plt.df[[i]]
  #df$x <- 1:nrow(df)
  df <- df[ord,]
  df$x <- 1:nrow(df)
  plt.df[[i]] <- df
}
plt.df <- do.call(rbind, plt.df)
ggplot(plt.df) + geom_point(aes(x=x, y=CN, color=factor(start))) + coord_trans(y='log') + xlab("Cells") + scale_color_manual(values = brewer.pal(n = 6, name = "Set1"))



  plt.df <- plt.df[order(plt.df$CN),]
  #plt.df$x <- cell.ord 
  plt <- ggplot(plt.df) + geom_point(aes(x=1:nrow(plt.df), y=CN))# + coord_trans(y='log')
  #plt <- ggplot(plt.df) + geom_point(aes(x=x, y=CN)) + coord_trans(y='log')
  plots[[i]] <- plt
}  
