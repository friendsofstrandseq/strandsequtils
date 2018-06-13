library(ggplot2)
library(reshape)
library(cowplot)

## Read in data
tab <- read.table("/media/porubsky/Elements/StrandSeqNation/C7/BFB_cell_CNs.table", header=T, stringsAsFactors = FALSE) 

## Cluster data
CN.mat <- split(tab$CN, tab$cell) 
CN.mat <- do.call(rbind, CN.mat)
euc.dist <- dist(CN.mat, method = "euclidean")
hc <- hclust(euc.dist, method = "ward.D2")

## Prepare data for plottting
CN.mat.df <- as.data.frame(CN.mat)
CN.mat.df$ID <- rownames(CN.mat.df)
rownames(CN.mat.df) <- NULL
CN.mat.df.ord <- CN.mat.df[hc$order,] #sort data based on clustering
plt.data <- melt(CN.mat.df.ord) #get long data format
plt.data$ID <- factor(plt.data$ID, levels=unique(plt.data$ID))
plt.data$y <- rep(0.5:length(levels(plt.data$ID)), length(unique(plt.data$variable)))

## Plot the data
#Plot Heatmap
pltlist <- list()
plt <- ggplot(plt.data) + geom_tile(aes(x=variable, y=y, fill=value), col='black') + theme_bw() + ylab("") + scale_y_continuous(breaks = plt.data$y, labels = plt.data$ID, expand = c(0,0))
plt <- plt + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=20), axis.line=element_blank(), axis.title.x=element_blank(), plot.margin=margin(l=-0.8,unit="cm"))
pltlist[['heatmap']] <- plt

#Plot Dendrogram
dhc <- stats::as.dendrogram(hc)
ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
ggdndr <- ggplot(ddata$segments) + geom_segment(aes(x=x, xend=xend, y=y, yend=yend)) + scale_y_continuous(breaks = plt.data$y, labels = plt.data$ID) + scale_y_reverse() + scale_x_continuous(expand=c(0,0))
ggdndr <- ggdndr + coord_flip()
ggdndr <- ggdndr + theme(panel.background=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.title=element_blank())
pltlist[['dendro']] <- ggdndr

## plot both dendrogram and heatmap together
cowplot::plot_grid(plotlist=rev(pltlist), align='h', ncol=length(pltlist), rel_widths = c(2,4))

## old dendrograms
#plot(hc)
#ggtree(as.phylo(hc))
