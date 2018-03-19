
plotReadDist <- function(count.tab=NULL, region=NULL, order=NULL){
  #order table
  #count.tab$cell <- as.numeric(gsub(count.tab$cell, pattern = "cell_", replacement = "")) + 1
  count.tab$cell <- bamNames$cell_id[match(count.tab$cell, bamNames$cell_name)]
  count.tab.ord <- count.tab[which(count.tab$cell %in% as.numeric(order)),]
  count.tab.ord <- count.tab.ord [order(match(count.tab.ord $cell, as.numeric(order))),] 
  
  count.tab.gr <- GRanges(seqnames=count.tab.ord$chrom, ranges=IRanges(start=count.tab.ord$start, end=count.tab.ord$end), W=count.tab.ord$w, C=count.tab.ord$c, cell=count.tab.ord$cell)
  region.lookup <- region
  start(region.lookup) <- start(region.lookup) - width(region)*2
  end(region.lookup) <- end(region.lookup) + width(region)*2
  
  count.tab.gr.sub <- subsetByOverlaps(count.tab.gr, region.lookup)
  
  count.tab.gr.sub.df <- as.data.frame(count.tab.gr.sub)
  count.tab.gr.sub.df$cell <- factor(count.tab.gr.sub.df$cell, levels=unique(count.tab.gr.sub.df$cell))
  plt.crick <- ggplot(count.tab.gr.sub.df) + geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=C), fill="paleturquoise4") + facet_grid(cell ~ .)
  plt.watson <- plt.crick + geom_rect(data=count.tab.gr.sub.df, aes(xmin=start, xmax=end, ymin=0, ymax=-W), fill="sandybrown")  + facet_grid(cell ~ .)
  plt.final <- plt.watson + geom_vline(xintercept = c(start(region), end(region)), color="red")
  
  return(plt.final)
}
