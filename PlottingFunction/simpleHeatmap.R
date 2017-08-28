plotHeatmapSegment <- function(dataFrame, plot.log=FALSE, file=NULL) {
  probs <- as.matrix(dataFrame[,c(8:ncol(dataFrame))])
  
  #cluter probs hclust
  #get order of rows
  #sort dataFrame rows based on hclust order
  ord <- order.dendrogram(as.dendrogram(hclust(dist(probs, method = "euclidean"), method = "ward.D")))
  dataFrame <- dataFrame[ord,]
  probs <- probs[ord,]
  
  dataFrame$cells <- factor(dataFrame$cells, levels=dataFrame$cells)
  
  tab.long <- melt(dataFrame, id.vars=c('cells', 'types', 'Wcount', 'Ccount','chr'), measure.vars=c("CN0","CN1","CN2","CN3","CN4","CN5","X00","X01","X10","X02","X11","X20","X03","X12","X21","X30","X04","X13","X22","X31","X40","X05","X14","X23","X32","X41","X50"))
  
  
  heatmap_theme <- theme(
    #legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    plot.margin = unit(c(-0.5,-0.5,-0.5,-0.5),"mm"),
    legend.position="bottom"
    #legend.text = element_text(size = 5),
    #legend.key.size = unit(0.5, "cm")
    #axis.ticks.x=element_blank()  
  )
  
  plt <- ggplot(tab.long) + geom_tile(aes(x=variable, y=cells, fill=as.numeric(value))) + heatmap_theme + scale_fill_continuous(name="")
  
  colColors <- brewer.pal(n=6, name="Set1")
  names(colColors) <- c("CN0","CN1","CN2","CN3","CN4","CN5")
  colAnnot.df <- data.frame(ID=factor(levels(tab.long$variable), levels=levels(tab.long$variable)), type = c("CN0","CN1","CN2","CN3","CN4","CN5", rep(c("CN0","CN1","CN2","CN3","CN4","CN5"), c(1,2,3,4,5,6))))
  
  header_theme <- theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.length = unit(0,"null"),
    plot.margin = unit(c(-0.5,-0.5,-0.5,-0.5),"mm"),
    legend.position="top"
  )
  
  header <- ggplot(colAnnot.df) + geom_tile(aes(x=ID, y=1, fill=type)) + scale_fill_manual(values = colColors) + header_theme + guides(fill = guide_legend(nrow = 1))
  
  g1 <- ggplotGrob(header)
  g2 <- ggplotGrob(plt)
  
  g <- rbind(g1, g2, size = "last")
  
  setPanelHeights <- function(g, heights){
    g$heights <- grid:::unit.list(g$heights)
    id_panels <- unique(g$layout[g$layout$name=="panel", "t"])
    g$heights[id_panels] <- heights
    g
  }
  
  g <- setPanelHeights(g, unit.c(unit(1,"line"), unit(nrow(g2),"line")))
  
  grid.newpage()
  grid.draw(g)
  
}  
