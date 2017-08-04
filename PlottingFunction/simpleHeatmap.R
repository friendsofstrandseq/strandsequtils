require("ComplexHeatmap")
require("RColorBrewer")

plotHeatmapSegment2 <- function(dataFrameSegm, log=FALSE) {
  
  probs <- as.matrix(tab[,c(8:ncol(dataFrame))])
  
  colColors <- brewer.pal(n=6, name="Set1")
  names(colColors) <- c("CN0","CN1","CN2","CN3","CN4","CN5")
  
  annot1 <- HeatmapAnnotation(df = data.frame(state=dataFrame$types), which = "row", col = list(state = c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown")) )
  colAnnot.df <- data.frame(type = c("CN0","CN1","CN2","CN3","CN4","CN5", rep(c("CN0","CN1","CN2","CN3","CN4","CN5"), c(1,2,3,4,5,6))))
  annot2 <- HeatmapAnnotation(df = colAnnot.df, col=list(type=colColors)) 
  
  if (log) {
    probs.log <- log10(probs)
    probs.log[is.infinite(probs.log)] <- NA
    plt <- Heatmap(probs.log, name = "Probs", cluster_columns = FALSE, cluster_rows = T, show_row_names = FALSE, top_annotation = annot2)
  } else {
    plt <- Heatmap(probs, name = "Probs", cluster_columns = FALSE, cluster_rows = T, show_row_names = FALSE, top_annotation = annot2) 
  }
  
  heat.plt <- annot1 + plt
  heat.plt <- draw(annot1 + plt, row_dend_side = "left", row_sub_title_side = "right")
  return(heat.plt)
} 