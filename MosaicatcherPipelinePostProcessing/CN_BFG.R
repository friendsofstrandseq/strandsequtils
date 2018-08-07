## Load required libraries
library(data.table)
library(GenomicAlignments)

bin.size <- 100000
maximumCN <- 1000
sample.id <- "C7_data"
regions <- "/media/porubsky/Elements/StrandSeqNation/C7/BFB_CN_estimates/BFBregion_C7_CNV_segments.txt"
dir <- "/media/porubsky/Elements/StrandSeqNation/C7/BFB_CN_estimates/"
outputfolder <- dir
bam.folder <- "/media/porubsky/Elements/StrandSeqNation/C7/selected_bam/"

# get the file names and the directory paths of input files (Mosaicatcher outputs)
countsFile <- file.path(dir, paste0(format(bin.size, scientific = F),"_fixed.txt.gz"))
infoFile <- file.path(dir, paste0(format(bin.size, scientific = F),"_fixed.info"))
stateFile <- file.path(dir, "final.txt")

# read the input files
counts <- fread(paste("zcat", countsFile))
segs <- getCountsPerRegions(bam.folder=bam.folder, outputfolder=outputfolder, regions=regions, sample.id=sample.id)
segs <- as(segs, 'data.table')
info <- fread(infoFile)
strand <- fread(stateFile)

# source mosaiclassifier or source manually if this doesn't work
file.sources = list.files(file.path(dir, "mosaiClassifier/"), pattern="*.R", full.names = TRUE)
for (file in file.sources) {
  source(file)
}

## Run get_CN_NB_prob function
get_CN_NB_prob(counts=counts, segs=segs, info=info, strand=strand, manual.segs=TRUE, confidence.interval=0.95)
## Load produced table
tab2 <- read.table("/media/porubsky/Elements/StrandSeqNation/C7/BFB_CN_estimates/100000_BFB_cell_confident_CN.table", header=T, stringsAsFactors = FALSE) 
## Plot the data
plotCNheatmap(tab2) -> hm.plot2 #CN estim


## FUNCTIONS to be loaded ##
############################

#' This funcion will go through all BAM files present in a 'bam.folder' and will
#' counts directional reads in specified 'regions'.
#' 
#' NOTE: Works only for paired-end reads!!!
#' 
#' @param bam.folder A folder where are stored BAMs.
#' @param outputfolder A folderr where the final counts will be stored.
#' @param regions A bed file with regions to count reads in.
#' @param sample.id A sample identifier.
#' @author David Porubsky
#' 
#' ## To run the function
#' sample.id <- "C7_data"
#' bam.folder <- "/media/porubsky/Elements/StrandSeqNation/C7/"
#' outputfolder <- "/media/porubsky/Elements/StrandSeqNation/C7/BFB_CN_estimates/"
#' regions <- "/media/porubsky/Elements/StrandSeqNation/C7/BFB_CN_estimates/BFBregion_C7_CNV_segments.txt"
#' getCountsPerRegions(bam.folder=bam.folder, outputfolder=outputfolder, regions=regions, sample.id=sample.id) -> count.table

getCountsPerRegions <- function(bam.folder=NULL, outputfolder=NULL, regions=NULL, sample.id=NULL) {
  
  ## Load required packages
  if ("GenomicAlignments" %in% rownames(installed.packages())) {
    require(GenomicAlignments)
  } else {
    stop("Package 'GenomicAlignments' cannot be loaded!!!")
  } 
  
  BFB.regions <- read.table(regions, header=FALSE, stringsAsFactors = FALSE)
  BFB.regions.gr <- GenomicRanges::GRanges(seqnames=BFB.regions$V1, ranges=IRanges(start=BFB.regions$V2, end=BFB.regions$V3))
  
  bams <- list.files(bam.folder, pattern = "\\.bam$", full.names = TRUE)
  
  file.header <- Rsamtools::scanBamHeader(bams[1])[[1]]
  chrom.lengths <- file.header$targets
  seqlengths(BFB.regions.gr) <- chrom.lengths[seqlevels(BFB.regions.gr)]
  
  #allReads <- GRangesList()
  binned.count.grl <- GRangesList()
  for (bam in bams) {
    filename <- basename(bam)
    cell <- unlist(strsplit(filename, "\\."))[1]
    message("Processing ", filename, " ...")
    ## Load data (paired-end only)
    suppressWarnings( data.raw <- GenomicAlignments::readGAlignmentPairs(file = bam, param=Rsamtools::ScanBamParam(which=range(BFB.regions.gr), what='mapq', flag=Rsamtools::scanBamFlag(isDuplicate=FALSE))) )
    ## Use only first mate of the pair
    data <- as(GenomicAlignments::first(data.raw), "GRanges")
    ## filter by mapping quality 10
    data <- data[mcols(data)$mapq >= 10]
    ## Count reads and creat bins object
    bins <- BFB.regions.gr
    bins$sample <- sample.id
    bins$cell <- cell
    bins$C <- countOverlaps(bins, data[strand(data) == "+"])
    bins$W <- countOverlaps(bins, data[strand(data) == "-"])
    bins$class <- "None"
    binned.count.grl[[bam]] <- bins
  } 
  
  binned.count.gr <- unlist(binned.count.grl, use.names = FALSE)
  binned.count.gr <- sort(binned.count.gr)
  count.table.df <- as.data.frame(binned.count.gr)
  colnames(count.table.df)[1] <- 'chrom'
  remove.cols <- which(colnames(count.table.df) %in% c('width','strand'))
  count.table.df <- count.table.df[,-remove.cols]
  
  destination <- file.path(outputfolder, paste0(sample.id, ".RegionCounts.txt"))
  write.table(count.table.df, file=destination, quote = FALSE, row.names = FALSE, sep = "\t")
  
  return(count.table.df)
}


#' This funcion calculate CN estimates for each user defined segment
#' @param ...
#' @author Maryam Ghareghani

get_CN_NB_prob <- function(counts, segs, info, strand, manual.segs=TRUE, confidence.interval=0.95) {
  
  # set the None classes to WW in the counts data table
  counts[class=="None", class:="WW"]
  # prepare the data to be used in SV classification
  probs <- mosaiClassifierPrepare(counts, info, strand, segs, manual.segs = manual.segs)
  probs[, nb_r:=expected*nb_p/(1-nb_p)]
  
  probs <- probs[, cbind(.SD, CN=1:maximumCN), 
                 by=.(sample, cell, chrom, start, end)]
  
  
  CNprobs <- probs[, .(cov=sum(W+C), 
                       expected=sum(expected), 
                       nb_p=nb_p[1], 
                       nb_r=sum(nb_r)), 
                   by=.(sample, chrom, start, end)]
  # first observation
  (CNprobs$cov / CNprobs$expected) * 2
  
  CNprobs <- CNprobs[, cbind(.SD, CN=1:maximumCN), 
                     by=.(sample, chrom, start, end)]
  
  
  # computing NB probs and calling most probable CN in probs table per segment per cell
  probs[, CN_ll:=dnbinom(x=C+W, size=CN*(nb_r/2), prob=nb_p)]
  cellCNcalls <- probs[, .SD[which.max(CN_ll)], 
                       by=.(sample, cell, chrom, start, end)]
  
  # computing NB probs and calling CN in CNprobs table per segment
  CNprobs[, CN_ll:=dnbinom(x=cov, size=CN*(nb_r/2), prob=nb_p)]
  CNcalls <- CNprobs[, .SD[which.max(CN_ll)], 
                     by=.(sample, chrom, start, end)]
  
  
  # removing extra columns for cleaning the data for plotting phylogenetic tree
  cellCNcalls <- cellCNcalls[start!=0]
  cellCNcalls[, `:=`(nb_p=NULL, class=NULL, C=NULL, W=NULL,expected=NULL, scalar=NULL, nb_r=NULL, CN_ll=NULL)]
  
  # normalize probs (all CN probs for a single-cell should sum up to 1)
  probs[, CN_ll:=CN_ll/sum(CN_ll), by=.(sample, chrom, start, end, cell)]
  
  # sort probs by CN_ll in each cell
  setorder(probs, sample, chrom, start, end, cell, -CN_ll)
  
  # add an extra column for the sum of CN probabilities
  probs[, cumsum_CN_ll:=cumsum(CN_ll),
        by = .(sample, chrom, start, end, cell)]
  
  # shrink the table to have only the parts of the table that fall into the confidence interval
  conf.probs <- probs[, num_conf_CNs:=min(length(which(cumsum_CN_ll < confidence.interval))+1, .N),
                           by = .(sample, chrom, start, end, cell)]
  
  conf.probs <- conf.probs[, head(.SD, num_conf_CNs[1]),
                           by = .(sample, chrom, start, end, cell)]
  
  # define the table that includes the most likely and the confidence interval range for each segment and cell
  conf.CNs <- conf.probs[, cbind(head(.SD,1), min_CN=min(CN), max_CN=max(CN)),
                         by = .(sample, chrom, start, end, cell)]
  
  conf.CNs[, `:=`(CN_ll=NULL, cumsum_CN_ll=NULL, num_conf_CNs=NULL)]
  
  write.table(conf.CNs, 
              file.path(dir, paste0(format(bin.size, scientific = F),"_BFB_cell_confident_CN.table")),
              sep = "\t",
              quote = F)
  
  write.table(cellCNcalls, 
              file.path(dir, paste0(format(bin.size, scientific = F),"_BFB_cell_CNs.table")),
              sep = "\t",
              quote = F)
  
  write.table(probs, 
              file.path(dir, paste0(format(bin.size, scientific = F),"_BFB_cell_CN_probs.table")),
              sep = "\t",
              quote = F)
  
}


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
  
  ## Reorder HC, put zero cells on the top of the figure
  new.order <- hc$order-2
  new.order[new.order==-1] <- 154
  new.order[new.order==0] <- 153
  hc <- reorder(hc, wts = new.order)
  
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
  if (continuous) {
    plt <- ggplot(plt.data) + geom_tile(aes(x=variable, y=y, fill=value)) + theme_bw() + ylab("") + scale_y_continuous(breaks = plt.data$y, labels = plt.data$ID, expand = c(0,0))
    plt <- plt + scale_fill_gradientn(colors = brewer.pal(name = 'BuPu', n = 9), name="CN")
    plt <- plt + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  } else {
    plt <- ggplot(plt.data) + geom_tile(aes(x=variable, y=y, fill=factor(ColCategs)), col='black') + theme_bw() + ylab("") + scale_y_continuous(breaks = plt.data$y, labels = plt.data$ID, expand = c(0,0))
    plt <- plt + scale_fill_manual(values = color.palette, name = "CN", labels = label.text)
    plt <- plt + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }
  plt <- plt + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(face="bold", size=10), axis.line=element_blank(), axis.title.x=element_blank(), plot.margin=margin(l=-0.8,unit="cm"))
  pltlist[['heatmap']] <- plt
  
  #Plot Dendrogram
  dhc <- stats::as.dendrogram(hc)
  ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
  ggdndr <- ggplot(ddata$segments) + geom_segment(aes(x=x, xend=xend, y=y, yend=yend)) + scale_y_continuous(breaks = plt.data$y, labels = plt.data$ID) + scale_x_continuous(expand=c(0,0.5))
  suppressMessages( ggdndr <- ggdndr + scale_y_reverse() )
  ggdndr <- ggdndr + coord_flip()
  ggdndr <- ggdndr + theme(panel.background=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank(), axis.title=element_blank())
  pltlist[['dendro']] <- ggdndr
  
  #Plot CN dotplot
  ord <- hc$order
  plt.df <- split(data.tab, data.tab$start) 
  for (i in 1:length(plt.df)) {
    df <- plt.df[[i]]
    df <- df[ord,]
    df$x <- 1:nrow(df)
    plt.df[[i]] <- df
  }
  plt.df <- do.call(rbind, plt.df)
  plt <- ggplot(plt.df)+ geom_errorbar(aes(ymin=min_CN, ymax=max_CN, x=plt.data$y))
  plt <- plt + geom_point(aes(x=plt.data$y, y=CN, color=factor(start))) + xlab("Cells")
  plt <- plt + scale_color_manual(values = brewer.pal(n = 6, name = "Set1"))
  plt <- plt + coord_flip() + scale_y_reverse() + facet_grid(. ~ factor(start), scales = 'free')
  plt <- plt + scale_x_continuous(breaks = unique(plt.data$y), expand = c(0,0))
  plt <- plt + theme(panel.grid.major.y = element_line(color = "grey80"), axis.text.y=element_blank(), strip.text.x = element_blank()) + xlab("")
  
  ## plot both dendrogram and heatmap together
  #complete.plt <- cowplot::plot_grid(plotlist=rev(pltlist), align='h', ncol=length(pltlist), rel_widths = c(2,4))
  complete.plt <- ggarrange(pltlist[[2]], pltlist[[1]], plt, ncol = 3, widths = c(1, 1.5, 4))
  return(complete.plt)
}