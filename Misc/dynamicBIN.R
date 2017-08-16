

bamFilePath <- "/media/daewoooo/5D7CBC85372E82B0/HGSVC_project/PUR/HG00733/"
chromosomes <- paste0('chr', c(1:22))

dynamicBin <- function(bamFilePath, min.mapq=10, bin.length=1000000, step=1000000, chromosomes = NULL) {
  
  bamFiles <- list.files(path, full.names=TRUE, recursive=FALSE, pattern='.bam$')
  
  file.header <- Rsamtools::scanBamHeader(bamFiles[1])[[1]]
  chrom.lengths <- file.header$targets
  chroms.in.data <- names(chrom.lengths)
  if (is.null(chromosomes)) {
    chromosomes <- chroms.in.data
  }
  
  chroms2use <- intersect(chromosomes, chroms.in.data)
  gr <- GenomicRanges::GRanges(seqnames=Rle(chroms2use), ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
  
  cumCov.perLib <- GRanges()
  for (bam in bamFiles) {
    bamName <- basename(bam)
    message("Loading data for lib: ", bamName)
    suppressWarnings( raw.reads <- GenomicAlignments::readGAlignmentPairs(bam, param=Rsamtools::ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F))) )

    first <- as(GenomicAlignments::first(raw.reads), 'GRanges')
    second <- as(GenomicAlignments::second(raw.reads), 'GRanges')
    data <- sort(c(first, second))
      
    if (!is.null(min.mapq)) {
      if (any(is.na(mcols(data)$mapq))) {
        warning(paste0(file,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
        mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
      }
      data <- data[mcols(data)$mapq >= min.mapq]
    }
  
    if (length(cumCov.perLib)==0) {
      cumCov.perLib <- reduce(data)
    } else {
      cumCov.perLib <- reduce(c(cumCov.perLib, data[,0]))
    }
  }
  cumCov.perLib <- keepSeqlevels(cumCov.perLib, chroms2use)
  cov <- coverage(cumCov.perLib)
  
  dynamic.bins <- GRangesList()
  for(chr in names(cov)) {
    chr.cov <- cov[[chr]]
    cs <- cumsum(chr.cov)
   
    bin.starts <- seq(from = 1, to = max(cs)-bin.length, by = step)
    bin.ends <- seq(from = bin.length, to = max(cs), by = step)
    bins <- c(rbind(bin.starts, bin.ends)) 
    
    gen.pos <- findInterval(bins, as.vector(cs), rightmost.closed = T)
    gen.bins <- reformat(gen.pos)
    gen.bins.gr <- GRanges(seqnames = chr, ranges=IRanges(start=gen.bins$start, end=gen.bins$end))
    dynamic.bins[[chr]] <- gen.bins.gr
  }
  seqlengths(dynamic.bins) <- chrom.lengths[chroms2use] 
  return(unlist(dynamic.bins))
}


reformat <- function(x) {
  out_list <- list() 
  for ( i in seq(1, length(x), 2) ) {
    out_list[[i]] <- c(x[i], x[i+1])
  }
  mt <- do.call("rbind",out_list)
  df <- data.frame(mt)
  colnames(df) <- c("start", "end")
  df
}
