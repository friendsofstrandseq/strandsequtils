#' Detect SCEs
#'
#' Function for the detection of sister chromatid exchanges. Build on the HMM classification by Sascha Meiers.
#'
#' @param d data table.
#' @param min_SCE size cutoff in MB. Only segments larger than the cutoff will be considered as potential SCEs. Defaults to 5MB.
#' @param rec_perc how many cells should a recurrent potential SV be to be counted an SV instead of an SCE percentage. 0 means at least one other cell, 100 means it should be present in all the cells. note that if there are several individuals in the data set, setting a high recurrence minimun might mean an SV would have to be shared between individuals
#' @param qc_fpb quality cutoff percentage for fragments per bin per segment. Only segments with medianan W + C reads/bin exceeding the cutoff percentage will be considered. Note that bin size affects fragments per bin, and the default value of 20 can be low for larger bin sizes.
#' @param rec_offset how many bin sizes can the potential SV/SCE be from the exact start and end. 0 means exactly the same start and end, 1 means 1 binsize either way for both the start and end, etc.
#' @param verbose print status messages. Defaults to FALSE.
#' @return data table with the locations of the potential SCEs.
#' @author Venla Kinanen
#' @keywords Strand-Seq, mosaiCatcher, SCE, chromosome, sample, cell
#' @seealso plot_all_into_pdf()
#' @export
#' @examples
#' detect_SCEs()
detect_SCEs <- function(d, min_SCE=5, rec_perc=2, qc_fpbp=50, rec_offset=2, verbose=FALSE) {

  require(assertthat)
  require(data.table)

  ## Check that correct data are given:
  assertthat::assert_that( "chrom" %in% colnames(d) )
  assertthat::assert_that( "start" %in% colnames(d) && is.integer( d$start ) )
  assertthat::assert_that( "end" %in% colnames(d) && is.integer( d$end ) )
  assertthat::assert_that( "sample" %in% colnames(d) )
  assertthat::assert_that( "cell" %in% colnames(d) )
  assertthat::assert_that( "w" %in% colnames(d) && is.integer( d$w ) )
  assertthat::assert_that( "c" %in% colnames(d) && is.integer( d$c ) )
  assertthat::assert_that(assertthat::is.number( rec_perc ) && rec_perc>=0 && rec_perc<=100 )
  assertthat::assert_that(assertthat::is.number( qc_fpbp ) && qc_fpbp>0 &&  qc_fpbp<=100)
  assertthat::assert_that(assertthat::is.number( rec_offset ) && rec_offset>=0 )

  if (!( data.table::is.data.table(d) )) {
    d = data.table::data.table(d)
  }

  # Not looking at empty bins
  d <- d[class != "None",]

  # Just for convenience. A unique identifier for each cell. Unnessary in cases where the 'cell' field is enough for identification.
  if(is.null(d$cell_name)) {
    d$cell_name = paste(d$sample, d$cell, sep = "_")
  }

  # Divide into segments, calculate segment length and distance from chromosome start/end, majority type
  #--- Sascha's code
  # Annotate majority class
  e = d
  majority_types = e[, .N, by = .(sample, cell, chrom, class)][order(sample, cell, chrom),][, .(majority = class[which.max(N)]), by = .(sample, cell, chrom)]

  e = merge(e, majority_types, by = c("sample", "cell", "chrom"))

  # combine consecutive bins into intervals
  e = e[order(sample, cell, chrom, start, end),]
  e$cnsc <- cumsum(e[, .(consecutive = c(1,abs(diff(as.numeric(factor(class)))))), by = .(sample, cell, chrom)]$consecutive)

  # determine length of each consecutive element
  lengths = e[, .(cnsc_len = max(end) - min(start), cnsc_class = class[1]), by = cnsc]
  e = merge(e, lengths, by = "cnsc")

  # determine distance of each consecutive element to the chromosome start / end
  chrom_sizes = e[, .(chr_start = min(start), chr_end = max(end)), by = chrom]
  e = merge(e, merge(e, chrom_sizes, by = "chrom")[, .(cnsc_dist = min(c(start - chr_start, chr_end - end))), by = cnsc], by = "cnsc")
  d = e
  #--- End Sascha's code

  ## Prepare data
  # Order chromosomes, initialize fields
  ##! this causes NAs in samples with non-standard chromosomes. Also not really necessary?
  ##!if (suppressWarnings(!is.na(as.numeric(as.character(d$chrom[1]))))) {
  ##!  d = d[, chrom := factor(chrom, levels=as.character(c(1:22,'X','Y')), ordered = T)]
  ##!} else {
  ##!  d = d[, chrom := factor(chrom, levels=as.character(c(paste0("chr", 1:22),'chrX','chrY')), ordered = T)]
  ##!}

  # These are not actually passed on to the user at the moment, could well be purged after debugging is done.
  d$recurrent = rep(FALSE)
  d$SCE = rep(NA)
  d$QC = rep("")

  # How many cells should a recurrent potential SV be to be counted an SV instead of an SCE
  # percentage. 0 means at least one other cell, 100 means it should be present in all the cells.
  # note that if there are several individuals in the data set, setting a high recurrence minimun
  # might mean an SV would have to be shared between individuals
  if (rec_perc > 0) {
    min_recurrence = round(length(unique(d$cell_name))*rec_perc/100)
  } else {
    min_recurrence = 1
  }

  bin_size = (d$end[1] - d$start[1])
  chrom_sizes = d[, .(start=min(start), end=max(end)), by=chrom]

  # For permitted strand changes for an SCE. Only one strand should change, so CC -> WW or WW -> CC would not be an SCE.
  permitted = F
  dimn = c("WW", "WC", "CC")
  states = matrix(c(0,1,0, 1,0,1, 0,1,0), nrow = 3, ncol = 3, byrow = T, dimnames = list(dimn, dimn))

  ## Create breakpoint table
  # Select first and last rows of every segment
  bkps = d[, .SD[c(1,.N)], by=cnsc]
  # Name first rows of segments "start", and last rows "end".
  bkps = bkps[, id:= c("start","end")]
  # Combine first and last rows of every segment into just one row per segment. 2 to 1.
  bkps = dcast(bkps, cnsc + sample + cell + class + majority + cnsc_len + chrom + cell_name ~ id, value.var =  c("start", "end"))
  # Acual segment start and segment end, instead of bin start and end
  bkps = bkps[, start:=start_start][, end:=end_end]
  bkps = bkps[, c("start_start", "start_end", "end_start", "end_end"):= NULL]
  # For ease of indexing
  bkps_chrom = data.table::data.table( bkps )
  data.table::setkey( bkps, cnsc )
  data.table::setkey( bkps_chrom, chrom, start )
  data.table::setkey( d, cnsc, cnsc_class )

  # Go through all the segments
  for (j in unique( bkps$cnsc )) {

    # This is the segment we are examining
    current = bkps[ cnsc == j ][1]

    # Skip if empty
    if ( dim( current )[1]==0 ) next

    # Skip if segment is the same class/type as the majority class/type of the chromosome
    if ( current$class == current$majority ) {
      data.table::set( d, which(d$cnsc==j), "SCE", FALSE )
      next
    }

    ## Length
    # Not an SCE if length less than parameter (in MB)
    if ( current$cnsc_len < ( min_SCE*1e6 ) ) {
      data.table::set( d, which( d$cnsc==j ), "cnsc_class", current$majority )
      data.table::set( d, which( d$cnsc==j ), "SCE", FALSE )
      next
    }

    ## Strand state change
    # Compare to next segment, unless this is the last segment
    if ( j < max( bkps$cnsc ) ) {
      fwd = bkps[ bkps$cnsc == (j+1) ][1]
      # Are this and next segment on the same chromosome
      if (!is.na(current$chrom) && !is.na(fwd$chrom) && current$chrom == fwd$chrom ) {
        # Is the strand state switch one that is possible for an SCE
        if ( (fwd$class == fwd$majority) && states[fwd$class, current$class] ) permitted=T
      }
    }

    # Compare to previous segment, unless this is the first segment
    if ( j > 1 ) {
      prev = bkps[ bkps$cnsc == (j-1) ][1]
      # Are this and the previous segment on the same chromosome
      if ( !is.na(current$chrom) && !is.na(prev$chrom) && current$chrom == prev$chrom ) {
        # Is the strand state switch one that is possible for an SCE
        if ( (prev$class == prev$majority) && states[prev$class, current$class] ) permitted=T
      }
    }

    # Strand state switches not compatible with SCE, probably an SV
    if ( !permitted ) {
      data.table::set( d, which( d$cnsc==j ), "cnsc_class", value=current$majority )
      data.table::set( d, which( d$cnsc==j ), "SCE", value=FALSE )
      # Compatible, leave class as is
    } else {
      permitted = F
      data.table::set( d, which( d$cnsc==j ), "SCE", value=TRUE )
    }


    ## Recurrence.
    reads_per_bin = d[ cnsc == current$cnsc ][, median(w + c) ]
    reads_per_chromosome = d[ chrom == d[cnsc == current$cnsc]$chrom[1] &
                                cell_name == d[cnsc == current$cnsc]$cell_name[1]][, median(w + c) ]

    # Require minimun fragments per bin.
    if ( reads_per_bin > (reads_per_chromosome * qc_fpbp/100) ) {
      # Subset segment table to segment start within +/- offset times bin size (wobble room)
      t = bkps[ chrom==current$chrom & start %in% seq(from=(current$start-rec_offset*bin_size), to=(current$start+rec_offset*bin_size), by=bin_size) ]
      # And further to segment end within +/- offset times bin size
      t = t[ end %in% seq(from=(current$end-bin_size), to=(current$end+bin_size), by=bin_size)]
      # Recurrence frequency over minimum limit, and class not the same as chrom majority
      if (( dim(t)[1]>min_recurrence ) && ( current$class != current$majority )) {
        data.table::set( d, which( d$cnsc==j ), "cnsc_class", value=current$majority )
        data.table::set( d, which( d$cnsc==j ), "SCE", value=FALSE )
        data.table::set( d, which( d$cnsc==j ), "recurrent", value=TRUE )
      }
      # Ignore segments with very low fragment count
    } else {
      data.table::set( d, which( d$cnsc==j ), "SCE", value=FALSE )
      data.table::set( d, which( d$cnsc==j ), "cnsc_class", value=current$majority )
      data.table::set( d, which( d$cnsc==j ), "QC", value="low" )
    }
  }

  ### Return original data table with additional columns.
  ### This was for debugging. Getting rid of debug bloat might be useful now.
  ###return( d )

  # Turn segments into breakpoints. This is hacky.
  # choose first and last rows of these columns (getting rid of bins)
  t = d[, .(start = min(start), end = max(end)), by = .(sample, cell, chrom, class = cnsc_class)]
  return(t)
}


library(data.table)

# Command Line
args = commandArgs(trailingOnly = T)
if(length(args)!=2) {
  warning("Usage: Rscript SCE.R count.table.gz SCEs.txt")
  stop()
}
f_in = args[1]
f_out = args[2]
zcat_command = "zcat"
if(grepl('\\.gz$', f_in))
  f_in = paste(zcat_command, f_in)
d = fread(f_in)

sces = detect_SCEs(d)
sces = sces[, .(chrom, start, end, sample, cell, class)]
write.table(sces, file = f_out, col.names = T, row.names = F, quote=F, sep="\t")
