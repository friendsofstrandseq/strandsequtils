#' A plotting function for entire cells
#'
#' This function allows you to plot a whole cell view of your mosaiCatcher binned Strand-Seq data.
#' Plots all the cells in a data table, one cell per page, into a pdf-file.
#' Optional breakpoint file can be given to plot breakpoints/segmentation/SCEs/structural variants.
#'
#' Based on the QC plotting code by Sascha Meiers.
#'
#' @param d data table of data to be plotted. For example output from run_GFLars.
#' @param pdf_out filename for the PDF. Note: if the file exists already, it will be overwritten without warnings. example "drive/folder/subfolder/example.pdf"
#' @param bkp_file optional breakpoint locations. For example output from 'detect_SCEs'.
#' @param verbose print status messages. Defaults to FALSE.
#' @param HMM use HMM strand state information for background colours. Defaults to TRUE.
#' @return
#' @keywords Strand-Seq, mosaiCatcher, chromosome, sample, cell, plot
#' @author Venla Kinanen
#' @author Sascha Meiers
#' @export
#' @examples
#' plot_all_into_pdf()

plot_all_into_pdf <- function(d, pdf_out, bkp_file = NULL, verbose = F, HMM = T) {

  require(assertthat)
  require(data.table)
  require(ggplot2)
  require(scales)
  require(cowplot)
  library(grDevices)

  # Check that correct data are given:
  assertthat::assert_that("chrom" %in% colnames(d))
  assertthat::assert_that("start" %in% colnames(d) &&
                            is.integer(d$start))
  assertthat::assert_that("end" %in% colnames(d) &&
                            is.integer(d$end))
  assertthat::assert_that("sample" %in% colnames(d))
  assertthat::assert_that("cell" %in% colnames(d))
  assertthat::assert_that("w" %in% colnames(d) && is.integer(d$w))
  assertthat::assert_that("c" %in% colnames(d) && is.integer(d$c))
  assertthat::assert_that("class" %in% colnames(d))

  if (!data.table::is.data.table(d)) {
    d = data.table::as.data.table(d)
  }

  # Re-name and -order chromosomes
  d = d[, chrom := sub('^chr', '', chrom)][]
  d = d[grepl('^([1-9]|[12][0-9]|X|Y)$', chrom), ]
  d = d[, chrom := factor(chrom, levels = as.character(c(1:22, 'X', 'Y')), ordered = T)]
  chrom_sizes = d[, .(start=min(start), end=max(end)), by=chrom]

  if (!is.null(bkp_file)) {
    # Re-name and -order chromosomes
    bkp_file = bkp_file[, chrom := sub('^chr', '', chrom)][]
    bkp_file = bkp_file[grepl('^([1-9]|[12][0-9]|X|Y)$', chrom), ]
    bkp_file = bkp_file[, chrom := factor(chrom, levels = as.character(c(1:22, 'X', 'Y')), ordered = T)]
  }
  # Plot all cells
  cairo_pdf(pdf_out,
            width = 14,
            height = 10,
            onefile = T)

  for (s in unique(d$sample)) {
    # one pdf per cell
    for (ce in unique(d[sample == s, ]$cell)) {
      if (verbose) {
        message(paste("Plotting sample", s, "cell", ce, "into", pdf_out))
      }

      e = d[sample == s & cell == ce, ]

      if (!is.null(bkp_file)) {
        g = bkp_file[sample == s & cell == ce, ]
      }


      # Calculate some information
      info_binwidth = stats::median(e$end - e$start)
      info_reads_per_bin = stats::median(e$w + e$c)
      info_chrom_sizes = e[, .(xend = max(end)), by = chrom]
      info_num_bins = nrow(e)
      info_total_reads = sum(e$c + e$w)
      info_y_limit = 2 * info_reads_per_bin + 1
      info_sample_name = paste(s, ce, sep = "_")
      info_cell_name   = ce

      # main plot:
      plt <- ggplot2::ggplot(e) +
        ggplot2::aes(x = (start + end) / 2)

      if (HMM) {
        # prepare consecutive rectangles for a better plotting experience
        consecutive = cumsum(abs(diff(as.numeric(
          as.factor(e$class)
        ))))
        consecutive = c(consecutive[1], consecutive)
        e$consecutive = consecutive
        f = e[, .(
          start = min(start),
          end = max(end),
          class = class[1]
        ), by = .(consecutive, chrom)][]

        plt <- plt +
          ggplot2::geom_rect(
            data = f,
            aes(
              xmin = start,
              xmax = end,
              ymin = -Inf,
              ymax = Inf,
              fill = class
            ),
            inherit.aes = F,
            alpha = 0.25
          )
      }

      plt <- plt +
        ggplot2::scale_fill_manual(
          values = c(
            WW = "sandybrown",
            CC = "paleturquoise4",
            WC = "yellow",
            None = NA
          )
        )

      # Watson/Crick bars
      plt <- plt +
        suppressWarnings(
          ggplot2::geom_bar(
            aes(y = -w, width = (end - start)),
            stat = 'identity',
            position = 'identity',
            fill = 'sandybrown'
          )
        ) +
        suppressWarnings(
          ggplot2::geom_bar(
            aes(y = c, width = (end - start)),
            stat = 'identity',
            position = 'identity',
            fill = 'paleturquoise4'
          )
        ) +
        # Trim image to 2*median cov
        ggplot2::coord_flip(expand = F,
                            ylim = c(-info_y_limit, info_y_limit)) +
        ggplot2::facet_grid(. ~ chrom, switch = "x") +
        ggplot2::ylab("Watson | Crick") + xlab(NULL) +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(12),
                                    labels = format_Mb) +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(3)) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          # deprecated in ggplot2 2.2.0
          #panel.margin = unit(0, "lines"),
          panel.spacing = unit(0, "lines"),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(fill = NA, colour = NA)
        ) +
        ggplot2::guides(fill = FALSE) +

        # Dotted lines at median bin count
        ggplot2::geom_segment(
          data = info_chrom_sizes,
          aes(
            xend = xend,
            x = 0,
            y = -info_reads_per_bin,
            yend = -info_reads_per_bin
          ),
          linetype = "dotted",
          col = "darkgrey",
          size = 0.5
        ) +
        ggplot2::geom_segment(
          data = info_chrom_sizes,
          aes(
            xend = xend,
            x = 0,
            y = +info_reads_per_bin,
            yend = +info_reads_per_bin
          ),
          linetype = "dotted",
          col = "darkgrey",
          size = 0.5
        ) +
        ggplot2::geom_segment(
          data = info_chrom_sizes,
          aes(xend = xend, x = 0),
          y = 0,
          yend = 0,
          size = 0.5
        )

      # Breakpoints
      if(!is.null(bkp_file) && nrow(g) > 0) {

        plt <- plt +
          geom_vline(data = g,  aes(xintercept = start), color="orangered", size=1, alpha=0.6) +
          geom_vline(data = g,  aes(xintercept = end), color="orangered", size=1, alpha=0.6)
      }

      # Rename classes:
      labels = e[, .N, by = class][, label := paste0(class, " (n=", N, ")")][]

      e[, class := factor(class,
                          levels = labels$class,
                          labels = labels$label)]

      all <- cowplot::ggdraw() + cowplot::draw_plot(plt) +
        cowplot::draw_label(
          paste("Sample:", info_sample_name),
          x = .29,
          y = .97,
          vjust = 1,
          hjust = 0,
          size = 14
        ) +
        cowplot::draw_label(
          paste(
            "Median binwidth:",
            format(
              round(info_binwidth / 1000, 0),
              big.mark = ",",
              scientific = F
            ),
            "kb"
          ),
          x = .29,
          y = .93,
          vjust = 1,
          hjust = 0,
          size = 10
        ) +
        cowplot::draw_label(
          paste("Number bins:", format(info_num_bins, big.mark = ",")),
          x = .29,
          y = .91,
          vjust = 1,
          hjust = 0,
          size = 10
        ) +
        cowplot::draw_label(
          paste(
            "Total number of reads:",
            format(info_total_reads, big.mark = ",")
          ),
          x = .29,
          y = .89,
          vjust = 1,
          hjust = 0,
          size = 10
        ) +
        cowplot::draw_label(
          paste("Median reads/bin (dotted):", info_reads_per_bin),
          x = .29,
          y = .87,
          vjust = 1,
          hjust = 0,
          size = 10
        ) +
        cowplot::draw_label(
          paste0("Plot limits: [-", info_y_limit, ",", info_y_limit, "]"),
          x = .29,
          y = .85,
          vjust = 1,
          hjust = 0,
          size = 10
        )

      print(all)
    }
  }
  grDevices::dev.off()
}

format_Mb <- function(x) {
  paste(scales::comma(x / 1e6), "Mb")
}



# Command Line
library(data.table)

args = commandArgs(trailingOnly = T)
if(length(args)!=3 || !grepl('\\.pdf$', args[3])) {
  warning("Usage: Rscript plot_all_into_pdf.R count.table.gz SCEs.txt output.pdf")
  stop()
}
f_in = args[1]
f_sce = args[2]
f_out = args[3]

zcat_command = "zcat"
if(grepl('\\.gz$', f_in))
  f_in = paste(zcat_command, f_in)
d = fread(f_in)

sce = fread(f_sce)

plot_all_into_pdf(d, f_out, bkp_file = sce)
