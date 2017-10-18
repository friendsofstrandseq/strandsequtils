# Shiny app for SCE classification
# global.R Helper functions

library(shiny)
library(shinythemes)
library(htmlwidgets)
library(data.table)
library(ggplot2)
library(shinyjs)

# Pretty printing, courtesy of Sascha Meiers
format_Mb <- function(x) {
  paste(scales::comma(x/1e6), "Mb")
}

# Data initialization, adapted from Sascha Meiers. Some redundance to handle slighly different mosaiCatcher runs.
initialize_data <- function(e) {
  e <- e[class != "None",]

  # Order chromosomes.
  e <- e[, chrom := sub('^chr', '', chrom)][] # get rid of 'chr' if is there
  e <- e[grepl('^([1-9]|[12][0-9]|X|Y)$', chrom), ] # only looking at standard chromosomes
  e$chrom <- paste0("chr", e$chrom) # add 'chr' for clarity on plots
  # Turn the 'chrom' field into a factor, and order it.
  e <- e[, chrom := factor(chrom, levels=as.character(c(paste0("chr", 1:22),'chrX','chrY')), ordered = T)]

  # A unique identifier for individual cells
  if (is.null(e$cell_name)) e$cell_name = paste(e$sample, e$cell, sep="_")

  # Find consecutive intervals and number them
  e <- e[order(cell_name, chrom, start, end),]
  e$cnsc <- cumsum(e[, .(consecutive = c(1,abs(diff(as.numeric(factor(class)))))), by = .(sample, cell, chrom)]$consecutive)

  return(e)
}

# Read result files for leaderboard.
get_results <- function(path, pattern) {
  files = NULL
  # List all result files for this dataset.
  files = list.files(path, pattern = paste0(pattern, ".txt"))
  results = data.table(cell_name = character(), cell = character(), sample = character(), chrom = character(),
                       cnsc = character(), end = integer(), user = character(), time = character(), counter = integer())
  if (!is.null(files) && length(files) > 0) {
    for (i in files) {
      a = fread(paste0(path, i), showProgress = F)
      results = rbind(results, a)
    }
  }
  return(results)
}

# Main plot with the watson and crick bars. Adapted from Sascha Meiers.
main_plot <- function(e) {
  if (the_end) {
    message <- paste0("\n   Congratulations!\n",
                 "       You have clicked through all of the chromosomes.\n")
    plt <- ggplot() +
      annotate("text", x=10, y=10, label = message) +
      theme_classic() +
      theme(panel.spacing = unit(0, "lines"),
            strip.background = element_rect(fill = NA, colour=NA),
            rect = element_blank()
      ) +
      guides(fill = FALSE)

    return(plt)
  }

  info_reads_per_bin = stats::median(e$w + e$c)
  info_y_limit = 2*info_reads_per_bin+1

  chr_end = max(e$end)

  plt <- ggplot(e) + aes(x = (start+end)/2)

  # prepare consecutive rectangles for a better plotting experience
  consecutive = cumsum(abs(diff(as.numeric(as.factor(e$class)))))
  consecutive = c(consecutive[1], consecutive)
  e$consecutive = consecutive
  f = e[, .(start = min(start), end = max(end), class = class[1]), by = .(consecutive, chrom, sample, cell)][]

  plt <- plt +
    geom_rect(data = f, aes(xmin = start, xmax=end, ymin=-Inf, ymax=Inf, fill=class),
              inherit.aes=F, alpha=0.25) +
    scale_fill_manual(values = c(WW = "sandybrown", CC = "paleturquoise4",
                                 WC = "yellow", None = NA))

  # Watson/Crick bars
  plt <- plt +
    suppressWarnings(geom_bar(aes(y = -w, width=(end-start)),
                              stat='identity', position = 'identity', fill='sandybrown')) +
    suppressWarnings(geom_bar(aes(y = c, width=(end-start)),
                              stat='identity', position = 'identity', fill='paleturquoise4')) +

    ylab("Watson | Crick") + xlab(NULL) +
    scale_x_continuous(breaks = scales::pretty_breaks(12), labels = format_Mb) +
    scale_y_continuous(breaks = scales::pretty_breaks(1)) +
    coord_cartesian(xlim = c(0, chr_end), ylim = c(-info_y_limit, info_y_limit)) +
    theme_classic() +
    theme(panel.spacing = unit(0, "lines"),
          strip.background = element_rect(fill = NA, colour=NA),
          plot.title = element_text(hjust = 0.5)) +
    guides(fill = FALSE)

  plt <- plt + facet_grid(sample+cell ~ .) +
    ggtitle(e$chrom[1], subtitle = NULL)

  return(plt)
}

# Overlay plot for just the breakpoints
bkp_plot <- function(e, bkps, highlighted) {
  if (the_end) {
    plt <- ggplot() + geom_blank() +
      theme_classic() +
      theme(panel.spacing = unit(0, "lines"),
            strip.background = element_rect(fill = NA, colour=NA),
            rect = element_blank()
      ) +
      guides(fill = FALSE)
    return(plt)
  }

  info_reads_per_bin <- round(median(e$w + e$c))
  info_y_limit <- 2*info_reads_per_bin+1

  chr_end <- max(e$end)

  bkps <- bkps[cnsc %in% e$cnsc]

  plt <- ggplot(bkps) +
    ylab("Watson | Crick") + xlab(NULL) +
    scale_x_continuous(breaks = scales::pretty_breaks(12), labels = format_Mb) +
    scale_y_continuous(breaks = scales::pretty_breaks(1)) +
    coord_cartesian(xlim = c(0, chr_end), ylim = c(-info_y_limit, info_y_limit)) +
    theme_classic() +
    theme(panel.spacing = unit(0, "lines"),
          strip.background = element_rect(fill = NA, colour=NA),
          plot.title = element_text(hjust = 0.5),
          rect = element_blank()
    ) +
    guides(fill = FALSE) +
    facet_grid(sample+cell ~ .) +
    geom_segment(aes(x=end, y=-info_y_limit, xend=end, yend=info_y_limit),
                 size=0.4, inherit.aes=F, linetype="dashed", col="dodgerblue") +
    ggtitle(e$chrom[1], subtitle = NULL)

  if ((!is.null(highlighted)) && (nrow(highlighted)>0)) {
    plt <- plt + geom_segment(data=highlighted, aes(x=end, y=-info_y_limit, xend=end, yend=info_y_limit), size=0.7, inherit.aes=F, linetype="solid", col="orangered")
  }

  return(plt)
}
