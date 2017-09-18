#' splits a read count data frame based on chromosomes and returns a list containing the splited dataframes
#' @param counts dataframe containing W and C read counts in some segments
#' @author Maryam Ghareghani
#' @export
#' 

split.chromosomes <- function(counts)
{
  chr = unique(counts$chromosome)
  list.counts = list()
  for (i in 1:length(chr))
  {
    list.counts[[i]] = counts[counts$chromosome == as.character(chr)[i],]
  }
  list.counts
}

# need to be changed... We shouldn't count first pairs only