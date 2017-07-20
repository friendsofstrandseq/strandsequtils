#' filters the list of dataframes based on the list od indexs
#' @param list.df list of dataframes
#' @param list.idx list of indexes for subsetting rows of the dataframe
#' @author Maryam Ghareghani
#' @export
#' 

filt <- function(list.df, list.idx)
{
  new.list.df = list()
  for (i in 1:length(list.df))
  {
    new.list.df[[i]] = list.df[[i]][list.idx[[i]],]
  }
  new.list.df
}

# need to be changed... We shouldn't count first pairs only
