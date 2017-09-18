# tested
#' Compute the segment type given a cell (majority) type and the segment status and returns a string containing a number of w characters followed by a number of c characters
#' @param cellType The (majority) type of a cell that can have one these possible values: "ww","cc","wc","cw", or "?"
#' @param r the dispersion parameter for the total (W+C) number of read counts (for CN = 2)
#' @param status a string of length 4 containing {CN in hap1, inv CN in hap1, CN in hap2, inv CN in hap2} respectively
#' @author Maryam Ghareghani
#' @export
#' 

getSegType = function(cellType, status)
{
  Wstatus = NULL
  if (cellType == "?")
  {
    segType = "?"
  }
  else
  {
    for (i in 1:2)
    {
      if (substr(cellType, i, i) == "w")
        Wstatus = paste0(Wstatus, "10")
      else
        Wstatus = paste0(Wstatus, "01")
    }
    Nw = Nc = 0
    
    for (i in 1:nchar(status))
    {
      Nw = Nw + as.integer(substr(status, i, i))*as.integer(substr(Wstatus, i, i))
      Nc = Nc + as.integer(substr(status, i, i))*(1-as.integer(substr(Wstatus, i, i)))
    }
    
    segType = paste0(str_dup("w",Nw),str_dup("c",Nc))
  }
  
  segType
}