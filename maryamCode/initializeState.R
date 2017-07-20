# tested
#' returns the first binary string with n 01 and m 1s
#' @param n the number of 0s
#' @param m the number of 1s
#' @author Maryam Ghareghani
#' @export
#' 

initialState = function(n, m)
{
  paste0(str_dup("0",n), str_dup("1",m))
}