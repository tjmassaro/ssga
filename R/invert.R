#' Invert a chromosome
#' 
#' Takes a binary chromsome as input and flips the ordering
#' 
#' @param p0 Binary string of 0s and 1s
#' @return Returns the reverse of \code{p0}
#' @examples
#' x <- '001101'
#' invert(x)
#' [1] "101100"
#' @export
invert <- function(p0){
  L <- nchar(p0)
  p.new <- NULL
  for (i in 1:L){
    g <- as.numeric(substr(p0, i, i))
    p.new <- paste0(p.new, toString((g+1)%%2))
  }
  return(p.new)
}