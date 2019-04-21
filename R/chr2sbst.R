#' Convert a chromosome to a subset
#'
#' Takes as input a binary string of 0s and 1s, and converts this to a subset for indexing
#'
#' @param p0 A string input
#' @return A vector of indices corresponding to \code{p0}
#' @examples
#' x <- '0011001'
#' chr2sbst(x)
#' [1] 3 4 7
#' @export
chr2sbst <- function(p0){
  L <- nchar(p0)
  index <- 1:L
  t.or.f <- as.logical(rep(1, L))
  for (i in 1:L){
    g <- substr(p0, i, i)
    if (identical(g, '0')){
      t.or.f[i] <- F
    }
  }
  return(index[t.or.f])
}