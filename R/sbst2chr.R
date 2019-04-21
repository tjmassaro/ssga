#' Convert a subset to a chromosome
#' 
#' Takes a vector of indices and converts it to a chromosome
#' @param sbst A vector of indices
#' @param k The desired length of the chromosome
#' @return Returns a binary string of length \code{k} corresponding to which features are present or absent in the model
#' @examples 
#' x <- c(1, 2, 5, 7)
#' k <- 9
#' sbst2chr(x, k)
#' [1] "110010100"
#' @export
sbst2chr <- function(sbst, k) {
  L <- is.element(1:k, sbst)
  return(paste(sprintf('%d', L), collapse = ''))
}