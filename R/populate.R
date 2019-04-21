#' Populate
#' 
#' Create an initial population of chromosomes
#' @param pop.size The number of chromosomes to generate
#' @param n.loci The number of loci on each chromosome, corresponding to the number of features
#' @return Returns a list of \code{pop.size} randomly generated \code{n.loci}-length chromosomes
#' @examples
#' pop.size <- 2
#' n.loci <- 10
#' populate(pop.size, n.loci)
#' [1] "0101101000" "1100010100"
#' @export
populate <- function(pop.size, n.genes){
  F0 <- NULL
  for (i in 1:pop.size){
    chr <- sample(0:1, n.genes, replace = T)
    F0[[i]] <- paste(chr, sep = '', collapse = '')
  }
  return(F0)
}