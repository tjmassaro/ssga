#' Cross two chromosomes
#' 
#' Simulates a random genetic crossover event between 2 binary chromosomes
#' 
#' @param p0 First chromosome
#' @param p1 Second chromosome
#' @return Returns a child string in which a section of p0 has been replaced by the corresponding section in p1
#' @examples
#' ## Example 1
#' p0 <- '0000000'
#' p1 <- '1111111'
#' crossover(p0, p1)
#' [1] "0001110" #<--This corresponds to a crossover occuring from the 4th through 6th loci
#' 
#' ## Example 2
#' p0 <- '0101'
#' p1 <- '00001101010'
#' crossover(p0, p1)
#' Error: L0 == L1 is not TRUE
#' @export
crossover <- function(p0, p1){
  L0 <- nchar(p0)
  L1 <- nchar(p1)
  stopifnot(L0 == L1)
  index <- 1:L0
  child <- NULL
  start <- 1
  which_parent <- 0
  n.cross.points <- 2
  cross.points <- sort(sample(index, n.cross.points, replace = F))
  for (i in 1:n.cross.points){
    which_parent <- which_parent%%2
    parent <- get(toString(paste0('p', toString(which_parent))))
    child <- paste0(child, substr(parent, start, cross.points[i]))
    start <- cross.points[i] + 1
    which_parent <- which_parent + 1
  }
  which_parent <- which_parent%%2
  parent <- get(toString(paste0('p', toString(which_parent))))
  child <- paste0(child, substr(parent, start, L0))
  return(child)
}