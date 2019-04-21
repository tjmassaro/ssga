#' Mutate a chromosome
#' 
#' Creates point mutations in a chromosome given a mutation rate
#' @param p0 Binary string of 0s and 1s
#' @param mut.rate Mutation rate
#' @return Returns p0 in which some loci have been randomly switched depending on the mutation rate
#' @examples
#' x <- '00000000'
#' m.rate <- 0.01
#' 
#' ## Example 1
#' mutate(x, m.rate)
#' [1] "00000000" #<--No mutation occurred
#' 
#' ## Example 2
#' mutate(x, m.rate)
#' [1] "00100000" #<--Mutation occurred at the 3rd locus
#' @export
mutate <- function(p0, m.rate){
  L <- nchar(p0)
  p.inv <- invert(p0)
  index <- 1:L
  p.mut <- p0
  n.mut.points <- floor(L*m.rate)
  if (n.mut.points == 0){
    coin.flip <- runif(1)
    if (coin.flip < m.rate){
      rand.gene <- sample(index, 1)
      substr(p.mut, rand.gene, rand.gene) <- substr(p.inv, rand.gene, rand.gene)
    }
  } else {
    mut.index <- sort(sample(index, n.mut.points, replace = F))
    for (i in 1:n.mut.points){
      substr(p.mut, mut.index[i], mut.index[i]) <- substr(p.inv, mut.index[i], mut.index[i])
    }
  }
  return(p.mut)
}