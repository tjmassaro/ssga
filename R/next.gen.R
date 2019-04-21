#' Create the next generation of chromosomes
#' 
#' Generate a subsequent population of chromosomes given a current population and their fitness values
#' @param F0 A list containing a set of binary strings
#' @param fitness A vector with the same length as \code{F0} containing the fitness values
#' @return Returns a list with the same dimensions as \code{F0} whose entries correspond to the chromosomes of the next generation of models
#' @examples
#' F0 <- c('00101', '11100', '01000')
#' fitness <- c(220, 180, 200)
#' m.rate <- 0.01
#' next.gen(F0, fitness, m.rate)
#' [1] "11100" "10000" "11101" #<--NB: the 'elite' model (i.e., model with the best fitness) reappears in the first position
#' @export
next.gen <- function(F0, fitness, m.rate) {
  pop.size <- length(F0)
  stopifnot(pop.size == length(fitness))
  n.var <- nchar(F0[[1]])
  F1 <- NULL
  
  
  ## Elitism
  best.chr <- F0[which(fitness == min(fitness))][1]
  F1[[1]] <- best.chr
  
  
  ## Create SSS neighborhoods
  g <- c(best.chr, sss.nbd(best.chr, '+'), sss.nbd(best.chr, '-'), sss.nbd(best.chr, 'o'))
  
  
  ## Create next generation
  n.nbd <- length(g)
  for (i in 2:pop.size) {
    rand.ints <- sample(1:n.nbd, 2, replace = FALSE)
    p1 <- g[rand.ints[1]]
    p2 <- g[rand.ints[2]]
    child <- crossover(p1, p2)
    F1[[i]] <- mutate(child, m.rate)
  }
  return(F1)
}