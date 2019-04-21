#' Shotgun Stochastic Search Neighborhood Generator
#' 
#' Create the 3 neighborhoods from a shotgun stochastic search algorithm (see Hans, et al. (2007))
#' @param p0 Binary string of 0s and 1s
#' @param token A character which is either '+', '-', or 'o', depending on which neighborhood is to be generated
#' @return Returns a list containing all chromosomes making up the given neighborhood requested using \code{token}.  Let k denote the number of 1s in \code{p0}, and k' denote the number of 0s.  If \code{token == '+'}, \code{nbd} returns the set of length k' in which each chromosome is sequentially \code{p0} augmented by one of its 0s switched to a 1; if \code{token == '-'}, \code{nbd} returns the set of length k in which each chromsome is sequentially \code{p0} decreased by one of its 1s switched to a 0; if \code{token == 'o'}, returns the set of length k*k' in which each chromosome is simultaneously \code{p0} decreased by one of its 1s switched to a 0 while each of its 0s is switched to a 1.
#' @examples
#' p0 <- '00101'
#' nbd(p0, '+')
#' [1] "10101" "01101" "00111"
#' nbd(p0, '-')
#' [1] "00001" "00100"
#' nbd(p0, 'o')
#' [1] "10001" "01001" "00011" "10100" "01100" "00110"
#' nbd(p0, 'z')
#' Error: (token == "+") || (token == "-") || (token == "o") is not TRUE 
#' @export
sss.nbd <- function(p0, token) {
  stopifnot((token == '+') || (token == '-') || (token == 'o'))
  k <- nchar(p0)
  sbst <- chr2sbst(p0)
  if (p0 == paste(rep(1, k), sep = '', collapse = '') || p0 == paste(rep(0, k), sep = '', collapse = '')) {
    nbd <- p0
    return(nbd)
  } else {
    nbd <- NULL
    if (token == '+') {
      sbst.x <- setdiff(1:k, sbst)
      for (i in 1:length(sbst.x)) {
        nbd[[i]] <- sbst2chr(sort(c(sbst.x[i], sbst)), k)
      }
    } else if (token == '-') {
      for (i in 1:length(sbst)) {
        nbd[[i]] <- sbst2chr(setdiff(sbst, sbst[i]), k)
      }
    } else if (token == 'o') {
      sbst.x <- setdiff(1:k, sbst)
      counter <- 1
      for (i in 1:length(sbst)) {
        s <- setdiff(sbst, sbst[i])
        for (j in 1:length(sbst.x)) {
          nbd[[counter]] <- sbst2chr(sort(c(sbst.x[j], s)), k)
          counter <- counter + 1
        }
      }
    }
    return(nbd)
  }
}