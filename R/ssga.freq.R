#' SSGA frequency analysis
#' 
#' Relative frequencies of variables in an object of class type 'ssga'
#' @param ssga.obj An object of class 'ssga' created using \code{ssga}
#' @param B Obtain frequencies using the top \code{B} models
#' @examples
#' ## Load ICU workspace (from aplore3 package)
#' load(system.file('extdata', 'ssga_icu_ws.RData', package = 'ssga'))
#' ssga.freq(ssga.obj, B = 500)
#' 
#' ## Load Swiss workspace
#' rm(list = ls(all = TRUE))
#' load(system.file('extdata', 'ssga_swiss_ws.RData', package = 'ssga'))
#' ssga.freq(ssga.obj, B = 10)
#' @export
ssga.freq <- function(ssga.obj, B = NULL, var.names = NULL) {
  
  
  ## Front matter
  if (!is(ssga.obj, 'ssga')) {
    warning('\nNeed to pass an SSGA object.\n')
  }
  if (!is.null(var.names)) {
    label.flag <- TRUE
    max.char <- max(nchar(var.names))
  } else {
    if (!is.null(colnames(ssga.obj$subset.dat))) {
      var.names <- colnames(ssga.obj$subset.dat)
      max.char <- max(nchar(var.names))
      label.flag <- TRUE
    } else {
      label.flag <- FALSE
    }
  }
  
  
  ## Parameters
  nmod <- length(ssga.obj$subset.list)
  if (is.null(B)) {
    B <- nmod
  }
  freq.table <- matrix(0, nrow = ssga.obj$p, ncol = 1)
  
  
  ## Compute frequencies
  for (i in 1:B) {
    sbst <- chr2sbst(ssga.obj$subset.list[i])
    freq.table[sbst] <- freq.table[sbst] + 1
  }
  
  
  ## Compute relative frequency
  rel.freq <- freq.table/B
  
  
  ## Create table of results
  for (j in 1:ssga.obj$p) {
    pct <- rel.freq[j]
    if (pct == 0) {
      star <- ''
    } else if ((0 <= pct) & (pct <= 0.20)) {
      star <- '*'
    } else if ((0.20 < pct) & (pct <= 0.40)) {
      star <- '**'
    } else if ((0.40 < pct) & (pct <= 0.60)) {
      star <- '***'
    } else if ((0.6 < pct) & (pct <= 0.80)) {
      star <- '****'
    } else if (0.8 < pct) {
      star <- '*****'
    }
    if (label.flag) {
      cat(sprintf(paste0('%', as.character(max.char), 's|  %3.2f  %s\n'), var.names[j], pct, star))
    } else {
      cat(sprintf('%3d|  %3.2f  %s\n', j, pct, star))
    }
  }
  cat(sprintf('\n\n'))
}