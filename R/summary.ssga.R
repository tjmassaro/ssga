#' SSGA summary
#'
#' Summarize an object of class type 'ssga'
#' @param ssga.obj An object of class 'ssga' created using \code{ssga}
#' @examples
#' ## Load ICU workspace (from aplore3 package)
#' load(system.file('extdata', 'ssga_icu_ws.RData', package = 'ssga'))
#' summary(ssga.obj)
#'
#' ## Load Swiss workspace
#' rm(list = ls(all = TRUE))
#' load(system.file('extdata', 'ssga_swiss_ws.RData', package = 'ssga'))
#' summary(ssga.obj)
#' @export
summary.ssga <- function(ssga.obj) {


  ## Front matter
  if (!is(ssga.obj, 'ssga')) {
    warning('\nNeed to pass an SSGA object.\n')
  }
  cat(sprintf('\n\n\nSUMMARY OF SSGA RESULTS:\n________________________\n'))
  cat(sprintf('\n\nTime elapsed (sec): %d', ceiling(as.numeric(ssga.obj$time.elapsed))))
  cat(sprintf('\nNo. unique models visited: %d', length(ssga.obj$subset.list)))
  cat(sprintf('\nNo. trials: %d', ssga.obj$no.trials))
  cat(sprintf('\nNo. obs.: %d', dim(ssga.obj$subset.dat)[1]))
  cat(sprintf('\nNo. subsettable variables: %d', ssga.obj$p))
  cat(sprintf('\nNo. forced variables: %d', ssga.obj$fp))
  cat(sprintf('\nModel type: %s', ssga.obj$model.type))
  cat(sprintf('\nFitness function: %s ', ssga.obj$fitness.fun))
  if (ssga.obj$fitness.fun == 'CUSTOM') cat(sprintf('(-2*logLike + %.2f*nparam)', ssga.obj$penalty))


  ## Display the 10 best models as chromosomes
  cat(sprintf('\n\n\nTop 10 models (chromosomes):\n'))
  for (i in 1:10) {
    cat(sprintf('%2d   %f   %s\n', i, ssga.obj$score.list[i], ssga.obj$subset.list[i]))
  }


  ## Create variable names
  if (is.null(colnames(ssga.obj$subset.dat))) {
    var.names <- NA
    for (j in 1:p) {
      var.names[j] <- paste0('V', j)
    }
  } else {
    var.names <- colnames(ssga.obj$subset.dat)
  }


  ## Display the 10 best models as vectors
  if (is.null(colnames(ssga.obj$subset.dat))) {
    cat(sprintf('\n\nTop 10 models (subsets):\n'))
    for (i in 1:10) {
      cat(sprintf('%2d   %f   %s\n', i, ssga.obj$score.list[i],
                  paste(chr2sbst(ssga.obj$subset.list[i]), sep = ' ', collapse = ' ')))
    }
  } else {
    cat(sprintf('\n\nTop 10 models (subsets):\n'))
    for (i in 1:10) {
      cat(sprintf('%2d   %f   %s\n', i, ssga.obj$score.list[i],
                  paste(var.names[chr2sbst(ssga.obj$subset.list[i])],
                        sep = ' ', collapse = ', ')))
    }
  }


  ## Frequency analysis of top 10 models
  cat(sprintf('\n\nRelative frequencies of variables\nappearing in the top 10 models (0-5 stars):\n'))
  ssga.freq(ssga.obj, B = 10)


  ## Display best model results

}
