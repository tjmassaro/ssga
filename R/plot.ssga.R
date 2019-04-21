#' SSGA violin plot
#' 
#' Create a violin plot based on an input object of class type 'ssga'
#' @param ssga.obj An object of class 'ssga' created using \code{ssga}
#' @param B Construct the violin plot using the top \code{B} models
#' @param nsamples Total number of samples to use when creating violin plot for each coefficient
#' @param var.names Labels to use for the violin plot; by default, column names from the subsettable data are used
#' @param color.flag A TRUE/FALSE flag indicating whether to create the violin plot in color or bw
#' @examples 
#' ## Load ICU workspace (from aplore3 package)
#' load(system.file('extdata', 'ssga_icu_ws.RData', package = 'ssga'))
#' plot(ssga.obj, B = 100, nsamples = 1000, color.flag = TRUE)
#' 
#' ## Load Swiss workspace
#' rm(list = ls(all = TRUE))
#' load(system.file('extdata', 'ssga_swiss_ws.RData', package = 'ssga'))
#' plot(ssga.obj)
#' @export
plot.ssga <- function(ssga.obj, B = NULL, nsamples = NULL, var.names = NULL, color.flag = FALSE) {
  
  
  ## Front matter
  if (!is(ssga.obj, 'ssga')) {
    warning('\nNeed to pass an SSGA object.\n')
  }
  p <- ssga.obj$p
  fp <- ssga.obj$fp
  X <- ssga.obj$subset.dat
  y <- ssga.obj$response
  glm.fam <- ssga.obj$glm.fam
  if (fp == 0) {
    forced.vars <- NULL
  } else {
    forced.vars <- ssga.obj$forced.dat
  }
  
  
  ## Set the number of models
  if (is.null(B)) B <- 10
  
  
  ## Set the number of samples
  if (is.null(nsamples)) nsamples <- 100
  
  
  ## Assign variable labels for the violin plot
  if (is.null(var.names) & is.null(colnames(ssga.obj$subset.dat))) {
    var.names <- NA
    for (j in 1:(p + fp)) {
      var.names[j] <- paste0('V', j)
    }
  } else if (is.null(var.names) & !is.null(colnames(ssga.obj$subset.dat))) {
    if (fp == 0) {
      var.names <- colnames(ssga.obj$subset.dat)
    }
    else {
      var.names <- c(rep('NA', fp), colnames(ssga.obj$subset.dat))
    }
  }
  
  
  ## Violin plot data
  beta.list <- list()
  for (i in 1:B) {
    chr <- ssga.obj$subset.list[i]
    sbst <- chr2sbst(chr)
    mod <- glm(y ~ cbind(forced.vars, X[ , sbst]), family = glm.fam)
    b <- mod$coefficients
    c <- summary(mod)$cov.scaled
    beta.sample <- rmvnorm(nsamples, mean = b, sigma = c)
    beta.list[[i]] <- beta.sample[ , 2:ncol(beta.sample)]
  }
  s <- ssga.obj$score.list[1:B]/sd(ssga.obj$score.list[1:B])
  w <- exp(min(s) - s)/sum(exp(min(s) - s))
  sample.order <- sample(nsamples)
  model.order <- sample(1:B, size = nsamples, prob = w, replace = TRUE)
  vs <- data.frame(matrix(0, nrow = nsamples, ncol = fp + p))
  for (i in 1:nsamples) {
    chr <- ssga.obj$subset.list[model.order[i]]
    sbst <- chr2sbst(chr)
    b <- as.matrix(beta.list[[model.order[i]]], nrow = nsamples, ncol = length(sbst))
    if (fp == 0) {
      vs[i, sbst] <- b[sample.order[i], ]
    } else {
      vs[i, c(1:fp, fp + sbst)] <- b[sample.order[i], ]
    }
  }
  names(vs) <- var.names
  
  
  ## Restructure data, create the violin plot
  if (color.flag) {
    vs.dat <- stack(vs)
    names(vs.dat) <- c('dat', 'label')
    p <- ggplot(vs.dat, aes(x = label, y = dat, fill = label)) + 
      geom_violin(trim = TRUE, scale = 'width') + 
      coord_flip() + 
      labs(x = '', y = '') + 
      stat_summary(fun.y = median, geom = 'point', size = 2) +
      theme(legend.position = 'none')
    print(p)
  } else {
    vs.dat <- stack(vs)
    names(vs.dat) <- c('dat', 'label')
    p <- ggplot(vs.dat, aes(x = label, y = dat)) + 
      geom_violin(trim = TRUE, scale = 'width') + 
      coord_flip() + 
      labs(x = '', y = '') + 
      stat_summary(fun.y = median, geom = 'point', size = 2) + 
      theme(legend.position = 'none')
    print(p)
  }
  
  
  ## Add bounds and weighted means to output
  output <- list(beta.est = vs)
  
  
  ## Return output
  return(output)
}