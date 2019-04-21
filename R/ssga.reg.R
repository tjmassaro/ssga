#' Shotgun Stochastic Genetic Algorithm
#'
#' Simultaneously perform variable subet selection and model averaging with the shotgun stochastic genetic algorithm regression
#' @param X A dataset whose rows represent a vector of p measurements across a set of independent variables
#' @param y A set of outcomes
#' @param model.type Specify the type of regression model ('linear', 'logistic', 'poisson')
#' @param fitness.fun Specify the fitness function ('AIC', 'BIC', 'CUSTOM')
#' @param forced.vars A set of 1 or more columns with the same number of rows as \code{X}, representing variables that are to be included in every regression model
#' @param print.flag When \code{print.flag} is set to TRUE, the current trial, iteration, best model fitness score, and chromosome corresponding to the best model will be shown in the console
#' @param penalty If \code{fitness.fun = 'CUSTOM'}, the fitness function takes the form -2*logLike + \code{penalty}*nparam, where logLike refers to the maximized log-likelihood of the regression model and nparam is the total number of free parameters in the model; note: AIC = -2*logLike + 2*nparam
#' @examples
#' ## Load ICU workspace for logistic regression (data were obtained from aplore3 package)
#' load(system.file('extdata', 'ssga_icu_ws.RData', package = 'ssga'))
#' no.trials <- 5
#' model.type <- 'logistic'
#' fitness.fun <- 'AIC'
#' ssga.obj <- ssga.reg(X, y,
#'                      no.trials = no.trials,
#'                      model.type = model.type,
#'                      fitness.fun = fitness.fun,
#'                      print.flag = TRUE)
#' summary(ssga.obj)
#'
#'
#' ## Poisson regression with the ICU data
#' model.type <- 'poisson'
#' y <- X[ , 'age']
#' X <- X[ , !colnames(X)%in%'age']
#' ssga.obj <- ssga.reg(X, y,
#'                      no.trials = no.trials,
#'                      model.type = model.type,
#'                      fitness.fun = fitness.fun)
#' summary(ssga.obj)
#'
#'
#' ## Load Swiss workspace for linear regression
#' rm(list = ls(all = TRUE))
#' load(system.file('extdata', 'ssga_swiss_ws.RData', package = 'ssga'))
#' no.trials <- 5
#' model.type <- 'linear'
#' fitness.fun <- 'AIC'
#' ssga.obj <- ssga.reg(X, y,
#'                      no.trials = no.trials,
#'                      model.type = model.type,
#'                      fitness.fun = fitness.fun,
#'                      print.flag = TRUE)
#' summary(ssga.obj)
#' @return Returns an object of class 'ssga'
#' @export
ssga.reg <- function(X, y, no.trials, model.type = 'linear', fitness.fun = 'AIC',
                 forced.vars = NULL, print.flag = TRUE, penalty = NULL) {


  ## Start the timer
  start.time <- Sys.time()


  ## Get subsettable data dimensions
  p <- dim(X)[2]
  n <- dim(X)[1]
  if (p < 4) stop('SSGA reqiures at least 4 subsettable variables.')


  ## Front matter
  if (is.null(forced.vars)) {
    fp <- 0
  } else {
    fp <- dim(forced.vars)[2]
  }


  ## Determine custom parameter
  if (fitness.fun == 'CUSTOM' & is.null(penalty)) {
    stop('User must provide a penalty for a custom fitness function.')
  }
  if (is.null(penalty)) {
    penalty <- 0
  }
  if (penalty < 0) {
    stop('The penalty cannot be a negative number.')
  }


  ## Create containers for fitness function and subset
  subset.list <- NULL
  score.list <- NULL


  ## Create parameters for algorithm
  N <- 10
  maxits <- 1E3
  thresh <- p
  m.rate <- 0.01


  ## Assign GLM family
  if (model.type == 'linear') {
    glm.fam <- as.function(alist(a =  , b = 'identity', a(link = b)))(gaussian)
  } else if (model.type == 'logistic') {
    glm.fam <- as.function(alist(a = , b = 'logit', a(link = b)))(binomial)
  } else if (model.type == 'poisson') {
    glm.fam <- as.function(alist(a = , b = 'log', a(link = b)))(poisson)
  }


  ## Run SSGA
  for (t in 1:no.trials) {


    ## Print trial
    if (print.flag) cat(sprintf('\n\n\nTRIAL %d\n\n', t))


    ## Local parameters
    nits <- 1
    F0 <- populate(N, p)
    subset <- NULL
    score <- NULL
    old.f1 <- 1E8
    f1 <- 0
    f2 <- 0
    phase.flag <- 0
    phase1.best <- 1E8
    phase2.best <- 1E8


    ## Run
    while (phase.flag == 0) {


      ## PHASE I
      phase1.fitness <- NULL
      phase1.aic <- NULL
      phase1.bic <- NULL
      phase1.custom <- NULL

      counter <- 1
      while (counter < thresh) {
        for (i in 1:N) {
          sbst <- chr2sbst(F0[i])
          if (length(sbst) == 0) {
            F0[i] <- paste(rep(1, p), sep = '', collapse = '')
            sbst <- chr2sbst(F0[i])
          }
          model <- glm(y ~ cbind(forced.vars, X[ , sbst]), family = glm.fam)
          nparam <- length(sbst) + fp
          phase1.aic[i] <- model$deviance + 2*nparam + 2
          phase1.bic[i] <- model$deviance + nparam*log(n)
          phase1.custom[i] <- model$deviance + penalty*nparam
        }
        if (fitness.fun == 'AIC') {
          phase1.fitness <- phase1.aic
        } else if (fitness.fun == 'BIC') {
          phase1.fitness <- phase1.bic
        } else if (fitness.fun == 'CUSTOM') {
          phase1.fitness <- phase1.custom
        }
        score <- c(score, phase1.fitness[!is.element(F0, subset)])
        subset <- c(subset, F0[!is.element(F0, subset)])
        F1 <- next.gen(F0, phase1.fitness, m.rate/counter)
        f1 <- min(phase1.fitness)
        if (f1 == old.f1) {
          counter <- counter + 1
        } else {
          counter <- 1
          old.f1 <- f1
        }
        if (print.flag) cat(sprintf('%3d%13.2f     %30s\n', nits, f1, F1[1]))
        nits <- nits + 1
        F0 <- F1
      }
      phase1.best <- f1


      ## Compare phases
      if (phase1.best == phase2.best) {
        phase.flag <- 1
      }


      ## PHASE II
      m <- F0[1]
      old.f2 <- f1
      flag <- 0
      while (flag == 0) {
        nbd <- c(sss.nbd(m, '+'), sss.nbd(m, '-'), sss.nbd(m, 'o'))
        nbd <- c(m, nbd[!is.element(nbd, subset)])
        new.N <- length(nbd)
        phase2.aic <- NULL
        phase2.bic <- NULL
        phase2.custom <- NULL

        for (i in 1:new.N) {
          sbst <- chr2sbst(nbd[i])
          model <- glm(y ~ cbind(forced.vars, X[ , sbst]), family = glm.fam)
          cov.mat <- summary(model)$cov.scaled
          nparam <- length(sbst) + fp
          phase2.aic[i] <- model$deviance + 2*nparam + 2
          phase2.bic[i] <- model$deviance + nparam*log(n)
          phase2.custom[i] <- model$deviance + penalty*nparam
        }
        if (fitness.fun == 'AIC') {
          new.fit <- phase2.aic
        } else if (fitness.fun == 'BIC') {
          new.fit <- phase2.bic
        } else if (fitness.fun == 'CUSTOM') {
          new.fit <- phase2.custom
        }
        f2 <- min(new.fit)
        if (f2 < old.f2) {
          best <- which(new.fit == min(new.fit))
          m <- nbd[best]
          score <- c(score, new.fit[!is.element(nbd, subset)])
          subset <- c(subset, nbd[!is.element(nbd, subset)])
          if (print.flag) cat(sprintf('%3d%13.2f     %30s\n', nits, f2, m))
          nits <- nits + 1
          old.f2 <- f2
        } else {
          flag <- 1
        }
        F0[1] <- m
      }
      phase2.best <- f2


      ## Compare phases
      if (phase1.best == phase2.best) {
        phase.flag <- 1
      }
    }

    subset <- subset[!duplicated(subset)]
    score <- score[!duplicated(score)]
    score.list <- c(score.list, score[!is.element(subset, subset.list)])
    subset.list <- c(subset.list, subset[!is.element(subset, subset.list)])
  }


  ## End time
  end.time <- Sys.time()


  ## Sort lists before exporting
  sort.list <- unlist(as.matrix(sort.int(score.list, index.return = TRUE))[2])
  output <- list(subset.list = subset.list[sort.list],
                 score.list = score.list[sort.list],
                 glm.fam = glm.fam,
                 model.type = model.type,
                 fitness.fun = fitness.fun,
                 p = p,
                 fp = fp,
                 subset.dat = X,
                 response = y,
                 time.elapsed = end.time - start.time,
                 no.trials = no.trials)
  if (!is.null(forced.vars)) output$forced.dat = forced.vars
  if (!is.null(penalty)) output$penalty = penalty
  class(output) <- append(class(output), 'ssga')


  ## Output results
  cat(sprintf('\n\n\n'))
  return(output)
}
