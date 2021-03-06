% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.ssga.R
\name{summary.ssga}
\alias{summary.ssga}
\title{SSGA summary}
\usage{
\method{summary}{ssga}(ssga.obj)
}
\arguments{
\item{ssga.obj}{An object of class 'ssga' created using \code{ssga}}
}
\description{
Summarize an object of class type 'ssga'
}
\examples{
## Load ICU workspace (from aplore3 package)
load(system.file('extdata', 'ssga_icu_ws.RData', package = 'ssga'))
summary(ssga.obj)

## Load Swiss workspace
rm(list = ls(all = TRUE))
load(system.file('extdata', 'ssga_swiss_ws.RData', package = 'ssga'))
summary(ssga.obj)
}
