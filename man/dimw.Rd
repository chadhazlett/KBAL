% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{dimw}
\alias{dimw}
\title{Difference in Means and Difference in Weighted Means}
\usage{
dimw(X, w, target)
}
\arguments{
\item{X}{matrix of data where rows are observations and columns are covariates.}

\item{w}{numeric vector of weights for each observation.}

\item{target}{numeric vector of length equal to the total number of units where population/treated units take a value of 1 and sample/control units take a value of 0.}
}
\value{
\item{dim}{the simple, unweighted difference in means.}
\item{dimw}{the weighted difference in means.}
}
\description{
Calculates the simple difference in means or weighted difference in means between the control or sample population and the treated or target population.
}
\examples{
\donttest{
#let's say we want to get the unweighted DIM and the weighted DIM using weights from the kbal
#function with the lalonde data:
#load and clean data a bit
set.seed(123)
data("lalonde")
# Select a random subset of 500 rows
lalonde_sample <- sample(1:nrow(lalonde), 500, replace = FALSE)
lalonde <- lalonde[lalonde_sample, ]

xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")

#get the kbal weights
kbalout= kbal(allx=lalonde[,xvars],
               sampledinpop=FALSE,
               treatment=lalonde$nsw)
 #now use dimw to get the DIMs
 dimw(X = lalonde[,xvars], w = kbalout$w, target = lalonde$nsw)}
}
