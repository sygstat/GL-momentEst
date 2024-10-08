% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glme.gev.R
\name{glme.gev}
\alias{glme.gev}
\title{Generalized L-moments estimation for generalized extreme value distribution}
\usage{
glme.gev(xdat, ntry = 10, pen = "beta")
}
\arguments{
\item{xdat}{A numeric vector of data to be fitted.}

\item{ntry}{Number of attempts for parameter estimation. Higher values increase the chance of finding a more accurate estimate by trying different initial conditions.}

\item{pen}{Penalization method used in estimation. Options are 'beta' and 'norm'.}
}
\value{
The glme.gev function returns a list containing the following elements:
\itemize{
 \item glme - The estimated parameters of the Generalized Extreme Value distribution.
 \item lme - The L-moment estimates of the parameters.
 \item covinv - The inverse of the covariance matrix of the L-moments.
 \item lcovdet - The log determinant of the covariance matrix.
 \item nllh.pref - The negative log-likelihood of the preferred solution.
 \item pen - The penalization method used ('beta' or 'norm').
}
}
\description{
This function estimates the Generalized L-moments of Generalized Extreme Value distribution using an updated algorithm.
}
\details{
The function uses an optimization approach to estimate the parameters of the Generalized Extreme Value distribution.
It implements a penalized likelihood method, where the penalization can be either 'beta' or 'norm'.
The function handles potential issues with covariance matrix calculation by using a bootstrap approach when necessary.
}
\author{
[Jeong-Soo Park]
}
