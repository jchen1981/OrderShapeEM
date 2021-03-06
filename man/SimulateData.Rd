% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OrderShapeEM.R
\name{SimulateData}
\alias{SimulateData}
\title{Simulate p-values and the auxiliary covariate under various scenarios.}
\usage{
SimulateData(
  prior.strength = c("Weak", "Moderate", "Strong"),
  feature.no = 10000,
  sig.dist = c("Normal", "Gamma"),
  sig.density = c("None", "Lower", "Low", "Medium", "High"),
  sig.strength = c("Weak", "Moderate", "Strong")
)
}
\arguments{
\item{prior.strength}{a character string from \code{'Weak', 'Moderate', 'Strong'} indicating the covariate informativeness.}

\item{feature.no}{an integer, the number of features to be simulated.}

\item{sig.dist}{a character string from \code{'Normal', 'Gamma'} indicating the distribution of the z-value under the alternative.}

\item{sig.density}{a character string from \code{'None', 'Lower', 'Low', 'Medium', 'High'} indicating the level of the signal density.}

\item{sig.strength}{a character string from \code{'Weak', 'Moderate', 'Strong'} indicating the level of the signal strength.}
}
\value{
A list with the elements
\item{pvalue}{a numeric vector of p-values.}
\item{prior}{a vector of covariate values reflecting the order of the prior null probabilities.}
\item{truth}{a vector of logical values indicating H0 (=0) or H1 (=1).}
}
\description{
The function simulates p-values and the auxiliary covariate under different signal structures (density and strength) and covariate informativeness.
}
\references{
Hongyuan Cao, Jun Chen, Xianyang Zhang. Optimal false discovery rate control for large-scale multiple testing with
 auxiliary information. Submitted.
}
\author{
Jun Chen
}
\keyword{Simulation}
