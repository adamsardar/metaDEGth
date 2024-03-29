% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitBetaUniformModel.R
\name{fitBetaUniformMixtureDistribution}
\alias{fitBetaUniformMixtureDistribution}
\alias{fitBetaUniformMixtureDistribution,Pvalues,ScalarInteger-method}
\alias{fitBetaUniformMixtureDistribution,numeric,ANY-method}
\alias{fitBetaUniformMixtureDistribution,Pvalues,numeric-method}
\title{Fit a Beta-Uniform decomposition model to the provided p-values}
\usage{
fitBetaUniformMixtureDistribution(pValues, nStarts = 10L)

\S4method{fitBetaUniformMixtureDistribution}{Pvalues,ScalarInteger}(pValues, nStarts = 10L)

\S4method{fitBetaUniformMixtureDistribution}{numeric,ANY}(pValues, nStarts = 10L)

\S4method{fitBetaUniformMixtureDistribution}{Pvalues,numeric}(pValues, nStarts = 10L)
}
\arguments{
\item{pValues}{A numeric vector of P-values for which to perform Beta-Uniform decomposition}

\item{nStarts}{How many repeats of the fitting routine should the routine run before the best LLH is picked? The default is 10, but in some cases (when there is no 'clean' P-value distribution), it might be worth running for more iterations.}
}
\value{
A BetaUniformModel object, detailign the fit
}
\description{
Following the approach of Morris & Pounds, a distribution of P-values can be decomposed as a mixture of two beta distributions,
one with shape paramter (\eqn{\alpha}) of 1 - modelling noise under the null-hypothesis - and another with a variable (a), modelling
signal. One can view this as modeling P-values as a random process where \eqn{P ~ (1-\lambda)\beta(a,1) + \lambda\beta(1,1)}
}
\details{
Analagous to BioNet::fitBumModel.
}
\section{Functions}{
\itemize{
\item \code{fitBetaUniformMixtureDistribution,Pvalues,ScalarInteger-method}: Typed method for dispatch

\item \code{fitBetaUniformMixtureDistribution,numeric,ANY-method}: Adapter method, attempts to convert numeric values to Pvalues

\item \code{fitBetaUniformMixtureDistribution,Pvalues,numeric-method}: fitting using nStarts defined
by the user as numeric
}}

\examples{
library(ggplot2)
data(siTAZdiffex_HFL1_diffexDT)

siTAZ1_BetaUniformMod <- fitBetaUniformMixtureDistribution(siTAZdiffex_HFL1_diffexDT$siTAZ1_pValue)

siTAZ1_BetaUniformMod

plot(siTAZ1_BetaUniformMod) # Produce some diagnostic plots

betaUnifScored <- betaUniformScore(siTAZ1_BetaUniformMod, FDR = 0.05)
qplot(betaUnifScored)
}
\references{
Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
}
\seealso{
betaUniformScore BioNet::fitBumModel
}
