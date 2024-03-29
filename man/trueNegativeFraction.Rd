% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BetaUniformModel_S4Class.R
\name{trueNegativeFraction}
\alias{trueNegativeFraction}
\alias{trueNegativeFraction,BetaUniformModel,ScalarNumeric-method}
\alias{trueNegativeFraction,ANY,missing-method}
\alias{trueNegativeFraction,ANY,numeric-method}
\title{True negative rate from M&P: D table 1 and equation 5}
\usage{
trueNegativeFraction(obj, pValueThreshold)

\S4method{trueNegativeFraction}{BetaUniformModel,ScalarNumeric}(obj, pValueThreshold)

\S4method{trueNegativeFraction}{ANY,missing}(obj, pValueThreshold)

\S4method{trueNegativeFraction}{ANY,numeric}(obj, pValueThreshold)
}
\arguments{
\item{obj}{data object}

\item{pValueThreshold}{p-value threshold}
}
\description{
True negative rate from M&P: D table 1 and equation 5
}
\section{Functions}{
\itemize{
\item \code{trueNegativeFraction,BetaUniformModel,ScalarNumeric-method}: method for dispatch

\item \code{trueNegativeFraction,ANY,missing-method}: compute true negative rate using a
threshold of 0.05

\item \code{trueNegativeFraction,ANY,numeric-method}: compute true negative rate using
numeric threshold defined by the user
}}

