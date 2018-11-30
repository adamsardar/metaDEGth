#' A variant of Fisher's method for P-value aggregation, but with a parametrically derived null-hypothesis that reflects the background signal in the data
#'
#' Authors in the literature (e.g. Luo et al.) have used a chi-squared test for gene set enrichment with the null-hypothesis that all genes
#' have P-values derived from a pure noise (uniform) distribution. This is typically demonstrably not true and will result in an inflated rate
#' of false positives
#' 
#' @param testStat The sum of the natural logarithm of P-values of differential expression of genes in the set of interest (no need to define if geneSet is provided)
#' @param nValues The number of values in the gene set of interest
#' @param geneSet A collection of gene names that will be used to construct the test statistic (no need to define if testStat is provided)
#' @param betaUniformFit An object of class 'bum' (produced using fitBetaUniformParameters or BioNet::fitBumModel) that contains the relevant parameters of a Beta-Uniform decomposition of P-values
#'
#' @return combinedPval A single p-value combining the results of multiple hypothesis tests into one
#' @seealso fishersCombindedProbabilityTest
#' @seealso fitBetaUniformParameters
#' @seealso BioNet::fitBumModel
#' @references Luo, W., Friedman, M. S., Shedden, K., Hankenson, K. D., & Woolf, P. J. (2009). GAGE: generally applicable gene set enrichment for pathway analysis. BMC Bioinformatics
#' @importFrom stats pgamma
#' @export
mixtureGammaTest <- function(testStat = NULL, nValues = NULL, geneSet = NULL, betaUniformFit = NULL){

  #TODO check that one of either 'geneSet' or 'testStat & nValues' alongside 'betaUniformFit' has been provided (but not both)

  #TODO extract out P-Values from betaUniformFit object to produce a testStatistic
  # testStat <- sum(-log(pVals)

  nValues <- length(geneSet)
  
  #TODO validate class of betaUniformFit
  fittedBetaShape <- betaUniformFit$alpha
  nonUniformProportion <- betaUniformFit$lambda
  
  combinedPval <- nonUniformProportion * pgamma(testStat, shape = nValues, rate = 1, lower.tail = FALSE) +
                  (1-nonUniformProportion) * pgamma(testStat, shape = nValues, rate = fittedBetaShape, lower.tail = FALSE)

  return(combinedPval)
}

#' Combine multiple P-values together using Fisher's method
#' 
#' Fisher's method for combining P-values together consists of studying a test statistic of \eqn{T = \sum_i -log(P_i)} and 
#' finding the corresponding point in the cumulative distribution of an appropriately parameterised chi-squared distribution.
#'  
#' The basic assumption is that P-values are uniformly distributed under the null-hypothesis, meaning that \eqn{-log(P) \sim exp(1)}. Recall that
#' if \eqn{X \sim exp(1)} then \eqn{\sum_{k} X_{k} \sim \Gamma (k,1) \sim \chi^2 (2k) } (note that the gamma and erlang distributions are identical for integer \eqn{k}).
#' We therefore simply calculate the point in the cumulative of the gamma distribution
#' 
#' @param pVals A numerical vector of P-values
#' @return combindedP A single p-value combining the results of multiple hypothesis tests into one
#' @references \url{https://en.wikipedia.org/wiki/Fisher%27s_method}
#' @references \url{https://en.wikipedia.org/wiki/Erlang_distribution#Related_distributions}
#' @importFrom stats pgamma
#' @export
fishersCombindedProbabilityTest <- function(pVals){
  
  validatePvalues(pVals)
  
  testStatistic <- sum(-log(pVals))
  nValues <- length(pVals)
  
  #Using the pgama parameterisation so as to be consistent with the mixtureGammaTest
  combinedPval <- pgamma(testStatistic, shape = nValues, rate = 1, lower.tail = FALSE) 
  
  return(combinedPval)
}