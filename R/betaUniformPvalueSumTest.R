#' A meta-analysis P-value summarising a set of hy
#' 
#' A variant of Fisher's method for P-value aggregation, but with a parametrically derived null-hypothesis that reflects the background signal in the data
#'
#' Authors in the literature (e.g. Luo et al.) have used a chi-squared test for gene set enrichment with the null-hypothesis that all genes
#' have P-values derived from a pure noise (uniform) distribution. This is typically demonstrably not true and will result in an inflated rate
#' of false positives
#' 
#' @param testPvalues A numerical vector of P-values to be tested.
#' @param betaUniformFit A beta-uniform model of the P-value distribution.
#' @param ... Unused. Included to allow for similar form across P-value sum tests.
#' 
#' @return combinedPval A single p-value combining the results of multiple hypothesis tests into one
#' @seealso fishersCombindedProbabilityTest
#' @seealso fitBetaUniformParameters
#' @seealso BioNet::fitBumModel
#' @family P-value sum tests
#' @references Luo, W., Friedman, M. S., Shedden, K., Hankenson, K. D., & Woolf, P. J. (2009). GAGE: generally applicable gene set enrichment for pathway analysis. BMC Bioinformatics
#' @importFrom stats pgamma
#' @export
setGeneric("betaUniformPvalueSumTest",
           signature = c("testPvalues", "betaUniformFit"),
           valueClass = "Pvalue",
           function(testPvalues, betaUniformFit, ...) standardGeneric("betaUniformPvalueSumTest"))

#' @describeIn betaUniformPvalueSumTest Convert to P-values on the fly
setMethod("betaUniformPvalueSumTest",
          signature = c(testPvalues = "numeric", betaUniformFit = "BetaUniformModel"),
          function(testPvalues, betaUniformFit) callGeneric(new("Pvalues", testPvalues), betaUniformFit))

#' @describeIn betaUniformPvalueSumTest Compute a P-value with the null-hypothesis being that they are drawn randomly from a mixture distribution
setMethod("betaUniformPvalueSumTest",
          signature = c(testPvalues = "Pvalues", betaUniformFit = "BetaUniformModel"),
          function(testPvalues, betaUniformFit) {
            
            nValues <- length(testPvalues)
            
            fittedBetaShape <- betaUniformFit@a
            nonUniformProportion <- betaUniformFit@lambda
            
            testStatistic <- sum(-log(testPvalues))
            
            combinedPval <- nonUniformProportion * pgamma(testStatistic, shape = nValues, rate = 1, lower.tail = FALSE) +
                        (1-nonUniformProportion) * pgamma(testStatistic, shape = nValues, rate = fittedBetaShape, lower.tail = FALSE)
            
            return( new("Pvalue", mkScalar(combinedPval)) )
          })