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
#' @references Stewart T, Strijbosch LWG, Moors H, van Batenburg P. A Simple Approximation to the Convolution of Gamma Distributions. 2007
#' @importFrom stats pgamma
#' @include Pvalues_S4Class.R
#' @export
setGeneric("betaUniformPvalueSumTest",
           signature = c("testPvalues", "betaUniformFit"),
           valueClass = "Pvalue",
           function(testPvalues, betaUniformFit, na.rm = TRUE) standardGeneric("betaUniformPvalueSumTest"))

#' @describeIn betaUniformPvalueSumTest Convert to P-values on the fly
setMethod("betaUniformPvalueSumTest",
          signature = c(testPvalues = "numeric", betaUniformFit = "BetaUniformModel"),
          function(testPvalues, betaUniformFit, na.rm = TRUE) callGeneric(new("Pvalues", testPvalues), betaUniformFit, na.rm))

#' @describeIn betaUniformPvalueSumTest Compute a P-value with the null-hypothesis being that they are drawn randomly from a mixture distribution
setMethod("betaUniformPvalueSumTest",
          signature = c(testPvalues = "Pvalues", betaUniformFit = "BetaUniformModel"),
          function(testPvalues, betaUniformFit, na.rm = TRUE) {
            
            if(na.rm){testPvalues <- na.omit(testPvalues)}
            
            nValues <- length(testPvalues)
            
            fittedBetaShape <- betaUniformFit@a
            uniformProportion <- betaUniformFit@lambda
            
            testStatistic <- sum(-log(testPvalues))
            
            mu <- randomSetSize*(uniformProportion +  (1-uniformProportion)*scaleParam)
            sigmaSq <-  randomSetSize*(uniformProportion +  (1-uniformProportion)*scaleParam^2)
            
            # Convolution of two gamma functions, taken from Stewart et al.
            combinedPval <- pgamma(testStatistic, shape = (mu^2)/sigmaSq, scale = sigmaSq/mu, lower.tail = F)
            
            return( new("Pvalue", mkScalar(combinedPval)) )
          })