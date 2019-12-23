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
#' @references \url{https://en.wikipedia.org/wiki/Phase-type_distribution}
#' @importFrom actuar pphtype
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
            if(length(testPvalues) == 0){stop("No non-NA P-values passed in to routine!")}
            if(any(testPvalues < .Machine$double.xmin)){warning("One or more test P-values are less than machine precision. This could easily result in inaccuracies.")}
            
            nValues <- length(testPvalues)
            fittedBetaShape <- betaUniformFit@a
            uniformProportion <- betaUniformFit@lambda
            
            testStatistic <- sum(-log(testPvalues))
            
            transitionMatrix <- constructHyperexponentialSumTransitionMatrix(nValues, uniformProportion, fittedBetaShape)
            
            initialProb <- Matrix(0,ncol = ncol(transitionMatrix), nrow = 1)
            initialProb[1:2] <- c(uniformProportion, 1-uniformProportion)
            
            # The inverse CDF of the phase-type distribution is sum( alpha * exp(x*S) ), where S is the transition matrix and alpha the initial probability
            phaseRateP <- sum(initialProb %*% Matrix::expm(testStatistic*transitionMatrix))
          
            return( new("Pvalue", mkScalar(phaseRateP)) )
          })

#' Model a sum of log-transformed beta-uniform mixture draws as a sum of hyperexponential distributions: a phase-type distribution
#' 
#' The distribution of a log-trasnformed beta-uniform draw is a hyperexponential distribution. The distribution of the sum of multiple
#' such hyperexponentials does nto have a name, but can be treated as a phase-type distribution with appropriate transition matrix.
#' 
#' @references \url{https://en.wikipedia.org/wiki/Phase-type_distribution}
#' @import Matrix
constructHyperexponentialSumTransitionMatrix <- function(nValues, uniformProportion, fittedBetaShape){

  recurrentStateMatrixTile <- Matrix(c(-1,0,0,-fittedBetaShape), ncol = 2, nrow = 2)
  transitionStateMatrixTile <- Matrix(c(uniformProportion, uniformProportion*fittedBetaShape,
                                       (1-uniformProportion),(1-uniformProportion)*fittedBetaShape), nrow = 2, ncol = 2)
  
  hyperexponentialSumTransitionMatrix <- kronecker(Matrix::diag(nrow = nValues, ncol = nValues), recurrentStateMatrixTile)
  
  # If the number of values is 1, then we have a hyperexponeitial distribution and this transition matrix is suitable
  
  # For nValues = 2, we can just insert the tile describing transition from state 1 to 2 only
  if(nValues == 2){  hyperexponentialSumTransitionMatrix[1:2,3:4] <- transitionStateMatrixTile }
  
  # for the general case, 
  if(nValues > 2){
    
    hyperexponentialSumTransitionMatrix <- hyperexponentialSumTransitionMatrix +
      # Create an off diagonal matrix and add the transition tiles in
      kronecker({M <- Matrix(0, nrow = nValues, ncol = nValues); diag(M[,-1]) <- 1; M}, transitionStateMatrixTile) }
  
  return(hyperexponentialSumTransitionMatrix)
}