#' Compute the adjusted log likelihood ratio score for beta-uniform model of P-values
#' 
#' Relative to some false discovery threshold, compute the relative ratio of probabilities and then take the natural logarithm.
#' 
#' This function is somewhat analagous to BioNet::scoreNodes
#' 
#' @param betaUniformFit The result of fitBetaUniformParameters method
#' @param FDR The tolerable false discovery rate. See Morris & Pounds (2003)
#' @param pVals Vector of P-values to score
#' 
#' @return betaUniformScores A vector of P-Value scores
#' @export
#' @importFrom stats pbeta
#' @include Pvalues-S4class.R
#' @seealso fitBetaUniformParameters BioNet::scoreNodes
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
#' @references Beisser D, Klau GW, Dandekar T, Müller T, Dittrich MT. BioNet: An R-Package for the functional analysis of biological networks. Bioinformatics. 2010
setGeneric("betaUniformScore",
           valueClass = "numeric",
           function(x, betaUniformFit = NULL, FDR = 5E-2) { standardGeneric("betaUniformScore")},
           signature = c("x", "betaUniformFit", "FDR") )


#' @describeIn betaUniformScore Score P-values against an explicit beta-uniform model
setMethod("betaUniformScore",
          c(x="Pvalues", betaUniformFit = "BetaUniformModel", FDR = "ScalarNumeric"),
          function(x, betaUniformFit, FDR) {
            
            FDRthreshold <- pValueThresholdAtConfidence(betaUniformFit, FDR)
            
            return( (betaUniformFit@a - 1)*(log(x) - log(FDRthreshold) ) )
          })



setOldClass("igraph")
#' @describeIn betaUniformScore Score nodes in a graph; there must be a 'pValue' attribute
#' @import igraph
setMethod("betaUniformScore",
          c(x="igraph", betaUniformFit = "ANY", FDR = "ANY"),
          function(x, betaUniformFit, FDR = 5E-2) { callGeneric(  V(x)$pValue, names=V(x)$pValue , betaUniformFit, FDR) })

#' @describeIn betaUniformScore Score numeric values, assumed to be P-values
setMethod("betaUniformScore",
          c(x="numeric", betaUniformFit = "ANY", FDR = "ANY"),
          function(x, betaUniformFit, FDR) { callGeneric( new("Pvalues", .Data = x) , betaUniformFit, FDR) })

#' @importFrom Biobase mkScalar
setMethod("betaUniformScore",
          c(x="ANY", betaUniformFit = "ANY", FDR = "numeric"),
          function(x, betaUniformFit, FDR) {  callGeneric(x, betaUniformFit, mkScalar(FDR))  })

#' @describeIn betaUniformScore Score P-values. Beta-uniform model fit on the fly.
#' @include Pvalues-S4class.R
setMethod("betaUniformScore",
          c(x="Pvalues", betaUniformFit = "missing", FDR = "ANY"),
          function(x, betaUniformFit, FDR = 5E-2) {
            
            callGeneric(x, fitBetaUniformMixtureDistribution(x.Data), FDR)
          })


setOldClass("bum") # from BioNet
#' @describeIn betaUniformScore Score P-values against an explicit beta-uniform model
#' @include Pvalues-S4class.R
setMethod("betaUniformScore",
          c(x="ANY", betaUniformFit = "bum", FDR = "ANY"),
          function(x, betaUniformFit, FDR){ callGeneric(x, convertFromBUM(betaUniformFit), FDR) })


