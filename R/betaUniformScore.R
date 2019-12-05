#' Fo each P-value, compute the log likelihood ratio of originating from signal or noise component
#' 
#' Relative to some false discovery threshold, compute the relative ratio of probabilities and then take the natural logarithm.
#' 
#' This function is somewhat analagous to BioNet::scoreNodes
#' 
#' @param betaUniformFit (optional) A model . Set to NULL to autofit on the fly.
#' @param FDR The tolerable fraction of false positives within the set of positive scoring values. See Morris & Pounds (2003)
#' @param pVals Vector of P-values to score
#' 
#' @return betaUniformScores A vector of P-Value scores
#' @export
#' @import igraph
#' @importFrom stats pbeta
#' @importFrom Biobase mkScalar
#' @include Pvalues-S4class.R
#' @seealso fitBetaUniformParameters BioNet::scoreNodes
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
#' @references Beisser D, Klau GW, Dandekar T, Müller T, Dittrich MT. BioNet: An R-Package for the functional analysis of biological networks. Bioinformatics. 2010
setGeneric("betaUniformScore",
           valueClass = "numeric",
           function(x, betaUniformFit = NULL, FDR = 5E-2) { standardGeneric("betaUniformScore")},
           signature = c("x", "betaUniformFit", "FDR") )


#' @describeIn betaUniformScore Score P-values against an explicit beta-uniform model object
setMethod("betaUniformScore",
          c(x="Pvalues", betaUniformFit = "BetaUniformModel", FDR = "ScalarNumeric"),
          function(x, betaUniformFit, FDR = 5E-2) {
            
            FDRthreshold <- pValueThresholdAtConfidence(betaUniformFit, FDR)
            
            return( (betaUniformFit@a - 1)*(log(x@.Data) - log(FDRthreshold) ) )
          })

setOldClass("igraph")
#' @describeIn betaUniformScore Score nodes in an igraph; there must be a 'pValue'-like attribute (pVal, p.value etc.)
setMethod("betaUniformScore",
          c(x="igraph", betaUniformFit = "ANY", FDR = "ANY"),
          function(x, betaUniformFit = NULL, FDR = 5E-2) { 
            
            x <- validateIgraphWithPvalues(x) #P-value attr is augmented as PVALUE vertex attribute
            
            callGeneric( structure(V(x)$PVALUE, names = V(x)$name), betaUniformFit, FDR) })

setMethod("betaUniformScore",
          c(x="numeric", betaUniformFit = "ANY", FDR = "ANY"),
          function(x, betaUniformFit = NULL, FDR = 5E-2) { callGeneric( new("Pvalues", .Data = x) , betaUniformFit, FDR) })

setMethod("betaUniformScore",
          c(x="ANY", betaUniformFit = "ANY", FDR = "numeric"),
          function(x, betaUniformFit = NULL, FDR) {  callGeneric(x, betaUniformFit, mkScalar(FDR))  })

#' @describeIn betaUniformScore Score P-values against a Beta-uniform model fit on the fly.
setMethod("betaUniformScore",
          c(x="Pvalues", betaUniformFit = "NULL", FDR = "ANY"),
          function(x, betaUniformFit, FDR = 5E-2) { callGeneric(x, fitBetaUniformMixtureDistribution(x@.Data), FDR) })

setOldClass("bum") # from BioNet
#' @describeIn betaUniformScore Score P-values against a beta-uniform model S3 object from BioNet
setMethod("betaUniformScore",
          c(x="ANY", betaUniformFit = "bum", FDR = "ANY"),
          function(x, betaUniformFit, FDR = 5E-2){ callGeneric(x, convertFromBUM(betaUniformFit), FDR) })
