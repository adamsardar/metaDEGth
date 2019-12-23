#' Fo each P-value, compute the log likelihood ratio of originating from signal or noise component
#' 
#' Relative to some false discovery threshold, compute the relative ratio of probabilities and then take the natural logarithm.
#' 
#' This function is somewhat analagous to BioNet::scoreNodes
#' 
#' @param betaUniformFit (optional) A beta-uniform model of the P-value distribution. Set to NULL to autofit on the fly.
#' @param FDR The tolerable fraction of false positives within the set of positive scoring values. See Morris & Pounds (2003)
#' @param pVals Vector of P-values to score
#' 
#' @return betaUniformScores A vector of P-Value scores
#' @export
#' @importFrom igraph V
#' @importFrom stats pbeta
#' @importFrom Biobase mkScalar
#' @include Pvalues_S4Class.R
#' @seealso fitBetaUniformParameters BioNet::scoreNodes
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
#' @references Dittrich MT, Klau GW, Rosenwald A, Dandekar T, Müller T. Identifying functional modules in protein-protein interaction networks: An integrated exact approach. Bioinformatics. 2008
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

#' @describeIn betaUniformScore Score p-values from a BetaUniformModel object against inself
setMethod("betaUniformScore",
          c(x="BetaUniformModel", betaUniformFit = "missing", FDR = "ANY"),
          function(x, betaUniformFit, FDR = 5E-2) {  callGeneric(x@pvalues, x, FDR)  })

#' @describeIn betaUniformScore Score P-values against a Beta-uniform model fit on the fly.
setMethod("betaUniformScore",
          c(x="Pvalues", betaUniformFit = "NULL", FDR = "ANY"),
          function(x, betaUniformFit=NULL, FDR = 5E-2) { callGeneric(x, fitBetaUniformMixtureDistribution(x@.Data), FDR) })

setOldClass("bum") # from BioNet
#' @describeIn betaUniformScore Score P-values against a beta-uniform model S3 object from BioNet
setMethod("betaUniformScore",
          c(x="ANY", betaUniformFit = "bum", FDR = "ANY"),
          function(x, betaUniformFit, FDR = 5E-2){ callGeneric(x, convertFromBUM(betaUniformFit), FDR) })
