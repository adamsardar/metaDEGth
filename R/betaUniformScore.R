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
#' @seealso fitBetaUniformParameters BioNet::scoreNodes
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
#' @references Beisser D, Klau GW, Dandekar T, Müller T, Dittrich MT. BioNet: An R-Package for the functional analysis of biological networks. Bioinformatics. 2010
setGeneric("betaUniformScore",
           function(x, betaUniformFit = NULL, FDR = 5E-2) { standardGeneric("betaUniformScore")},
           signature = c("x", "betaUniformFit", "FDR") )

#' @describeIn betaUniformScore Score numeric values, assumed to be P-values
setMethod("betaUniformScore",
          c(x="numeric", betaUniformFit = "ANY", FDR = "numeric"),
          function(x, betaUniformFit = NULL, FDR = 5E-2) {
            # validate is pvalue, then convert
            
            callGeneric(structure(x, class = "pvalues"), betaUniformFit, FDR)
          })

setOldClass("igraph")

#' @describeIn betaUniformScore Score nodes in a graph; there must be a 'pValue' attribute
#' @import igraph
setMethod("betaUniformScore",
          c(x="igraph", betaUniformFit = "ANY", FDR = "numeric"),
          function(x, betaUniformFit = NULL, FDR = 5E-2) {
            # validate is pvalue, then convert
            
            callGeneric(structure(V(x)$pValue, names=V(x)$pValue, class = "pvalues"), betaUniformFit, FDR)
          })

setOldClass("pvalues")

#' @describeIn betaUniformScore Score P-values. Beta-uniform model fit on the fly.
setMethod("betaUniformScore",
          c(x="pvalues", betaUniformFit = "missing", FDR = "numeric"),
          function(x, betaUniformFit, FDR = 5E-2) {})

setOldClass("bum")

#' @describeIn betaUniformScore Score P-values against an explicit beta-uniform model
setMethod("betaUniformScore",
          c(x="pvalues", betaUniformFit = "bum", FDR = "numeric"),
          function(x, betaUniformFit, FDR) {

            return(log(pbeta(x, betaUniformFit$a, 1)) - rep(log(pbeta(FDR, betaUniformFit$a, 1)), times = length(x)) )
          })