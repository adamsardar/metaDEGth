#' Object representation of a beta-uniform model fit to a P-value distribution
#' 
#' Details the model P ~ (1-λ)β(a,1) + λβ(1,1), fit to a P-value set
#' 
#' @slot a Shape parameter of distribution β(a,1)
#' @slot lambda Mixture parameter
#' @slot pvalues A numeric vector of the values fit by the model
#' @slot negLL The negative log-likelihood of the model
#' 
#' @importClassesFrom  Biobase ScalarNumeric
#' @export
setClass("BetaUniformModel",
         representation(a = "ScalarNumeric",
                        lambda = "ScalarNumeric",
                        pvalues = "numeric",
                        negLL = "ScalarNumeric") )


#' Simple summary printer for a beta-uniform fit
#' 
#' @export
setMethod("show", 
          c(object="BetaUniformModel"),
          function(object) {
            
            cat("Beta-Uniform-Mixture (BUM) model\n\n")
            cat(paste(length(object@pvalues), "pvalues fitted\n\n"))
            cat(sprintf("Mixture parameter (lambda):\t%1.3f\n", object@lambda))
            cat(sprintf("shape parameter (a): \t\t%1.3f\n", object@a))
            cat(sprintf("log-likelihood:\t\t\t%.1f\n", -object@negLL)) 
            
            invisible(object)
          })

# A simple adapter from the class produced by BioNet
#' @export 
setGeneric("convertFromBUM", valueClass = "BetaUniformModel", function(obj) standardGeneric("convertFromBUM"))

setOldClass("bum")
setMethod("convertFromBUM", signature = c(obj = "bum"),
          function(obj){ new("BetaUniformModel", 
                             a = mkScalar(obj$a), 
                             lambda = mkScalar(obj$lambda),
                             pvalues = obj$pvalues,
                             negLL = mkScalar(obj$negLL))} )

## Distribution methods ##

## Group generic for extracting protions of a distribution at particular cutoffs
setGroupGeneric("densityPartitioning",
                knownMembers  = c("truePositiveFraction", "falsePositiveFraction", "trueNegativeFraction", "falseNegativeFraction"),
                def = function(obj, pValueThreshold){ NULL})

# Upper bound on the the noise component of mixture distribution [Equation (3)]: π in the equations below
setGeneric("noiseFractionUpperBound", valueClass = "ScalarNumeric", function(obj) standardGeneric("noiseFractionUpperBound"))

setMethod("noiseFractionUpperBound",
          signature(obj = "BetaUniformModel"),
          function(obj){ return(obj@lambda + (1-obj@lambda)*obj@a) } )


# True positive rate from M&P: A in table 1 and equation 5
#' @export 
setGeneric("truePositiveFraction",
           group = "densityPartitioning",
           function(obj, pValueThreshold){ standardGeneric("truePositiveFraction") } )

setMethod("truePositiveFraction",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) {
            
            Ftau <- obj@lambda*pValueThreshold + (1-obj@lambda)*pValueThreshold^obj@a
            return( Ftau - noiseFractionUpperBound(obj)*pValueThreshold ) }) # p_A(τ)= F(τ)−πτ
      

# FALSE positive rate from M&P: C table 1 and equation 5
#' @export 
setGeneric("falsePositiveFraction",  valueClass = "ScalarNumeric",
           group = "densityPartitioning",
           function(obj, pValueThreshold){ standardGeneric("falsePositiveFraction") } )

setMethod("falsePositiveFraction",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) { return( noiseFractionUpperBound(obj)*pValueThreshold ) }) # p_C= πτ 


# False negative rate from M&P: B table 1 and equation 5
#' @export 
setGeneric("falseNegativeFraction",  valueClass = "ScalarNumeric",
           group = "densityPartitioning",
           function(obj, pValueThreshold){ standardGeneric("falseNegativeFraction") } )

setMethod("falseNegativeFraction",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) {
            
            Ftau <- obj@lambda*pValueThreshold + (1-obj@lambda)*pValueThreshold^obj@a
            
            return(1 - Ftau - (1-pValueThreshold)*noiseFractionUpperBound(obj)   ) }) # p_B(τ)=1−F(τ)−(1−τ)π


# True negative rate from M&P: D table 1 and equation 5
#' @export 
setGeneric("trueNegativeFraction",  valueClass = "ScalarNumeric",
           group = "densityPartitioning",
           function(obj, pValueThreshold){ standardGeneric("trueNegativeFraction") } )

setMethod("trueNegativeFraction",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) {
            
            Ftau <- obj@lambda*pValueThreshold + (1-obj@lambda)*pValueThreshold^obj@a
            
            return( (1-pValueThreshold)*noiseFractionUpperBound(obj) ) }) # p_D(τ)=(1−τ)π

# False discovery rate - of BH fame. Equation 6 from M&P

setGeneric("FDRestimate",  valueClass = "ScalarNumeric",
           group = "densityPartitioning",
           function(obj, pValueThreshold){ standardGeneric("FDRestimate") } )

setMethod("FDRestimate",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) {
            return( falsePositiveFraction(obj, pValueThreshold)/( truePositiveFraction(obj, pValueThreshold)  + falsePositiveFraction(obj, pValueThreshold) )  ) }) # P_c/(P_a + P_c)

# pValueThreshold estimate - estimates a threshdold such that FDRestimate(obj, pValueThreshold) <= confidenceLevel. Equatuion 7 in Morris & Pounds
setGeneric("pValueThresholdAtConfidence",  valueClass = "ScalarNumeric",
           function(obj, confidenceLevel){ standardGeneric("pValueThresholdAtConfidence") } )

setMethod("pValueThresholdAtConfidence",
          signature = signature(obj = "ANY", confidenceLevel = "numeric"),
          function(obj, confidenceLevel) { return(callGeneric(obj, mkScalar(confidenceLevel)))  })

setMethod("pValueThresholdAtConfidence",
          signature = signature(obj = "BetaUniformModel", confidenceLevel = "ScalarNumeric"),
          function(obj, confidenceLevel) {
            
            return( (( noiseFractionUpperBound(obj) - confidenceLevel*obj@lambda )/(confidenceLevel*(1-obj@lambda)))^(1/(obj@a - 1))  ) }) 


### Group-level general methods available across all density paritioning methods. This allows us to have missing values and for it to auto-convert to a NumericScalar
setMethod("densityPartitioning",
          signature = signature(obj = "ANY", pValueThreshold = "missing"),
          function(obj, pValueThreshold) { message("Using a threshold of 0.05"); return(callGeneric(obj, 0.05 )) })

#' @importFrom Biobase mkScalar
setMethod("densityPartitioning",
          signature = signature(obj = "ANY", pValueThreshold = "numeric"),
          function(obj, pValueThreshold) { return(callGeneric(obj, mkScalar(pValueThreshold)))  })

