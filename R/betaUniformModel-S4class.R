
#' @importClassesFrom  Biobase ScalarNumeric
#' @importFrom Biobase mkScalar
setClass("BetaUniformModel",
         representation(a = "ScalarNumeric",
                        lambda = "ScalarNumeric",
                        pvalues = "numeric",
                        negLL = "ScalarNumeric") )


#' Simple summary printer for a beta-uniform fit
#' 
#' 
#' 
#' @export
#' @importClassesFrom methods show
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

setGeneric("convertFromBUM", valueClass = "BetaUniformModel", function(obj) standardGeneric("convertFromBUM"))

setOldClass("bum")
setMethod("convertFromBUM", signature = c(obj = "bum"),
          function(obj){ new("BetaUniformModel", 
                             a = mkScalar(obj$a), 
                             lambda = mkScalar(obj$lambda),
                             pvalues = obj$pvalues,
                             negLL = mkScalar(obj$negLL))} )

# Distribution methods

# Upper bound on the the noise component of mixture distribution [Equation (3)]: π in the equations below
setGeneric("noiseFractionUpperBound", valueClass = "ScalarNumeric", function(obj) standardGeneric("noiseFractionUpperBound"))

setMethod("noiseFractionUpperBound",
          signature(obj = "BetaUniformModel"),
          function(obj){ return(obj@lambda + (1-obj@lambda)*obj@a) } )


# True posive rate from M&P: A in table 1 and equation 5
setGeneric("truePositiveFraction",
           function(obj, pValueThreshold){ standardGeneric("truePositiveFraction") } )

setMethod("truePositiveFraction",
          signature = signature(obj = "ANY", pValueThreshold = "missing"),
          function(obj, pValueThreshold) { return(callGeneric(obj, 0.05)) })

setMethod("truePositiveFraction",
          signature = signature(obj = "ANY", pValueThreshold = "numeric"),
          function(obj, pValueThreshold) { return(callGeneric(obj, mkScalar(pValueThreshold)))  })

setMethod("truePositiveFraction",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) {
            
            Ftau <- obj@lambda*pValueThreshold + (1-obj@lambda)*pValueThreshold^obj@a
            return( Ftau - noiseFractionUpperBound(obj)*pValueThreshold ) # p_A(τ)= F(τ)−πτ
          })


# True posive rate from M&P: C table 1 and equation 5
setGeneric("falsePositiveFraction",  valueClass = "ScalarNumeric",
           function(obj, pValueThreshold){ standardGeneric("falsePositiveFraction") } )

setMethod("falsePositiveFraction",
          signature = signature(obj = "ANY", pValueThreshold = "missing"),
          function(obj, pValueThreshold) { return(callGeneric(obj, 0.05)) })

setMethod("falsePositiveFraction",
          signature = signature(obj = "ANY", pValueThreshold = "numeric"),
          function(obj, pValueThreshold) { return(callGeneric(obj, mkScalar(pValueThreshold)))  })

setMethod("falsePositiveFraction",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) { return( noiseFractionUpperBound(obj)*pValueThreshold ) }) # p_C= πτ 


# False negative rate from M&P: B table 1 and equation 5
setGeneric("falseNegativeFraction",  valueClass = "ScalarNumeric",
           function(obj, pValueThreshold){ standardGeneric("falseNegativeFraction") } )

setMethod("falseNegativeFraction",
          signature = signature(obj = "ANY", pValueThreshold = "missing"),
          function(obj, pValueThreshold) { return(callGeneric(obj, 0.05)) })

setMethod("falseNegativeFraction",
          signature = signature(obj = "ANY", pValueThreshold = "numeric"),
          function(obj, pValueThreshold) { return(callGeneric(obj, mkScalar(pValueThreshold)))  })

setMethod("falseNegativeFraction",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) {
            
            Ftau <- obj@lambda*pValueThreshold + (1-obj@lambda)*pValueThreshold^obj@a
            
            return(1 - Ftau - (1-pValueThreshold)*noiseFractionUpperBound(obj)   ) }) # p_B(τ)=1−F(τ)−(1−τ)π


# True negative rate from M&P: D table 1 and equation 5
setGeneric("trueNegativeFraction",  valueClass = "ScalarNumeric",
           function(obj, pValueThreshold){ standardGeneric("trueNegativeFraction") } )

setMethod("trueNegativeFraction",
          signature = signature(obj = "ANY", pValueThreshold = "missing"),
          function(obj, pValueThreshold) { return(callGeneric(obj, 0.05)) })

setMethod("trueNegativeFraction",
          signature = signature(obj = "ANY", pValueThreshold = "numeric"),
          function(obj, pValueThreshold) { return(callGeneric(obj, mkScalar(pValueThreshold)))  })

setMethod("trueNegativeFraction",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) {
            
            Ftau <- obj@lambda*pValueThreshold + (1-obj@lambda)*pValueThreshold^obj@a
            
            return( (1-pValueThreshold)*noiseFractionUpperBound(obj) ) }) # p_D(τ)=(1−τ)π

