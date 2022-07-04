#' Object representation of a beta-uniform model fit to a P-value distribution
#'
#' Details the model P ~ (1-λ)β(a,1) + λβ(1,1), fit to a P-value set
#'
#' @param pValueThreshold p-value threshold
#' @param object beta-uniform model
#' @slot a Shape parameter of distribution β(a,1)
#' @slot lambda Mixture parameter
#' @slot pvalues A numeric vector of the values fit by the model
#' @slot negLLH The negative log-likelihood of the model
#' @export
setClass("BetaUniformModel",
         representation(a = "ScalarNumeric",
                        lambda = "ScalarNumeric",
                        pvalues = "numeric",
                        negLLH = "ScalarNumeric") )


#' Simple summary printer for a beta-uniform fit
#'
#' @param object BetaUniformModel object
#' @export
setMethod("show",
          c(object = "BetaUniformModel"),
          function(object) {

            cat("Beta-Uniform-Mixture (BUM) model\n\n")
            cat(paste(length(object@pvalues), "pvalues fitted\n\n"))
            cat(sprintf("Mixture parameter (lambda):\t%1.3f\n", object@lambda))
            cat(sprintf("Beta shape parameter (a): \t%1.3f\n", object@a))
            cat(sprintf("Negative log-likelihood:\t%.1f\n", object@negLLH))

            invisible(object)
          })

#' A simple adapter from the class produced by BioNet
#'
#' @param obj data object
#' @rdname convertFromBUM
#' @export
setGeneric("convertFromBUM", valueClass = "BetaUniformModel", function(obj) standardGeneric("convertFromBUM"))

setOldClass("bum")
#' @describeIn convertFromBUM typed method for dispatch
setMethod("convertFromBUM", signature = c(obj = "bum"),
          function(obj){
            new("BetaUniformModel",
            a = mkScalar(obj$a),
            lambda = mkScalar(obj$lambda),
            pvalues = obj$pvalues,
            negLLH = mkScalar(obj$negLL))})

## Distribution methods ##


#' True positive rate from M&P: A in table 1 and equation 5
#'
#' @param obj data object
#' @param pValueThreshold p-value threshold
#' @rdname truePositiveFraction
#' @export
setGeneric("truePositiveFraction",
           function(obj, pValueThreshold){ standardGeneric("truePositiveFraction") } )
#' @describeIn truePositiveFraction method for dispatch
setMethod("truePositiveFraction",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) {

            Ftau <- obj@lambda * pValueThreshold + (1-obj@lambda)*pValueThreshold^obj@a
            return(Ftau - noiseFractionUpperBound(obj)*pValueThreshold)}) # p_A(τ)= F(τ)−πτ
#' @describeIn truePositiveFraction using a threshold of 0.05
setMethod("truePositiveFraction",
          signature = signature(obj = "ANY", pValueThreshold = "missing"),
          function(obj, pValueThreshold) { message("Using a threshold of 0.05"); return(callGeneric(obj, 0.05 )) })

#' @importFrom Biobase mkScalar
#' @describeIn truePositiveFraction threshold defined by the user
setMethod("truePositiveFraction",
          signature = signature(obj = "ANY", pValueThreshold = "numeric"),
          function(obj, pValueThreshold) {
            return(callGeneric(obj, mkScalar(pValueThreshold)))})


#'  FALSE positive rate from M&P: C table 1 and equation 5
#' @inheritParams falseNegativeFraction
#' @rdname falsePositiveFraction
#' @export
setGeneric("falsePositiveFraction",  valueClass = "ScalarNumeric",
           function(obj, pValueThreshold){
            standardGeneric("falsePositiveFraction")})
#' @describeIn falsePositiveFraction method for dispatch
setMethod("falsePositiveFraction",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) {
            return( noiseFractionUpperBound(obj)*pValueThreshold)}) # p_C= πτ
#' @describeIn falsePositiveFraction compute false positive rate using a
#' threshold of 0.05
setMethod("falsePositiveFraction",
          signature = signature(obj = "ANY", pValueThreshold = "missing"),
          function(obj, pValueThreshold) {
            message("Using a threshold of 0.05"); return(callGeneric(obj, 0.05 ))})

#' @importFrom Biobase mkScalar
#' @describeIn falsePositiveFraction compute false positive rate using threshold
#' defined by the user, numeric
setMethod("falsePositiveFraction",
          signature = signature(obj = "ANY", pValueThreshold = "numeric"),
          function(obj, pValueThreshold) {
            return(callGeneric(obj, mkScalar(pValueThreshold)))})


#' False negative rate from M&P: B table 1 and equation 5
#'
#' @param obj data object
#' @param pValueThreshold p-value threshold
#' @rdname falseNegativeFraction
#' @export
setGeneric("falseNegativeFraction",  valueClass = "ScalarNumeric",
           function(obj, pValueThreshold){
            standardGeneric("falseNegativeFraction")})
#' @describeIn falseNegativeFraction method for dispatch
setMethod("falseNegativeFraction",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) {

            Ftau <- obj@lambda * pValueThreshold + (1-obj@lambda)*pValueThreshold^obj@a

            return(1 - Ftau - (1-pValueThreshold) * noiseFractionUpperBound(obj) ) }) # p_B(τ)=1−F(τ)−(1−τ)π
#' @describeIn falseNegativeFraction compute false negative rate using a
#' threshold of 0.05
setMethod("falseNegativeFraction",
          signature = signature(obj = "ANY", pValueThreshold = "missing"),
          function(obj, pValueThreshold) { message("Using a threshold of 0.05"); return(callGeneric(obj, 0.05 )) })

#' @importFrom Biobase mkScalar
#' @describeIn falseNegativeFraction compute false negative rate using
#' numeric threshold defined by the user
setMethod("falseNegativeFraction",
          signature = signature(obj = "ANY", pValueThreshold = "numeric"),
          function(obj, pValueThreshold) { return(callGeneric(obj, mkScalar(pValueThreshold)))  })



#' True negative rate from M&P: D table 1 and equation 5
#' @inheritParams falseNegativeFraction
#' @rdname trueNegativeFraction
#' @export
setGeneric("trueNegativeFraction",  valueClass = "ScalarNumeric",
           function(obj, pValueThreshold){
            standardGeneric("trueNegativeFraction")})
#' @describeIn trueNegativeFraction method for dispatch
setMethod("trueNegativeFraction",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) {
            Ftau <- obj@lambda * pValueThreshold + (1-obj@lambda)*pValueThreshold^obj@a

            return((1-pValueThreshold)*noiseFractionUpperBound(obj))}) # p_D(τ)=(1−τ)π
#' @describeIn trueNegativeFraction compute true negative rate using a
#' threshold of 0.05
setMethod("trueNegativeFraction",
          signature = signature(obj = "ANY", pValueThreshold = "missing"),
          function(obj, pValueThreshold) {
            message("Using a threshold of 0.05"); return(callGeneric(obj, 0.05 ))})

#' @importFrom Biobase mkScalar
#' @describeIn trueNegativeFraction compute true negative rate using
#' numeric threshold defined by the user
setMethod("trueNegativeFraction",
          signature = signature(obj = "ANY", pValueThreshold = "numeric"),
          function(obj, pValueThreshold) {
            return(callGeneric(obj, mkScalar(pValueThreshold)))})


#' Upper bound on the the noise component of mixture distribution
#'
#' π in equation 3 from Morris & Pounds
#' @param obj data object
#' @rdname noiseFractionUpperBound
#' @export
setGeneric("noiseFractionUpperBound", valueClass = "ScalarNumeric", function(obj) standardGeneric("noiseFractionUpperBound"))
#' @describeIn noiseFractionUpperBound method for dispatch
setMethod("noiseFractionUpperBound",
          signature(obj = "BetaUniformModel"),
          function(obj){
            return(obj@lambda + (1-obj@lambda)*obj@a)})



#' False discovery rate - of BH fame. Equation 6 from M&P
#'
#' @inheritParams falseNegativeFraction
#' @rdname FDRestimate
#' @export
setGeneric("FDRestimate",  valueClass = "ScalarNumeric",
           function(obj, pValueThreshold){
            standardGeneric("FDRestimate")} )
#' @describeIn FDRestimate method for BetaUniformModel object
setMethod("FDRestimate",
          signature = signature(obj = "BetaUniformModel", pValueThreshold = "ScalarNumeric"),
          function(obj, pValueThreshold) {
            return( falsePositiveFraction(obj, pValueThreshold)/( truePositiveFraction(obj, pValueThreshold)  + falsePositiveFraction(obj, pValueThreshold) )  ) }) # P_c/(P_a + P_c)



#' pValueThreshold estimate - estimates a threshdold such that FDRestimate(obj, pValueThreshold) <= confidenceLevel. Equatuion 7 in Morris & Pounds
#'
#' @param obj data object
#' @param confidenceLevel confidence level
#' @rdname pValueThresholdAtConfidence
#' @export
setGeneric("pValueThresholdAtConfidence",  valueClass = "ScalarNumeric",
           function(obj, confidenceLevel){
            standardGeneric("pValueThresholdAtConfidence")})
#' @describeIn pValueThresholdAtConfidence method for dispatch
setMethod("pValueThresholdAtConfidence",
          signature = signature(obj = "ANY", confidenceLevel = "numeric"),
          function(obj, confidenceLevel) {
            return(callGeneric(obj, mkScalar(confidenceLevel)))})
#' @describeIn pValueThresholdAtConfidence p-value estimate for scalar numeric
setMethod("pValueThresholdAtConfidence",
          signature = signature(obj = "BetaUniformModel", confidenceLevel = "ScalarNumeric"),
          function(obj, confidenceLevel) {

            return( (( noiseFractionUpperBound(obj) - confidenceLevel*obj@lambda )/(confidenceLevel*(1-obj@lambda)))^(1/(obj@a - 1)))})



