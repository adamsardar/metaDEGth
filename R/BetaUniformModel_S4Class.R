#' Object representation of a beta-uniform model fit to a P-value distribution
#' 
#' Details the model P ~ (1-λ)β(a,1) + λβ(1,1), fit to a P-value set
#' 
#' @slot a Shape parameter of distribution β(a,1)
#' @slot lambda Mixture parameter
#' @slot pvalues A numeric vector of the values fit by the model
#' @slot negLLH The negative log-likelihood of the model
#' 
#' @importClassesFrom  Biobase ScalarNumeric
#' @export
setClass("BetaUniformModel",
         representation(a = "ScalarNumeric",
                        lambda = "ScalarNumeric",
                        pvalues = "numeric",
                        negLLH = "ScalarNumeric") )


#' Simple summary printer for a beta-uniform fit
#' 
#' @export
setMethod("show", 
          c(object="BetaUniformModel"),
          function(object) {
            
            cat("Beta-Uniform-Mixture (BUM) model\n\n")
            cat(paste(length(object@pvalues), "pvalues fitted\n\n"))
            cat(sprintf("Mixture parameter (lambda):\t%1.3f\n", object@lambda))
            cat(sprintf("Beta shape parameter (a): \t%1.3f\n", object@a))
            cat(sprintf("Negative log-likelihood:\t%.1f\n", object@negLLH)) 
            
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
                             negLLH = mkScalar(obj$negLL))} )

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

globalVariables(c("..density.."))
#' Plot summary for fitted Beta-Uniform models
#' 
#' Produce a list of summary ggplots. The default print method produces a grid panel, or you can extract and use them at your wish.
#' 
#' These should assist in diagnosing good or poor model fit.
#' 
#' @param x A BetaUniformModel fit of data
#' @param ... Ignored
#' @param showFit Flag to specify if you would like the Beta-Uniform model formula outputted as well (default: TRUE)
#' @param studyName Subtitle for plots (default: variable name of x passed in)
#' @param outputFormula Flag to specify if you would like the Beta-Uniform model formula outputted as well (default: TRUE)
#' @param outputParameters Flag to specify if the plot should list the fitted parameters for the model (defualt:TRUE)
#' @return A named list with three ggplots: p-value distribution; QQplot of P-values against a the fitted mixture model and liklihood surface of paramter estimates.
#' @import ggplot2
#' @import latex2exp
#' @export
plot.BetaUniformModel <- function(x, showFit = TRUE, outputParameters = TRUE, outputFormula = TRUE, studyName = NULL, ...){
  
  if(is.null(studyName)){studyName <- deparse(substitute(x))}
  
  validateBooleanFlag(outputParameters)
  validateBooleanFlag(outputFormula)
  validateBooleanFlag(showFit)
  validateSingleString(studyName)
  
  pVal_plot <- ggplot(data = NULL, aes(x = x@pvalues)) + 
    geom_histogram(aes(y=..density..),colour="black",fill="azure3",binwidth = 0.025, center = 0.0125) +
    labs(title  = "P-value distribution",
         subtitle = studyName,
         x = "P-Value",
         y = "Density") +
    theme_bw() + 
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) 
  
  if(showFit){
    
    pVal_plot <-   pVal_plot +
      stat_function(fun = dbetauniform, args = list(a = x@a, lambda = x@lambda), aes(colour = "Beta-Uniform Mixture"),size=1.5) +
      stat_function(fun = uniformComponentBound, args = list(a = x@a, lambda = x@lambda),  aes(colour = "Upper Bound On\nUniform Component"),size=1.5) +
      scale_colour_manual("Model Density", values = c("red", "black")) +
      theme(legend.position=c(0.8,0.7), plot.title = element_text(size=20))
  }
  
  if(outputFormula){
    
    FBmodelEqn <- TeX("Fitted Model: $f(p|\\alpha,\\lambda) = \\lambda + (1-\\lambda)p^{\\alpha-1}$",output = "character")
    pVal_plot <- pVal_plot + annotate("text",x=0.5,y=dbetauniform(0.35, x@a, x@lambda)+2,label = FBmodelEqn, parse=TRUE, size = 4)
  }
  
  pVal_qqPlot <- ggplot(data = NULL, aes(sample = x@pvalues@.Data)) +
    stat_qq(distribution = qbetauniform, dparams = list(a = x@a, lambda =x@lambda)) +
    theme_bw() + 
    labs(title  = "QQ plot of observations against fitted Beta-Uniform model",
         subtitle = studyName,
         x = "Theoretical Quantile",
         y = "Observed Quantile") +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted")
  
  paramSpace <- expand.grid(a = seq(0,1,0.1), lambda = seq(0,1,0.1)) %>% data.table
  
  pVals <- x@pvalues@.Data #slot access is expensive
  paramSpace[, LLH := sum(log( dbetauniform(pVals, a, lambda) )), by = .(a, lambda)]

  LLsurface_plot <- ggplot(paramSpace[LLH > 0], aes(a, lambda)) +
    geom_tile(aes(fill=LLH), colour = "grey50") +
    stat_contour(aes(z=LLH), colour = "lightgrey") +
    scale_fill_gradientn(colours = terrain.colors(20)) +
    geom_point(x = x@a, y = x@lambda, size = 3, colour = "lightgrey") +
    labs(title = "Likelihood surface and ML parameter estimate",
         subtitle = studyName,
         x = "a", y = "lambda") +
    theme_bw()
  
  if(outputParameters){
    
    eq <- TeX(str_c("$\\alpha =",round(x@a, digits = 2),
                    ", \\lambda =",round(x@lambda, digits = 2),
                    ", LLH =",round(-x@negLLH,digits = 2),"$"), output = "character")
    
    LLsurface_plot <- LLsurface_plot + annotate("text", x=0.5, y=0.5, label = eq, parse = TRUE, size = 4)
    pVal_plot <- pVal_plot + annotate("text", x=0.5, y=dbetauniform(0.35, x@a, x@lambda)+1, label = eq, parse = TRUE, size = 4)
  }
    
  return(structure(list(pValHist = pVal_plot,
                        qqPlot = pVal_qqPlot,
                        LLsurface = LLsurface_plot), class = "BetaUniformPlots"))
}

#' @import patchwork
print.BetaUniformPlots <- function(x, ...){
  
  return( plot(x$pValHist / (x$qqPlot | x$LLsurface)) )
}
