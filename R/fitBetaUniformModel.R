globalVariables(c("..density.."))
#' Fit a Beta-Uniform decomposition model to the provided p-values
#' 
#' Following the approach of Morris & Pounds, a distribution of P-values can be decomposed as a mixture of two beta distributions,
#' one with shape paramter (\eqn{\alpha}) of 1 - modelling noise under the null-hypothesis - and another with a variable (a), modelling
#' signal. One can view this as modeling P-values as a random process where \eqn{P ~ (1-\lambda)\beta(a,1) + \lambda\beta(1,1)}
#' 
#' Analagous to BioNet::fitBumModel.
#' 
#' @param pValues A numeric vector of P-values for which to perform Beta-Uniform decomposition
#' @param nStarts How many repeats of the fitting routine should the routine run before the best LLH is picked? The default is 10, but in some cases (when there is no 'clean' P-value distribution), it might be worth running for more iterations.
#' @return A BetaUniformModel object, detailign the fit
#' @importFrom MASS fitdistr
#' @importFrom stats na.omit coef logLik sd
#' @importFrom purrr map possibly
#' @export
#' @seealso betaUniformScore BioNet::fitBumModel
#' 
#' @examples 
#' library(ggplot2)
#' data(siTAZdiffex_HFL1_diffexDT)
#' 
#' siTAZ1_BetaUniformMod <- fitBetaUniformMixtureDistribution(siTAZdiffex_HFL1_diffexDT$siTAZ1_pValue)
#' 
#' siTAZ1_BetaUniformMod
#' 
#' plot(siTAZ1_BetaUniformMod) # Produce some diagnostic plots
#' 
#' betaUnifScored <- betaUniformScore(siTAZ1_BetaUniformMod, FDR = 0.05)
#' qplot(betaUnifScored)
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
setGeneric("fitBetaUniformMixtureDistribution", 
           signature = signature("pValues","nStarts"),
           valueClass = "BetaUniformModel",
           function(pValues, nStarts = 10L){ standardGeneric("fitBetaUniformMixtureDistribution") })

#' @describeIn fitBetaUniformMixtureDistribution Typed method for dispatch
setMethod("fitBetaUniformMixtureDistribution",
          signature(pValues = "Pvalues", nStarts="ScalarInteger"),
          function(pValues, nStarts = 10L){
            
            epsilon <- 1E-5
            
            #Fit beta-uniform distribution from a number of uniform random starting points. We need to wrap the method safely  
            optimiseBetaUniformParams <- possibly(
              ~ fitdistr(x = na.omit(pValues), densfun = dbetauniform,
                         start = list(a = runif(1, min = epsilon, max = 1), 
                                      lambda = runif(1, min = epsilon, max = 1-epsilon) ),
                         lower = c(epsilon, epsilon),
                         upper = c(3, 1-epsilon),
                         method = "L-BFGS-B") ,
              otherwise = NA_real_ )
            
            fittedBetaUniformModelParams <- map(1:nStarts, optimiseBetaUniformParams)
            
            failedOptimisationRuns <- ( sapply(fittedBetaUniformModelParams, class) != "fitdistr" )
            
            modelLLHs <- sapply(fittedBetaUniformModelParams[!failedOptimisationRuns], logLik)
            
            bestFit <- which.max(modelLLHs)
            if (length(bestFit) == 0) { stop("Beta-Uniform model could not be fitted to data") }
            if(sd(modelLLHs) > 3) warning("log-likelihoods are very diverse - it might be worth increasing the 'nStarts' parameter, else accepting that parameter estimates are a little unstable")
            
            bestFittedBetaUniformModelParam <- fittedBetaUniformModelParams[!failedOptimisationRuns][[bestFit]]
            
            betaUnifMod <- new("BetaUniformModel", 
                               a = mkScalar(coef(bestFittedBetaUniformModelParam)["a"]),
                               lambda = mkScalar(coef(bestFittedBetaUniformModelParam)["lambda"]),
                               pvalues = pValues,
                               negLLH = mkScalar(-logLik(bestFittedBetaUniformModelParam)[1]))
            
            #Inspect fb for whether any parameters lie on extremas
            if( isTRUE(all.equal(betaUnifMod@a, epsilon, tolerance = epsilon)) |
                isTRUE(all.equal(betaUnifMod@a, 1, tolerance = epsilon)) |
                isTRUE(all.equal(betaUnifMod@lambda, epsilon, tolerance = epsilon)) |
                isTRUE(all.equal(betaUnifMod@lambda, 1-epsilon, tolerance = epsilon)) ){
              
              warning("One or both parameters are on the limit of the defined parameter space")
            }
            
            return(betaUnifMod)
          })

#' @describeIn fitBetaUniformMixtureDistribution Adapter method, attempts to convert numeric values to Pvalues
setMethod("fitBetaUniformMixtureDistribution",
          signature(pValues = "numeric", nStarts="ANY"),
          function(pValues, nStarts=10L){callGeneric( new("Pvalues", .Data = pValues), nStarts)})


setMethod("fitBetaUniformMixtureDistribution",
          signature(pValues = "Pvalues", nStarts="numeric"),
          function(pValues, nStarts){callGeneric(pValues, new("ScalarInteger",nStarts))})


## Density functions used in fits

#' Density function for Beta-uniform Model
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
#' @importFrom stats dbeta
dbetauniform <- function(x, a, lambda, log = FALSE){
  
  # The beta-uniform decomposition can be written as a mixture distribution:
  # (1-lambda)*dbeta(x,a,1) + lambda*dbeta(x,1,1)
  
  #However, dbeta(1,1) = duniform() = 1 everywhere, so there is no point wasting computation.
  # (1-lambda)*dbeta(x,a,1) + lambda
  
  # Further, the beta-distribution simplifies to a Kumaraswamy distribution when one of the Beta distribution parameters is unity. This is an order of magnitude faster than the first mixture distribution
  # See https://en.wikipedia.org/wiki/Kumaraswamy_distribution
  
  if(log){
    log((1-lambda)* a * x^(a-1)  + lambda)
  }else{
    (1-lambda)* a * x^(a-1)  + lambda
  }
}



#' Quantile function for Beta-uniform model
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
#' @importFrom stats qbeta
qbetauniform <- function(p, a, lambda){
  
  # A simple translation of the beta-uniform mixture distribution 
  # (1-lambda)*qbeta(p,a,1) + lambda*qbeta(p,1,1) 
  
  # This simplifies substatially. Beta(1,1) = Unif(), and the CDF of uniform at x is x. Beta(a,1) is Kumaraswamy(a,1), whose CDF is 1-(1-x^a).
  return((1-lambda)*qbeta(p,a,1) + lambda*p )
}


#' Random draw function for Beta-uniform model
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
#' @importFrom stats runif rbeta rbinom
#' @export
rbetauniform <- function(n, a, lambda){
  
  validateSingleInteger(n)
  validateSinglePositiveDefinite(n)
  
  validateSinglePositiveDefinite(a)
  
  validateNumericCutoff(lambda)
  validatePvalues(lambda)
  validateSinglePositiveDefinite(lambda)
  
  # A simple translation of the beta-uniform mixture distribution 
  # (1-lambda)*qbeta(p,a,1) + lambda*qbeta(p,1,1) 
  nUnifDraws <- rbinom(1, n, lambda) # The number of draws from unfirom will be binomially distributed
  
  # Create a random permutation of the signal and noise components stapled together
  betaUniformDraws <- sample(c(runif(nUnifDraws),
                               rbeta(n = (n-nUnifDraws),  shape1 = a, shape2 = 1)))
  
  return(betaUniformDraws)
}

#' Upper bound on Uniform component of Beta-uniform model
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
uniformComponentBound <- function(x, a, lambda){
  
  # A literal translation of the maths from Morris & Pounds would be:
  # ((1-lambda)*a + lambda)*dunif(x)
  
  # But dunif is 1 everywhere (in the bounds) by definition.
  return( ((1-lambda)*a + lambda) )
}

