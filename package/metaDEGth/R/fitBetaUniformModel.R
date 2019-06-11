#' Fit a Beta-Uniform decomposition model to the provided p-values
#' 
#' Following the approach of Morris & Pounds, a distribution of P-values can be decomposed as a mixture of two beta distributions,
#' one with shape paramter (\eqn{\alpha}) of 1 - modelling noise under the null-hypothesis - and another with a variable (a), modelling
#' signal. One can view this as modeling P-values as a random process where \eqn{P ~ (1-\lambda)\beta(a,1) + \lambda\beta(1,1)}
#' 
#' Designed to be a drop-in replacement for BioNet::fitBumModel.
#' 
#' @param pVals A numeric vector of P-values for which to perform Beta-Uniform decomposition
#' @param nStarts How many repeats of the fitting routine should the routine run before the best LLH is picked? The default is 10, but in some cases (when there is no 'clean' P-value distribution), it might be worth running for more iterations.
#' @return fb A list of three components: a - the shape parameter, lambda - the mixture parameter and negLL - the negative log-likelihood of the data under the model
#' @importFrom MASS fitdistr
#' @importFrom stats na.omit coef logLik
#' @importFrom purrr map possibly
#' @export
#' @seealso betaUniformScore BioNet::fitBumModel
#' 
#' @examples 
#' data(pvalues)
#' 
#' fbMod <- fitBetaUniformMixtureDistribution(pvalues$pValue_diffex)
#' 
#' #Note that you must subtract one from the bayes-factor in order to achieve equivilence with BioNet::scoreNodes
#' betaUnifScored <- betaUniformScore(pvalues$pValue_diffex, fbMod, FDR = 0.01) - 1
#' 
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
fitBetaUniformMixtureDistribution <- function(pVals, nStarts = 10L){
  
  validateSingleInteger(nStarts)
  assert_that(all(is.pval(na.omit(pVals))), msg = "Expecting P-values as input") #TODO - replace with a proper validator
  
  epsilon <- 1E-5

  #Fit beta-uniform distribution from a number of uniform random starting points. We need to wrap the method safely  
  optimiseBetaUniformParams <- possibly(
    ~ fitdistr(x = na.omit(pVals), densfun = betaUniformDensity,
               start = list(a = runif(1, min = epsilon, max = 3), 
                            lambda = runif(1, min = epsilon, max = 1-epsilon) ),
               lower = c(epsilon, epsilon), upper = c(3, 1-epsilon),
               method = "L-BFGS-B") ,
    otherwise = NA_real_ )
  
  fittedBetaUniformModelParams <- map(1:nStarts, optimiseBetaUniformParams)

  failedOptimisationRuns <- ( sapply(fittedBetaUniformModelParams, class) != "fitdistr" )
  
  modelLLHs <- sapply(fittedBetaUniformModelParams[!failedOptimisationRuns], logLik)

  bestFit <- which.max(modelLLHs)
  if (length(bestFit) == 0) { stop("Beta-Uniform model could not be fitted to data") }
  if(sd(modelLLHs) > 5) warning("log-likelihoods are very diverse - it might be worth increasing the 'nStarts' parameter.")
  
  bestFittedBetaUniformModelParam <- fittedBetaUniformModelParams[!failedOptimisationRuns][[bestFit]]
  
  #Build an object compatible with BioNet's beta-uniform model class
  fb <- list(a = coef(bestFittedBetaUniformModelParam)["a"],
             lambda = coef(bestFittedBetaUniformModelParam)["lambda"],
             pvalues = na.omit(pVals),
             negLL = logLik(bestFittedBetaUniformModelParam)[1]) 
  
  fb %<>% lapply(as.numeric)
  class(fb) <- "bum"
  
  #Inspect fb for whether any parameters lie on extremas
  if( isTRUE(all.equal(fb$a, epsilon, tolerance = epsilon)) |
      isTRUE(all.equal(fb$a, 3, tolerance = epsilon)) |
      isTRUE(all.equal(fb$lambda, epsilon, tolerance = epsilon)) |
      isTRUE(all.equal(fb$lambda, 1-epsilon, tolerance = epsilon)) ){
    
    warning("One or both parameters are on the limit of the defined parameter space")
  }
  
  return(fb)
}

#' @export
print.bum <- function (x, ...) {
  
  cat("Beta-Uniform-Mixture (BUM) model\n\n")
  cat(paste(length(x$pvalues), "pvalues fitted\n\n"))
  cat(sprintf("Mixture parameter (lambda):\t%1.3f\n", x$lambda))
  cat(sprintf("shape parameter (a): \t\t%1.3f\n", x$a))
  cat(sprintf("log-likelihood:\t\t\t%.1f\n", -x$negLL))
}

#' Compute the adjusted log likelihood ratio score for beta-uniform model of P-values
#' 
#' Relative to some false discovery threshold, compute the relative ratio of probabilities and then take the natural logarithm.
#' 
#' Designed to be a drop-in replacement for the BioNet package scoreNodes function
#' 
#' @param fb The result of fitBetaUniformParameters method
#' @param FDR The tolerable false discovery rate. See Morris & Pounds (2003)
#' @param pVals Vector of P-values to score
#' 
#' @return betaUniformScores A vector of P-Value scores
#' @export
#' @importFrom stats pbeta
#' @seealso fitBetaUniformParameters BioNet::scoreNodes
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
betaUniformScore <- function(pVals, betaUniformFit, FDR = 5E-2){
  
  #TODO check betaUniformFit object via validator
  
  return(log( pbeta(pVals, betaUniformFit$a, 1), pbeta(rep(FDR, length(pVals)), betaUniformFit$a, 1)) - 1)
}
