#' Fit a Beta-Uniform decomposition model to the provided p-values
#' 
#' Following the approach of Morris & Pounds, a distribution of P-values can be decomposed as a mixture of two beta distributions,
#' one with shape paramter (\eqn{\alpha}) of 1 - modelling noise under the null-hypothesis - and another with a variable (a), modelling
#' signal. One can view this as modeling P-values as a random process where \eqn{P ~ (1-\lambda)\beta(a,1) + \lambda\beta(1,1)}
#' 
#' Designed to be a drop-in replacement for the BioNet fitBumModel
#' 
#' @param pVals A numeric vector of P-values for which to perform Beta-Uniform decomposition
#' @param nStarts How many repeats of the fitting routine should the routine run?
#' @return fb A list of three components: a - the shape parameter, lambda - the mixture parameter and negLL - the negative log-likelihood of the data under the model
#' @importFrom MASS fitdistr
#' @importFrom stats na.omit coef logLik
#' @importFrom purrr map
#' @export
#' @seealso betaUniformScore BioNet::fitBumModel
#' 
#' @examples 
#' data(pvalues)
#' 
#' fbMod <- fitBetaUniformParameters(pvalues$pValue_diffex)
#' 
#' #Note that you must subtract one from the bayes-factor in order to achieve equivilence with BioNet::scoreNodes
#' betaUnifScored <- betaUniformScore(pvalues$pValue_diffex, fbMod, FDR = 0.01) - 1
#' 
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
fitBetaUniformParameters <- function(pVals, nStarts = 10L){
  
  validateSingleInteger(nStarts)
  assert_that(all(is.pval(na.omit(pVals))), msg = "Expecting P-values as input")
  
  #Fit beta-uniform distribution from a number of random starting points
  fittedBetaUniformModelParam_manyRun <-  map(1:nStarts,
                                              ~ tryCatch({ fitMod <- fitdistr(na.omit(pVals),
                                                                              betaUniformDensity,
                                                                              start = list(a = runif(1, min = 1E-5, max = 3), 
                                                                                           lambda = runif(1, min = 1E-5, max = 1-1E-5) ),
                                                                              lower = c(1E-5, 1E-5),
                                                                              upper = c(3, 1-1E-5),
                                                                              method = "L-BFGS-B")},
                                                         error = function(e){
                                                           NAlogLik <- list(estimate = as.numeric(c(a = NA, lambda = NA)),
                                                                            sd = as.numeric(c(a = NA, lambda = NA)),
                                                                            vcov = matrix(c(0,0,0,0),nrow = 2, dimnames = list(c("a","lambda"),c("a","lambda"))),
                                                                            loglik = as.numeric(NA),
                                                                            n = as.integer(NA))
                                                           class(NAlogLik) <- "fitdistr"
                                                           return(NAlogLik)} ) )
  
  bestFit <- which.max(map_dbl(fittedBetaUniformModelParam_manyRun, logLik))
  
  if (length(bestFit) == 0) { stop("Beta-Uniform model could not be fitted to data") }
  
  fittedBetaUniformModelParam <- fittedBetaUniformModelParam_manyRun[[bestFit]]
  
  fb <- list(a = coef(fittedBetaUniformModelParam)["a"] %>% as.numeric,
             lambda = coef(fittedBetaUniformModelParam)["lambda"] %>% as.numeric,
             pvalues = pVals,
             negLL = logLik(fittedBetaUniformModelParam)[1] %>% as.numeric)
  class(fb) <- "bum"
  
  if( isTRUE(all.equal(fb$a, 1E-5, tolerance = 1E-4)) | isTRUE(all.equal(fb$a, 3, tolerance = 1E-4)) |
      isTRUE(all.equal(fb$lambda, 1E-5, tolerance = 1E-4)) | isTRUE(all.equal(fb$lambda, 1-1E-5, tolerance = 1E-4)) ){
    
    warning("One or both parameters are on the limit of the defined parameter space")
  }
  
  return(fb)
}

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
  
  #TODO check betaUniformFit object
  
  return(log( pbeta(pVals, betaUniformFit$a, 1), pbeta(rep(FDR, length(pVals)), betaUniformFit$a, 1)) - 1)
}
