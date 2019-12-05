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
#' betaUnifScored <- betaUniformScore(pvalues$pValue_diffex, fbMod, FDR = 0.01)
#' 
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
              ~ fitdistr(x = na.omit(pValues), densfun = betaUniformDensity,
                         start = list(a = runif(1, min = epsilon, max = 3), 
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
            if(sd(modelLLHs) > 5) warning("log-likelihoods are very diverse - it might be worth increasing the 'nStarts' parameter.")
            
            bestFittedBetaUniformModelParam <- fittedBetaUniformModelParams[!failedOptimisationRuns][[bestFit]]
            
            betaUnifMod <- new("BetaUniformModel", 
                               a = mkScalar(coef(bestFittedBetaUniformModelParam)["a"]),
                               lambda = mkScalar(coef(bestFittedBetaUniformModelParam)["lambda"]),
                               pvalues = pValues,
                               negLL = mkScalar(logLik(bestFittedBetaUniformModelParam)[1]))
            
            #Inspect fb for whether any parameters lie on extremas
            if( isTRUE(all.equal(betaUnifMod@a, epsilon, tolerance = epsilon)) |
                isTRUE(all.equal(betaUnifMod@a, 3, tolerance = epsilon)) |
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

