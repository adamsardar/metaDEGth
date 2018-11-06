globalVariables(c("..density.."))
#' Plot histograms of p-values and, optionally, fit a Beta-Uniform model as well 
#' 
#' A helpful utility function that will plot histograms of P-values and, optionally, fit a beta-uniform model
#' 
#' @import ggplot2
#' @import latex2exp
#' @importFrom stringr str_c
#' 
#' @param DT A data.table containing the P-values of interest
#' @param pValAttr The column name in the data.table containing P-values (default: "P.Value)
#' @param fitBetaUniform An option logical flag specifying whter or not to fit and plot the beta uniform model of Paunds & Morris (default: TRUE)
#' @param outputFormula Flag to specify if you would like the Beta-Uniform model formula outputted as well (default: FALSE)
#' @param outputParameters Flag to specify if the plot should list the fitted parameters for the model (defualt:TRUE)
#' @param FBmodel Allows a user to provide a pre-fitted beta-uniform model of Pvalues (default:NULL)
#' 
#' @return A list of plots and, if requested, the Beta-Uniform model fit
#' 
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
#' @references Beisser, D., Klau, G. W., Dandekar, T., Muller, T., & Dittrich, M. T. (2010). BioNet: an R-Package for the functional analysis of biological networks. Bioinformatics
#' 
#' @importFrom stats na.omit
#' @export
pvalHist <- function(DT,pValAttr = "P.Value", fitBetaUniform = TRUE, outputParameters = TRUE, outputFormula = FALSE, FBmodel = NULL){
  
  studyName <- deparse(substitute(DT))
  
  assert_that(is.data.table(DT),
              pValAttr %in% colnames(DT),
              all(is.pval(DT[,na.omit(get(pValAttr))])),
              msg = "DT muse be a data.table with the specified P-value attribute")
  
  validateBooleanFlag(fitBetaUniform)
  
  #Need to pass the enviroment around so that the "get" command below works
  .e <- environment()
  pVal_plot <- DT[!is.na(get(pValAttr))] %>% ggplot(aes(x = get(pValAttr)),y = ..density.., environment = .e) + 
    geom_histogram(aes(y=..density..),colour="black",fill="azure3",binwidth = 0.01) +  
    ggtitle(studyName) + 
    theme_bw() + 
    xlab(pValAttr) 
  
  fb <- NULL
  
  #User provided model
  if(!is.null(FBmodel)){
    if(class(FBmodel) != "bum") stop("Model of unusable class provided to routine: ", class(FBmodel))
    #TODO a more thorough type check here please
    
    fb <- FBmodel
    fitBetaUniform <- TRUE
  }
  
  if(fitBetaUniform){
    
    if(is.null(fb)){ fb <- fitBetaUniformParameters(DT[!is.na(get(pValAttr)),get(pValAttr)]) }
    
    #Density function for Beta-uniform Model
    betaUniformConvFunc <- function(x){betaUniformDensity(x,fb$a,fb$lambda)}
    #Quantile function for Beta-uniform model
    qbetaUniformConvFunc <- function(p){ qbetaUniformFunc(p,fb$a,fb$lambda) }
    #Upper bound on Uniform component of Beta-uniform model
    uniformComponentConvFunc <- function(x){uniformComponentFunc(x,fb$a,fb$lambda) }
    
    pVal_plot <-   pVal_plot +
      stat_function(fun = betaUniformConvFunc, aes(colour = "Beta-Uniform Mixture"),size=1.5) +
      stat_function(fun = uniformComponentConvFunc, aes(colour = "Upper Bound On\nUniform Component"),size=1.5) +
      scale_colour_manual("Model Density", values = c("red", "black")) +
      scale_x_continuous(limits = c(0.0,1)) + 
      theme(legend.position=c(0.8,0.7), plot.title = element_text(size=22))
    
    if(outputParameters){
      
      fb_eqn <- function(fb){
        
        eq <- TeX(str_c("$\\alpha =",round(fb$a, digits = 2),
                        ", \\lambda =",round(fb$lambda, digits = 2),
                        ", log(LH) =",round(fb$negLL,digits = 2),"$"), output = "character")
        return(eq)
      }
      
      pVal_plot <- pVal_plot + annotate("text",x=0.5,y=betaUniformConvFunc(0.35)+1,label = fb_eqn(fb),parse = TRUE, size = 6)
    }
    
    if(outputFormula){
      FBmodelEqn <- TeX("Fitted Model: $f(p|\\alpha,\\lambda) = \\lambda + (1-\\lambda)p^{\\alpha-1}$",output = "character")
      pVal_plot <- pVal_plot + annotate("text",x=0.5,y=betaUniformConvFunc(0.35)+1.4,label = FBmodelEqn, parse=TRUE, size = 6)
    }
    
    pVal_qqPlot <- DT[!is.na(get(pValAttr))] %>% ggplot(aes(sample = get(pValAttr)), environment = .e) + 
      stat_qq(distribution = qbetaUniformConvFunc) + 
      ggtitle(studyName) + 
      theme_bw() + 
      xlab("Theoretical Quantile") +
      ylab("Observed Quantile") +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted")
    
    return(list(pValHist = pVal_plot, qqPlot = pVal_qqPlot, fb = fb))
    
  }else{
    
    return(list(pValHist = pVal_plot))
  }
}

#' Density function for Beta-uniform Model
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
#' @importFrom stats dbeta
betaUniformDensity <- function(x, a, lambda){
  
  (1-lambda)*dbeta(x,a,1) + lambda*dbeta(x,1,1)
}

#' Quantile function for Beta-uniform model
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
#' @importFrom stats qbeta
qbetaUniformFunc <- function(p, a, lambda){
  
  (1-lambda)*qbeta(p,a,1) + lambda*qbeta(p,1,1)
}

#' Upper bound on Uniform component of Beta-uniform model
#' @references Pounds, S., & Morris, S. W. (2003). Estimating the occurrence of false positives and false negatives in microarray studies by approximating and partitioning the empirical distribution of p-values. Bioinformatics
#' @importFrom stats dunif
uniformComponentFunc <- function(x, a, lambda){
  
  (lambda +    (1-lambda)*a)*dunif(x)
}
