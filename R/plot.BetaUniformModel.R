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
      geom_line(stat = "function", fun = dbetauniform, args = list(a = x@a, lambda = x@lambda), aes(colour = "Beta-Uniform Mixture"), size=1.5) +
      geom_line(stat = "function", fun = uniformComponentBound, args = list(a = x@a, lambda = x@lambda),  aes(colour = "Upper Bound On\nUniform Component"), size=1.5) +
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

    eq <- TeX(paste0("$\\alpha =",round(x@a, digits = 2),
                    ", \\lambda =",round(x@lambda, digits = 2),
                    ", LLH =",round(-x@negLLH,digits = 2),"$", collapse = ""), output = "character")

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
