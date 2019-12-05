# Simple classes that enforce 

setClass("Pvalues",
         representation(.Data = "numeric"),
         validity = function(object){ return( all(is.pval( na.omit(object@.Data)) ) ) })

setClass("Pvalue",
         contains = "ScalarNumeric",
         validity = function(object){ return( is.pval( na.omit(object@.Data)) ) })
