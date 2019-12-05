# Simple classes that enforce P-value type and logic
setClass("Pvalues",
         contains = "numeric",
         representation(.Data = "numeric",
                        names = "character"),
         validity = function(object){ return( all(is.pval( na.omit(object@.Data)) ) ) })

setClass("Pvalue",
         contains = "ScalarNumeric",
         validity = function(object){ return( is.pval( na.omit(object@.Data)) ) })
