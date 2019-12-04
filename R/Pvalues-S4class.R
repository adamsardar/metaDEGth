

setClass("Pvalues",
         representation(.Data = "numeric"),
         validity = function(object){ return( all(is.pval(object@.Data)) ) })


setClass("Pvalue",
         contains = "ScalarNumeric",
         validity = function(object){ return( is.pval(object@.Data) ) })
