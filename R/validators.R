#Uniprot regexp taken from http://www.uniprot.org/help/accession_numbers
uniprotRegexp <- "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"

validateSANTAresultDT <- function(candidateDT){
  
    assert_that(is.data.table(candidateDT),
                nrow(candidateDT) > 0,
                all( c("geneSet","geneSetMembers","SANTAobject","SANTAnodeScores","SANTA.Pvalue","SANTA.Qvalue") %in% colnames(candidateDT)),
                nrow(candidateDT[duplicated(geneSet)]) == 0)
  
  invisible(candidateDT)
}

on_failure(validateSANTAresultDT) <- function(call, env){ paste0(deparse(call$x), " should be a data.table resulting from running SANTAnetEnrichment function, with no duplciated geneSet names") }


############

validateAnnotationDT <- function(candidateDT){
  
  assert_that(is.data.table(candidateDT),
              "geneSet" %in% colnames(candidateDT),
              ncol(candidateDT) >= 2)
  
  invisible(candidateDT)
}

on_failure(validateAnnotationDT) <- function(call, env){ paste0(deparse(call$x), ": you've provided an (optional) annotation field - this must be a data.table with a geneSet column (used for join)") }

############

#' @importFrom igraph is.igraph
validateNetwork <- function(candidateNet, directed = FALSE){
  
  assert_that(is.igraph(candidateNet),
              is.directed(candidateNet) == directed,
              vcount(candidateNet) > 0)
  
  invisible(candidateNet)
}

on_failure(validateNetwork) <- function(call, env){ paste0(deparse(call$x), ":  expecting an igraph as input")}

############

validateGeneSetOfInterest <- function(candidateGeneSet, enrichDT){
  
  enrichDT %>% validateSANTAresultDT
  
  assert_that(is.character(candidateGeneSet),
              length(candidateGeneSet) > 0,
              all(candidateGeneSet %in% enrichSANTAdt$geneSet) )
  
  invisible(candidateGeneSet)
}

on_failure(validateGeneSetOfInterest) <- function(call, env){ paste0(deparse(call$x), ": expecting a vector of gene set names that are all in the supplied gene set data.table ...") }


############

validateGeneSetList <- function(candidateGeneSetList){
  
  assert_that(is.list(candidateGeneSetList),
              length(candidateGeneSetList) >= 1,
              !is.null(names(candidateGeneSetList)),
              !any(duplicated(names(candidateGeneSetList))),
              all(sapply(candidateGeneSetList,is.character)) )
  
  invisible(candidateGeneSetList)
}

on_failure(validateGeneSetList) <- function(call, env){ paste0(deparse(call$x), ": expecting a uniquely named list of gene sets") }


############


validateGeneSetOfInterest <- function(candidateGeneSet, enrichDT){
  
  enrichDT %>% validateSANTAresultDT
  
  assert_that(is.character(candidateGeneSet),
              length(candidateGeneSet) > 0,
              all(candidateGeneSet %in% enrichDT$geneSet) )
  
  invisible(candidateGeneSet)
}

on_failure(validateGeneSetOfInterest) <- function(call, env){ paste0(deparse(call$x), ": expecting a vector of gene set names that are all in the supplied gene set data.table ...") }

############

validateGeneScores <- function(candidateGeneScores){
  
  allNamed <- function(v){ !is.null(names(v))} #Little utility function for validation step
  
  assert_that(is.list(candidateGeneScores),
              all(sapply(candidateGeneScores,is.numeric)),
              all(sapply(candidateGeneScores,allNamed)),
              all(!is.null(names(candidateGeneScores))),
              !any(duplicated(names(candidateGeneScores))) )
  
  invisible(candidateGeneScores)
}

on_failure(validateGeneScores) <- function(call, env){ paste0(deparse(call$x), ": geneScores should be a list of named numeric values with (unduplicated) names") }

###############


validateBooleanFlag <- function(candidateBool){
  
  assert_that(is.logical(candidateBool),
              length(candidateBool) == 1)
  
  invisible(candidateBool)
}

on_failure(validateBooleanFlag) <- function(call, env){ paste0(deparse(call$x), " must be a single boolean value") }


############


validateSingleString <- function(candidateString){
  
  assert_that(is.character(candidateString),
              length(candidateString) == 1)
  
  invisible(candidateString)
}

on_failure(validateSingleString) <- function(call, env){ paste0(deparse(call$x), " must be a single character value") }

############

validateNumericCutoff <- function(candidateCutoff){
  
  assert_that(is.numeric(candidateCutoff),
              length(candidateCutoff) == 1)
  
  invisible(candidateCutoff)
}

on_failure(validateNumericCutoff) <- function(call, env){ paste0(deparse(call$x), " must be a single numeric value") }

############

validatePvalues <- function(candidatePvalues){
  
  assert_that( all( is.pval(candidatePvalues), na.rm = TRUE) )
  
  invisible(candidatePvalues)
}

on_failure(validatePvalues) <- function(call, env){ paste0(deparse(call$x), "  should be a P-value, hence must be between 0 and 1") }


############

validateIgraphWithPvalues <- function(candidateNetWithPvals){
  
  validateNetwork(candidateNetWithPvals)
  
  pValAttr <- str_extract(vertex_attr_names(candidateNetWithPvals),
                          regex("(\\bpval\\b|\\bpvalue\\b|\\bp\\.value\\b|\\bp\\.value\\b|\\bp\\b)", ignore_case = TRUE) ) %>% na.omit
  
  assert_that(length(pValAttr) == 1)
  
  V(candidateNetWithPvals)$PVALUE <- vertex_attr(candidateNetWithPvals, "pValue")
  
  invisible(candidateNetWithPvals)
}

on_failure(validateIgraphWithPvalues) <- function(call, env){ paste0(deparse(call$x), " must have a single P-value attribute. Found:", vertex_attr_names(call$x)  ) }

############

validateSingleInteger <- function(candidateInteger){
  
  assert_that(length(candidateInteger) == 1,
              is.numeric(candidateInteger), # We don't want to restrict ourselves to the integer class
              !is.infinite(candidateInteger),
              round(candidateInteger) == candidateInteger) # Is it a whole number
  
  invisible(candidateInteger)
}

on_failure(validateSingleInteger) <- function(call, env){ paste0(deparse(call$x), "  should be a single whole number") }


############

validateSingleInteger <- function(candidateInteger){
  
  assert_that(length(candidateInteger) == 1,
              is.numeric(candidateInteger), # We don't want to restrict ourselves to the integer class
              !is.infinite(candidateInteger),
              round(candidateInteger) == candidateInteger) # Is it a whole number
  
  invisible(candidateInteger)
}

on_failure(validateSingleInteger) <- function(call, env){ paste0(deparse(call$x), ":  should be a single whole number") }


############

validateSinglePositiveDefinite <- function(candidatePosDef){
  
  assert_that(is.numeric(candidatePosDef),
              length(candidatePosDef) == 1,
              candidatePosDef > 0)
  
  invisible(candidatePosDef)
}

on_failure(validateSinglePositiveDefinite) <- function(call, env){ paste0(deparse(call$x), ": should be a single positive definite number")}


############

#In particular for use with validator functions
is.pval <- function(x){ is.numeric(x) & x >=  0 & (x-1) <= 1E-12} # allow small epsilon