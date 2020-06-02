

#' meta differential gene expression (DEG) analysis
#'
#' Given a tabular collection of gene sets and a table of
#'
#'
#'
setGeneric("metaDEG",
           valueClass = "data.frame",
           function(pValueSet, geneSet, betaUniformFit, ...) { standardGeneric("metaDEG") },
           signature = c("pValueSet", "geneSet", "betaUniformFit") )

# table and table
# The real method
#' @describeIn metaDEG Build a betaUniform model on the fly
#' @importFrom utils capture.output
setMethod("metaDEG",
          c(pValueSet="data.table", geneSet="data.table", betaUniformFit = "BetaUniformModel"),
          function(pValueSet, geneSet, betaUniformFit, pValAttr = "pVal", geneSetNameAttr = "geneSet", geneSetMembersAttr = "gene"){

            #TODO Check that geneSetMembersAttr is in pValueFrame

            # A little ugly, but we need to ensure that all the fields are the correct type (often comes up with entrez gene ID's)
            pValueSetDT <- unique(pValueSet[, .( as.character(.SD[[1]] ), as.numeric(.SD[[2]]) ),  .SDcols = c(geneSetMembersAttr, pValAttr)])
            setnames(pValueSetDT, c(geneSetMembersAttr, pValAttr) )
            geneSetDT <- unique(geneSet[, lapply(.SD, as.character), .SDcols = c(geneSetNameAttr, geneSetMembersAttr)])

            geneSetPvalsDT <- pValueSetDT[geneSetDT, , on = geneSetMembersAttr, allow.cartesian=TRUE, nomatch = NULL]

            if( any(duplicated(geneSetPvalsDT, by = c(geneSetNameAttr, geneSetMembersAttr))) ){

              duplicatedEntries <- geneSetPvalsDT[ duplicated(geneSetPvalsDT, by = c(geneSetNameAttr, geneSetMembersAttr) ) |
                                                   duplicated(geneSetPvalsDT, by = c(geneSetNameAttr, geneSetMembersAttr), fromLast = TRUE) ]

              examplesOfDuplication <- paste0( capture.output( head(duplicatedEntries,  n = 15)), collapse = "\n")

              warning("There are multiple ", geneSetMembersAttr ," entries for some ", geneSetNameAttr, " sets upon merging with the frame of P-values.",
                      " This is usually because of repeated entries in the pValueSet and can cause complications. This is usually the result of a many-to-one gene mapping.",
                      "Please inspect the input and aggregate (across splice-isoforms/microarray probes etc.) appropriately - or just choose the smallest p-value.",
                      " Here are some examples:\n",
                      examplesOfDuplication, "\n",
                      "(",nrow(duplicatedEntries)," such duplicates in total ...)")
            }

            if(any(geneSetPvalsDT[,any(is.na(get(pValAttr)))])){

              warning("Pvalue attribute ",pValAttr ," contans NA values. These will be filtered out, but it is better practice for it to be done elsewhere explicitly")
              geneSetPvalsDT <- geneSetPvalsDT[!is.na( get(pValAttr) )]
            }

            geneSetMetaDEGthPvals <- geneSetPvalsDT[ ,
                  .(pValueSet = list( sort(structure( .SD[, get(pValAttr)], names = .SD[, get(geneSetMembersAttr)]) ) ) ),
                  by = geneSetNameAttr]

            geneSetMetaDEGthPvals[,size := sapply(pValueSet, length)]
            geneSetMetaDEGthPvals[,members := sapply(pValueSet, names)]

            geneSetMetaDEGthPvals[, betaUniformMixtureP :=
                                      betaUniformPvalueSumTest(pValueSet[[1]], betaUniformFit),
                                    by = geneSetNameAttr]

            geneSetMetaDEGthPvals[, betaUniformMixtureQ := p.adjust(betaUniformMixtureP, method = "fdr")]

            return(geneSetMetaDEGthPvals[order(betaUniformMixtureP)])
          }
)


#' @describeIn metaDEG Build a betaUniform model on the fly
setMethod("metaDEG",
          c(pValueSet="data.table", geneSet="data.table", betaUniformFit = "missing"),
          function(pValueSet, geneSet, betaUniformFit, pValAttr, plot = TRUE, ...) {

            freshBetaUniformFit <- fitBetaUniformMixtureDistribution(unique(pValueSet)[,get(pValAttr)])

            if(plot){ print(plot(freshBetaUniformFit)) }

            callGeneric(pValueSet = pValueSet, geneSet = geneSet,  betaUniformFit = freshBetaUniformFit, pValAttr = pValAttr, ...) })



# named list (geneSets) and named vector - cast to below

# geneSetList (class) and named pValues (class)


