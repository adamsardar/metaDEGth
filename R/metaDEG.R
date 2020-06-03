#' meta differential gene expression (DEG) analysis
#'
#' Given a collection of gene sets and
#'
#' @return A frame of gene sets, their members, their sizes and their meta-analysis P-values
#' @param pValueSet P-values from experimental assay
#' @param geneSet Gene set collections
#' @param betaUniformFit A beta-uniform model of the P-value distribution [optional - can be estimated from the P-Values provided]
#' @param ... Additional parameters for methods
#' @export
setGeneric("metaDEG",
           valueClass = "data.frame",
           function(pValueSet, geneSet, betaUniformFit, ...) { standardGeneric("metaDEG") },
           signature = c("pValueSet", "geneSet", "betaUniformFit") )


#' @describeIn metaDEG Build a beta-uniform model on the fly
#' @importFrom utils capture.output
#' @param pValAttr The name of the P-value attribute in the pValueSet frame
#' @param geneSetNameAttr The name of the gene set identifier attribute in the geneSet frame
#' @param geneSetMembersAttr The name of the gene set members attribute in the geneSet frame
setMethod("metaDEG",
          c(pValueSet="data.frame", geneSet="data.frame", betaUniformFit = "BetaUniformModel"),
          function(pValueSet, geneSet, betaUniformFit, pValAttr = "pVal", geneSetNameAttr = "geneSet", geneSetMembersAttr = "gene"){

            # TODO Check that geneSetMembersAttr is in pValueFrame
            # TODO stopifnot checks

            setDT(pValueSet)
            setDT(geneSet)

            # A little ugly, but we need to ensure that all the fields are the correct type (often comes up with entrez gene ID's)
            pValueSetDT <- unique(pValueSet[, .( as.character(.SD[[1]] ), as.numeric(.SD[[2]]) ),  .SDcols = c(geneSetMembersAttr, pValAttr)])
            setnames(pValueSetDT, c(geneSetMembersAttr, pValAttr) )
            geneSetDT <- unique(geneSet[, lapply(.SD, as.character), .SDcols = c(geneSetNameAttr, geneSetMembersAttr)])

            geneSetPvalsDT <- pValueSetDT[geneSetDT, , on = geneSetMembersAttr, allow.cartesian=TRUE, nomatch = NULL]

            if( any(duplicated(geneSetPvalsDT, by = c(geneSetNameAttr, geneSetMembersAttr))) ){

              duplicatedEntries <- geneSetPvalsDT[ duplicated(geneSetPvalsDT, by = c(geneSetNameAttr, geneSetMembersAttr) ) |
                                                   duplicated(geneSetPvalsDT, by = c(geneSetNameAttr, geneSetMembersAttr), fromLast = TRUE) ]

              examplesOfDuplication <- paste0( capture.output( head(duplicatedEntries,  n = 15)), collapse = "\n")

              warning("There are multiple ", geneSetMembersAttr ," entries for some ", geneSetNameAttr, " sets upon merging with the frame of P-values, which is ambiguous.",
                      " This is usually because of repeated entries in the pValueSet, typically the result of a many-to-one gene mapping.",
                      " Please inspect the input and aggregate (across splice-isoforms/microarray probes etc.) appropriately - or just choose the smallest p-value per gene.",
                      " Here are some examples of the ",nrow(duplicatedEntries), " duplicated rows:\n",
                      examplesOfDuplication, "\n")
            }

            if(any(geneSetPvalsDT[,any(is.na(get(pValAttr)))])){

              warning("Pvalue attribute ",pValAttr ," contans NA values. These shall be filtered out, but it is better practice for it to be done elsewhere explicitly")
              geneSetPvalsDT <- geneSetPvalsDT[!is.na( get(pValAttr) )]
            }

            geneSetMetaDEGthPvals <- geneSetPvalsDT[ ,
                  .(pValueSet = list( sort(structure( .SD[, get(pValAttr)], names = .SD[, get(geneSetMembersAttr)]) ) ) ),
                  by = geneSetNameAttr]

            geneSetMetaDEGthPvals[,size := sapply(pValueSet, length)]
            geneSetMetaDEGthPvals[,members := list(list(names(pValueSet[[1]]))), by = geneSetNameAttr]

            geneSetMetaDEGthPvals[, betaUniformMixtureP :=
                                      betaUniformPvalueSumTest(pValueSet[[1]], betaUniformFit),
                                    by = geneSetNameAttr]

            geneSetMetaDEGthPvals[, betaUniformMixtureQ := p.adjust(betaUniformMixtureP, method = "fdr")]

            return(geneSetMetaDEGthPvals[order(betaUniformMixtureP)])
})


#' @describeIn metaDEG Build a betaUniform model on the fly
setMethod("metaDEG",
          c(pValueSet="data.table", geneSet="data.table", betaUniformFit = "missing"),
          function(pValueSet, geneSet, betaUniformFit, pValAttr, plot = TRUE, ...) {

            freshBetaUniformFit <- fitBetaUniformMixtureDistribution(unique(pValueSet)[,get(pValAttr)])

            if(plot){ print(plot(freshBetaUniformFit)) }

            callGeneric(pValueSet = pValueSet, geneSet = geneSet,  betaUniformFit = freshBetaUniformFit, pValAttr = pValAttr, ...) })


#' @describeIn metaDEG Provide geneSet as a named list of members
#' @importFrom purrr map_lgl map2
setMethod("metaDEG",
          c(pValueSet="ANY", geneSet="list", betaUniformFit = "ANY"),
          function(pValueSet, geneSet,  geneSetNameAttr = "geneSet", geneSetMembersAttr = "gene", ...) {

            tryCatch(
              stopifnot(
                !is.null(names(geneSet)),
                !duplicated(names(geneSet))),
                all(map_lgl(geneSet, is.vector)),
              error = function(e){stop("geneSet must be a uniquely named list of members (as a vector) to use this method.")} )

            geneSetDT <- map2(names(geneSet), geneSet,
                                ~{ data.table(geneSet = .x, gene = unique(as.character(.y)) ) }) %>% #Cast is for tricky ID's like entrez (integers)
                           rbindlist

            setnames(geneSetDT, c(geneSetNameAttr, geneSetMembersAttr))

            callGeneric(pValueSet = pValueSet, geneSet = geneSetDT, geneSetNameAttr = geneSetNameAttr, geneSetMembersAttr = geneSetMembersAttr, ...)
          })


#' @describeIn metaDEG Provide pValueSet as a uniquely named numeric vector of pvalues
setMethod("metaDEG",
          c(pValueSet="numeric", geneSet="ANY", betaUniformFit = "ANY"),
          function(pValueSet, geneSet,  pValAttr = "pVal", geneSetMembersAttr = "gene", ...) {

            tryCatch(
              stopifnot(
                !is.null(names(pValueSet)),
                !duplicated(names(pValueSet)),
                all(is.pval(pValueSet))
              ), error = function(e){stop("pValueSet must be a uniquely named vector of P-values (no NA's) to use this method.")} )

            pValueSetDT <- data.table(gene = names(pValueSet), pVal = pValueSet)

            setnames(pValueSetDT, c(geneSetMembersAttr, pValAttr))

            callGeneric(pValueSet = pValueSetDT, geneSet = geneSet, pValAttr = pValAttr, geneSetMembersAttr = geneSetMembersAttr, ...)
          })
