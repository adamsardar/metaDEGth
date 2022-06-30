#' Differential expression from microarray assay of ATF3 knockouts in HCT116
#' 
#' GSE74362 used recombinant AAV mediated genome engineering to knockdown the ATF3 transcription factor in a HCT116 cell line. These cells were
#' assayed using Illumina HumanHT-12 V4.0 expression beadchips and expression abundance measurements collected.
#' 
#' Generated 9 December 2019 from raw expression values dowloaded from GEO. Values were normalised using limma::neqc and processed with further limma functions.
#' 
#' It would appear that there is almost little-to-no differential expression in this study (inspect the  P-value distribution - almost all uniform). it is incldued
#' as a decoy set to test methods against.
#' @docType data
#' @usage data(koATF3_HCT116_diffexDT)
#' @examples
#' library(data.table)
#' library(ggplot2)
#' library(magrittr)
#' data(koATF3_HCT116_diffexDT)
#' 
#' koATF3_HCT116_diffexDT %>% ggplot(aes(x = koATF3_logFC, y = -log10(koATF3_pValue))) + 
#'     geom_point(aes(colour = p.adjust(koATF3_pValue, "fdr") < 0.05)) +
#'     scale_colour_manual(values = c("darkgrey","dodgerblue")) +
#'     theme_bw()
#' 
#' koATF3_geneIDExpress <- koATF3_HCT116_diffexDT[!is.na(geneSymbol), .SD[koATF3_pValue == min(koATF3_pValue, na.rm = TRUE)][1] , by = geneID]
#' 
#' koATF3_geneIDExpress[,qplot(koATF3_pValue, bins = 100)]
#' 
#' ATF3knockout_betaUniformModel <- fitBetaUniformMixtureDistribution(koATF3_geneIDExpress$koATF3_pValue, nStarts = 20)
#' # LLH fluctates a bit as we are on the bounds of parameter space!
#' 
#' ATF3knockout_betaUniformModel
#' 
#' noiseFractionUpperBound(ATF3knockout_betaUniformModel)
#' 
#' # If we were to use a P-value cutoff of 0.01, what would the FP rate be?
#' TP <- truePositiveFraction(ATF3knockout_betaUniformModel, pValueThreshold = 0.01)
#' FP <- falsePositiveFraction(ATF3knockout_betaUniformModel, pValueThreshold = 0.01)
#' round(100*FP/(TP+FP), digits = 1) # Over 30% would be false positives!
#' # Perform meta-analysis using the TTRUST TF->target sets
#' 
#' koATF3_geneIDExpress[, betaUnifScore_FDR0.05 := betaUniformScore(koATF3_pValue, ATF3knockout_betaUniformModel, FDR = 0.05)]
#' 
#' TTRUST_TF2targets_DT[, geneID := as.integer(geneID)]
#' 
#' TTRUST_TF_pvalues <- koATF3_geneIDExpress[ unique(TTRUST_TF2targets_DT[!is.na(geneID), .(TF,geneID)]), , on = "geneID"][!is.na(koATF3_pValue),.(pValueSet = list(koATF3_pValue), geneSet = list(geneSymbol), scoreSum = sum(betaUnifScore_FDR0.05)), by = TF]
#' 
#' 
#' TTRUST_TF_pvalues[, betaUniformMixtureP := betaUniformPvalueSumTest(pValueSet[[1]], ATF3knockout_betaUniformModel), by = TF]
#' TTRUST_TF_pvalues[, fishersP := fishersPvalueSumTest(pValueSet[[1]]), by = TF]
#' 
#' TTRUST_TF_pvalues[p.adjust(betaUniformMixtureP) < 0.05] # No significant gene sets
#' TTRUST_TF_pvalues[p.adjust(fishersP) < 0.05] # Quite a few ... too many for so little signal?
#' 
#' TTRUST_TF_pvalues[scoreSum > 0]
#' 
#' TTRUST_TF_pvalues[p.adjust(betaUniformMixtureP) < 0.05] # No significant gene sets
#' TTRUST_TF_pvalues[p.adjust(fishersP) < 0.05] # Quite a few ... too many for so little signal?
#' 
#' TTRUST_TF_pvalues %>% 
#'     ggplot(aes(x = -log(betaUniformMixtureP), y = -log(fishersP))) + 
#'     geom_point() +
#'     labs(title = "Plot to show how much more conservative beta-uniform P-Values are than Fisher's combined P-values",
#'          subtitle = "ATF3 knockout study (which has very little signal).\nBonferroni (conservative) multiple hypothesuis correction cut-off shown.") +
#'     geom_vline(xintercept = -log(0.05/800), linetype = "dashed", colour = "red") +
#'     geom_hline(yintercept = -log(0.05/800), linetype = "dashed", colour = "red")
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74362}
#' @references \url{https://www.genecards.org/cgi-bin/carddisp.pl?gene=ATF3}
#' @references Zhao J, Li X, Guo M, Yu J, Yan C. The common stress responsive transcription factor ATF3 binds genomic sites enriched with p300 and H3K27ac for transcriptional regulation. BMC Genomics. 2016
#' @references Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015
#' @format A data.table with 47300 rows, documenting differential gene expression between wildtype cells and knockdown (E11) cells, both treated with DMSO.
"koATF3_HCT116_diffexDT"