% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cas9OCT4_hsZygote_diffexDT.R
\docType{data}
\name{cas9OCT4_hsZygote_diffexDT}
\alias{cas9OCT4_hsZygote_diffexDT}
\title{Differential expression from RNAseq assay of OCT4(POU5F1) CRISPR-Cas9 knockouts in human zygote}
\format{
A data.table with 15,028 rows, documenting differential gene expression (logFC and unadjusted P-values) of 15,003 ENSEMBL gene entities in zygotes grown with and wihout the OCT4 transcription factor gene. Addition gene symbols, entrez gene IDs and uniprot annotations also included (duplicated ENSG entries result from many-to-one mappings)
}
\source{
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100118}
}
\usage{
data(cas9OCT4_hsZygote_diffexDT)
}
\description{
GSE100118 is a CRISPR knockdown study of the pluripotency transcription factor OCT4 during human embryogenesis.

Generated 10 December 2019 from the trimmed reads associated with project on GEO. These were downloaded and aligned against the human genome using
Kallisto before transcript abundances were estimated using limma-voom. Finally, differential gene expression values were computed.
}
\examples{
library(ggplot2)
library(magrittr)
data(cas9OCT4_hsZygote_diffexDT)

cas9OCT4_hsZygote_diffexDT \%>\% ggplot(aes(x = koOCT4_logFC, y = -log10(koOCT4_pValue))) + 
    geom_point(aes(colour = p.adjust(koOCT4_pValue, "fdr") < 0.01)) +
    scale_colour_manual(values = c("darkgrey","dodgerblue")) +
    theme_bw()
    
   

cas9OCT4_hsZygote_diffexDT[,qplot(koOCT4_pValue, bins = 100)]

OCT4knockout_betaUniformModel <- fitBetaUniformMixtureDistribution(cas9OCT4_hsZygote_diffexDT$koOCT4_pValue, nStarts = 20)


OCT4knockout_betaUniformModel

noiseFractionUpperBound(OCT4knockout_betaUniformModel)

# If we were to use a P-value cutoff of 0.01, what would the FP rate be?
TP <- truePositiveFraction(OCT4knockout_betaUniformModel, pValueThreshold = 0.01)
FP <- falsePositiveFraction(OCT4knockout_betaUniformModel, pValueThreshold = 0.01)
round(100*FP/(TP+FP), digits = 1) # Only 2.8\% FPs


# Perform meta-analysis using the TTRUST TF->target sets    
   
cas9OCT4_hsZygote_diffexDT[, betaUnifScore_FDR0.05 := betaUniformScore(koOCT4_pValue, OCT4knockout_betaUniformModel, FDR = 0.05)]

TTRUST_TF2targets_DT[, geneID := as.character(geneID)]

TTRUST_TF_pvalues <- cas9OCT4_hsZygote_diffexDT[ unique(TTRUST_TF2targets_DT[!is.na(geneID), .(TF,geneID)]), , on = 'geneID'][!is.na(koOCT4_pValue),.(pValueSet = list(koOCT4_pValue), geneSet = list(geneSymbol), scoreSum = sum(betaUnifScore_FDR0.05)), by = TF]

TTRUST_TF_pvalues[, betaUniformMixtureP := betaUniformPvalueSumTest(pValueSet[[1]], OCT4knockout_betaUniformModel), by = TF]
TTRUST_TF_pvalues[, fishersP := fishersPvalueSumTest(pValueSet[[1]]), by = TF]

TTRUST_TF_pvalues[p.adjust(betaUniformMixtureP, method = 'fdr') < 0.05] # No significant gene sets
TTRUST_TF_pvalues[p.adjust(fishersP, method = 'fdr') < 0.05] # Quite a few ... too many for so little signal?

TTRUST_TF_pvalues[scoreSum > 0]

}
\references{
Fogarty NME, McCarthy A, Snijders KE, Powell BE et al. Genome editing reveals a role for OCT4 in human embryogenesis. Nature 2017

Pimentel H, Bray NL, Puente S, Melsted P, Pachter L. Differential analysis of RNA-seq incorporating quantification uncertainty. Nat Methods. 2017

Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 2014
}
\keyword{datasets}
