#' Differential expression from RNAseq assay of TAZ siRNA knockdowns in HFL-1
#' 
#' GSE73555 studied genes regulated by the transcription co-factor TAZ in a lung fibroblast cell line (HFL-1). The authors used siRNA
#' to knock down the TAZ gene (aka WWTR1); two different siRNAs for TAZ and a negative control (NTC) siRNA. It serves 
#' as a perfect 'spike-in' study for us to test gene set enrichment techniques against.
#' 
#' Generated 6 December 2019 from the trimmed reads associated with project. These were downloaded and aligned against the human genome using
#' Kallisto before transcript abundances were estimated using limma-voom. Finally, differential gene expression values were computed.
#' 
#' @docType data
#' @usage data(siTAZdiffex_HFL1_DT)
#' 
#' @example
#' 
#' library(ggplot2)
#' data(siTAZdiffex_HFL1_DT)
#' 
#' siTAZdiffex_HFL1_DT %>% ggplot(aes(x = siTAZ1_logFC, y = -log10(siTAZ1_pValue))) + 
#'     geom_point(aes(colour = p.adjust(siTAZ1_pValue, "fdr") < 0.01)) +
#'     scale_colour_manual(values = c("darkgrey","dodgerblue")) +
#'     theme_bw()
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73555}
#' @references \url{https://www.genecards.org/cgi-bin/carddisp.pl?gene=WWTR1}
#' @references Noguchi S, Saito A, Mikami Y, Urushiyama H, Horie M, Matsuzaki H, et al. TAZ contributes to pulmonary fibrosis by activating profibrotic functions of lung fibroblasts. Sci Rep. 2017
#' @references Pimentel H, Bray NL, Puente S, Melsted P, Pachter L. Differential analysis of RNA-seq incorporating quantification uncertainty. Nat Methods. 2017
#' @references Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 2014
#' @format A data.table with 14,501 rows, documenting differential gene expression (logFC and unadjusted P-values) of 14,475 ENSEMBL gene entities in two experimental conditions (siTAZ1 and siTAZ2). Addition gene symbols, entrez gene IDs and uniprot annotations also included (duplicated ENSG entries result from many-to-one mappings)
"siTAZdiffex_HFL1_DT"