#' Microarray assay studying aberrant MAPK1/ERK signaling, causing resistance to EGFR kinase inhibitors
#' 
#' GSE37699: WZ4002 is a second generation TKI, effective at treating cells with the most common T790M mutation of aquired drug resistance. This study
#' inspects cell lines which have evolved resitance to WZ4002 by repeated exposure, concluding that it is a result of upregulated MAPK1 and/or ERK pathway
#' activity. Expression captured using Illumina HumanHT-12 V3.0 expression beadchip.
#' 
#' Generated 10 December 2019 from raw expression values dowloaded from GEO. Values were normalised using limma::neqc and processed with further limma functions.
#' 
#' @docType data
#' @usage data(WZ4002_adrH1975_diffexDT)
#' @examples
#' library(ggplot2)
#' data(WZ4002_adrH1975_diffexDT)
#' 
#' WZ4002_adrH1975_diffexDT %>%
#' ggplot(aes(x = H1975_logFC, y = -log10(H1975_pValue))) +
#'     geom_point(aes(colour = p.adjust(H1975_pValue) < 0.01 )) +
#'     scale_colour_manual(values = c("darkgrey","dodgerblue")) +
#'     theme_bw() 
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37699}
#' @references Ercan D, Xu C, Yanagita M, Monast CS, Pratilas CA, Montero J, et al. Reactivation of ERK signaling causes resistance to EGFR kinase inhibitors. Cancer Discov. 2012
#' @format A data.table with 22,277 rows, documenting probe-level differential expression between unchallenged and immune (evolved resistance) NCI-H1975 cell lines. Addition gene symbols and entrez gene IDs are also included.
"WZ4002_adrH1975_diffexDT"