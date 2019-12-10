#' Microarray assay to understand acquired resistance to EGFR-targeted therapy (erlotinib) in cell lines
#' 
#' GSE38310: HCC827 cell lines, which are susceptible to erlotinib, were treated to induce resistence in vivo by repeated exposure. Gene expression of cell lines
#' treated with erlotibin was collected (alongside controls) in both susceptible and resistant cell lines. Expression captured using Illumina HumanHT-12 V3.0 expression beadchip.
#' 
#' The study authors cocluded that AXL kinase caused resistance to EGFR therapy. A clean peturbagenic experiment.
#' 
#' Generated 10 December 2019 from raw expression values dowloaded from GEO. Values were normalised using limma::neqc and processed with further limma functions.
#' 
#' @docType data
#' @usage data(erlotinib_ADRcellLines)
#' 
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38310}
#' @references \url{https://en.wikipedia.org/wiki/Erlotinib}
#' @references Zhang Z, Lee JC, Lin L, Olivas V et al. Activation of the AXL kinase causes resistance to EGFR-targeted therapy in lung cancer. Nat Genet 2012
#' @format A data.table with 48,803 rows, documenting probe-level differential expression three erlotinib cell line assays: against HCC827 (susceptible) and against ER3/T15.2 (evolved resistance). Addition gene symbols, entrez gene IDs and uniprot annotations also included.
"erlotinib_ADRcellLines"