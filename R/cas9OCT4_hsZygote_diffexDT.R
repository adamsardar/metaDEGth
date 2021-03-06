#' Differential expression from RNAseq assay of OCT4(POU5F1) CRISPR-Cas9 knockouts in human zygote
#' 
#' GSE100118 is a CRISPR knockdown study of the pluripotency transcription factor OCT4 during human embryogenesis.
#'  
#' Generated 10 December 2019 from the trimmed reads associated with project on GEO. These were downloaded and aligned against the human genome using
#' Kallisto before transcript abundances were estimated using limma-voom. Finally, differential gene expression values were computed.
#' 
#' @docType data
#' @usage data(cas9OCT4_hsZygote_diffexDT)
#' 
#' @examples
#' library(ggplot2)
#' data(cas9OCT4_hsZygote_diffexDT)
#' 
#' cas9OCT4_hsZygote_diffexDT %>% ggplot(aes(x = koOCT4_logFC, y = -log10(koOCT4_pValue))) + 
#'     geom_point(aes(colour = p.adjust(koOCT4_pValue, "fdr") < 0.01)) +
#'     scale_colour_manual(values = c("darkgrey","dodgerblue")) +
#'     theme_bw()
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100118}
#' @references Fogarty NME, McCarthy A, Snijders KE, Powell BE et al. Genome editing reveals a role for OCT4 in human embryogenesis. Nature 2017
#' @references Pimentel H, Bray NL, Puente S, Melsted P, Pachter L. Differential analysis of RNA-seq incorporating quantification uncertainty. Nat Methods. 2017
#' @references Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 2014
#' @format A data.table with 15,028 rows, documenting differential gene expression (logFC and unadjusted P-values) of 15,003 ENSEMBL gene entities in zygotes grown with and wihout the OCT4 transcription factor gene. Addition gene symbols, entrez gene IDs and uniprot annotations also included (duplicated ENSG entries result from many-to-one mappings)
"cas9OCT4_hsZygote_diffexDT"