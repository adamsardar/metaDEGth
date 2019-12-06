#' A Manually Curated Dataset Of Human Transcription Factor - Target Relationships
#' 
#' The TRRUST (Transcriptional Regulatory Relationships Unraveled by Sentence-based Text mining) resource is a semi-automated database of 
#' transcription factor relationships. It has been augmented with ensembl gene identifiers.
#' 
#' TTRUST is distributed under a Creative Commons Attribution-ShareAlike license. You are free to distribute this data independently of the package license, in keeping with that of creative commons.
#' 
#' @docType data
#' @usage data(TTRUST_TF2targets_DT)
#' 
#' @source \url{https://www.grnpedia.org/trrust/}
#' @references \url{https://creativecommons.org/licenses/by-sa/4.0/}
#' @references Han H, Cho J-W, Lee S, Yun A, Kim H, Bae D, et al. TRRUST v2: an expanded reference database of human and mouse transcriptional regulatory interactions. Nucleic Acids Res. 2018
#' @format A data.table detailing 795 transcription factors and 2,492 targets. Some target entries are duplicated on account of many-to-one mappings to ensembl identifiers. 
"TTRUST_TF2targets_DT"