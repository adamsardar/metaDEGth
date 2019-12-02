#' Log fold change and significance data from 59 siRNA knockdowns of known transcription factors
#' 
#' Cusanovich et al. assayed many siRNA knockdowns of trancription factors in a lymphoblastoid cell line. 
#' Each knockdown was performed in triplicate (except for SP1 and IRF5 which were in hexaplet). 18 controls (siRNAs
#' that were not known to bind to a gene product) were also collected. This constitutes an excellent testing and validation
#' set for active regulon detection.
#' 
#' Generated 20th November 2018 from the unnormalised data downloaded from GEO 
#' 
#' @docType data
#' @usage data(TF_siRNA_knockdowns_GSE50588_summary)
#' 
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50588}
#' @references Cusanovich DA, Pavlovic B, Pritchard JK, Gilad Y. The functional consequences of variation in transcription factor binding. PLoS Genet 2014
#' 
#' @format A data.table with 123 columns: 1 detailing probe identifiers (HumanHT12v4probeID), 4 more annotating gene identifiers and descriptions (RefSeq_ID, Entrez_Gene_ID, Symbol, Definition)
#' and then 2 columns for each of the 59 transcription factors knocked down: logFC and P-value, refering to limma estimates of each statistic.
#' 
"TF_siRNA_knockdowns_GSE50588_summary"