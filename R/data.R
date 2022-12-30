#package "frlib"
#
#data.R module.

#NOTE: how to save an R data object?
#	#save the data object at the R package module directory
#	usethis::use_data(data_name1, data_name2,...)
#	#then log the data object in the data.R file.

#first data object 
#' Example gene expression data analyzed by DEseq2 
#'
#' @source data set for gene expression analyzed by DESeq2
#'  Human data from PJ01
#' @format A data frame with columns:
#' \describe{
#'  \item{baseMean}{mean of the base across samples.}
#'  \item{log2FoldChange}{between groups}
#'  \item{lfcSE}{standard error}
#'  \item{stat}{negative binomial statq}
#'	\item{pvalue}{nominal p value}
#'	\item{padj}{FDR adjusted}
#' }
#' @examples
#' \dontrun{
#'  gexp_dds
#' }
"gexp_dds"