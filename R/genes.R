#'@import biomaRt


#'@title add annototations to data
#'@description add external gene names and decriptions to the data 
#'	on the basis of the ensembl_gene_id. We also assume an opened 
#'	biomaRt database with correct database (attributes).
#'@param df data frame contains the ensembl_gene_id was the row names
#'@param mart the biomaRt object contains the annotations
#'@param attributes a vector of fields that will be passed to getBM in order to
#'		get the information
#'@param ... extra parameters. Currently we don't use any. 
#'
#'@return the data frame combining both original dataframe with the annotatoins. The
#'		first several fields are the annotations and the rest are from the original
#'		data frame.
#'
#'@examples
#'	#load the example data
#'	data(gexp_dds)
#'	
#'	#get data frame.
#'	df.gexp<-data.frame(gexp_dds)
#'	
#'	#now call the function
#'	library(biomaRt)
#'	ensembl <- useMart("ensembl", host="useast.ensembl.org")#call by host a different ensemble server
#'	ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
#'	df.gexp<-add_anns(df.gexp, ensembl)
#'
#'@export
add_anns <- function(df, mart, attributes=c("ensemble_gene_id", "hgnc_symbol", "description")

			,...)
{
    nm <- rownames(df)
    anns <- getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol",
"description"),
        filters = "ensembl_gene_id", values = nm, mart = mart)
    anns <- anns[match(nm, anns[, 1]), ]
    colnames(anns) <- c("ID", "Gene Symbol", "Gene Description")
    df <- cbind(anns, df[, 1:ncol(df)])
    rownames(df) <- nm
    df
}
 
