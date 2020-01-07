#R code module for gene length
#the RNASeq data analysis needs information on the gene length. 
# For example, it needs the gene length to calculate FPMK/TPM.
#

#'@import GenomicRanges
#'@import rtracklayer
#'@import Rsamtools

#In this module, we pre-compute the lengths and save it 
#then each time, we simply load it.
#
#There are clearly different ways to "compute" the lengths, we actually use 
#the salmone way: we borrow it from the salmon quantification. Basically read
#in the quantification results by salmon and then only use the length information.
#The length was one field when we read/import using tximport.
#The other way is to calculate from the gtf file directly; check this
#at   https://bioinformatics.stackexchange.com/questions/2567/how-can-i-calculate-gene-length-for-rpkm-calculation-from-counts-data
#Here you can find some example R code to compute the gene length given a GTF file (
#	it computes GC content too, which you don't need). This uses one of 
#a number of ways of computing gene length, in this case 
#the length of the "union gene model". In this method, 
#the non-duplicated exons for each gene are simply 
#summed up ("non-duplicated" in that no genomic base is double counted). 
#This is a very simple way of getting a gene length.\cr

#There are alternative methods that you should be aware of, among which are:
#
#    Median transcript length: That is, the exonic lengths in 
#each transcript are summed and the median across transcripts is used. 
#This is probably a little more valid than the code that I linked to.
#
##    Per-sample effective gene lengths: the optimal method, 
#though it requires using something like RSEM, 
#which will give you an effective gene length.\cr
#    Using the length of the "major isoform" in your tissue of interest. 
#This isn't as good as method 2, but is more accurate than all of the others.

#At the end of the day, you're just coming up with a scale factor for each gene,
# so unless you intend to compare values across genes 
#(this is problematic to begin with) then it's questionable if using some of the more correct 
#but also more time-involved methods are really getting you anything.

#Edit: Note that if you want to plug these values into some sort of subtyping tool (TNBC in your case), 
#you should first start with some samples for which you know the subtype. 
#Then you can at least see if you're getting reasonable results. 
#After that, do read up on how the method works and 
#see if there's anything about RNAseq that makes it incompatible


#'@title obtain the gene lengths from gtf file
#'@description calculate the gene lengths based on the information in the gtf files.
#'@details  (to be added)
#'@return an integer vector holding the gene length. 
#'@examples
#'#DON'T RUN
#'#GTFfile<-"/home/feng/Windows/D/feng/LAB/MSI/mouse_genome_MM9/Homo_sapiens.GRCh38.91.gtf"
#'#ldata<-getGeneLengthGff(GTFfile, format="gtf", genome="GRCh38.91") 
#'@export
getGeneLengthGff<-function( gff.file,format="gtf", genome="GRCm38.71"
		
		)
{
	if(missing(gff.file))
	{
		stop("ERROR***: no gff.file specified");
	}
	#Load the annotation and reduce it
	cat("** ** Start reading the input file........\n"); flush.console();
	
	GTF <- import.gff(GTFfile, format=format, genome=genome, #asRangedData=F, 
				feature.type="exon")
	cat("\t.......Done\n")
	cat("** **Computing the lengths........\n"); flush.console();
	grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
	reducedGTF <- unlist(grl, use.names=T)
	elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))

	#Open the fasta file
	#FASTA <- FaFile(FASTAfile)
	#open(FASTA)

	#Add the GC numbers
	#elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
	elementMetadata(reducedGTF)$widths <- width(reducedGTF)
	cat("\t.......Done\n")
	cat("** **Generating the output........\n");flush.console();
	#Create a list of the ensembl_id/GC/length
	calc_GC_length <- function(x) {
		#nGCs = sum(elementMetadata(x)$nGCs)
		width = sum(elementMetadata(x)$widths)
		c(width)#, nGCs/width)
	}
	output <- (sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
#	colnames(output) <- c("Length")#, "GC")
	cat("\t........Done\n");
	return(output)
}
#'@import tximport
#'@import readr

#'@title get the gene lengths from salmon quantification
#'@description a wraper to call the tximport to get the gene lengths.
#'@param files file object containing the salmon quantification data.
#'@param tx2gene file mapping transcript to gene 
#'@return Gene length data
#'@examples 
#' qfile<-system.file("extdata", package="frLib")
#' qfile<-file.path(qfile, "quant.sf")
#' #GTFfile<-"/home/feng/Windows/D/feng/LAB/MSI/mouse_genome_MM9/Homo_sapiens.GRCh38.91.gtf"
#'# format="gtf"
#' #tx2gene	= getTx2Gene(GTFfile, format=format);
#' #lens<-getGeneLengthSalmon(qfile, tx2gen=tx2gene);
#' #
#'@export
getGeneLengthSalmon<-function( file, tx2gene
    )
{
	tx.salmon <- tximport(file, type = "salmon", tx2gene = tx2gene, 
                      #reader = read_tsv, 
					  countsFromAbundance = "no")
					  
	#now get the reads as list
	tx.salmon$length
}
#'@import rjson
#'@import GenomicFeatures

#'@import readr 

#'@export
getTx2Gene<-function ( gff.file, format="gtf"

		)
{
	cat("** **Reading input file and making transcript database........\n"); flush.console();
	txdb.ens <- makeTxDbFromGFF(gff.file,format=format) #from ensembl ftp

	k <- keys(txdb.ens, keytype = "GENEID")
	df <- select(txdb.ens, keys = k, keytype = "GENEID", columns = "TXNAME")
	cat("\t.......Done\n");
	cat("** **Selecting the appropriate fields for output.............\n"); flush.console();
	Tx.ensemble <- transcripts(txdb.ens,#edb,
			  columns = c("TXID", "GENEID", "TXNAME")#, "gene_name")#,
			  #return.type = "DataFrame"
			  )
	Tx.ens<-mcols(Tx.ensemble)
	Tx.ens<-as.data.frame(Tx.ens, stringsAsFactors=F)
	
	#nrow(Tx.ens)

	tx2gene<- Tx.ens[,c(3,2)]
	tx2gene[,2]<-as.character(tx2gene[,2])
	cat("\t........Done\n");
	#============above code used to build the tx2gene file. and also very weird, it seems we need to build and then save 
	#		and then read it. otherwise the tximport will run into errors??????

	#file.dir<-"E:\\feng\\LAB\\MSI\\singleCellRNASeq\\20170203_plateI_SCRS01\\salmon"
	#file.dir<-"/home/feng/Windows/D/feng/LAB/MSI/singleCellRNASeq/20170203_plateI_SCRS01/salmon/"
	#setwd(file.dir)
	#tx2gene<-read.table("tx2gene.txt", sep="\t", header=F)
	#colnames(tx2gene)<-c("tx_id", "gene_id")
	return (tx2gene)
}	
#'@title show available gene length data
#'@description show avaiable gene length data for loading
#'@details for users, this is the one on the system.dir
#'@param db.path the path to the data base.
#'@examples
#' showGeneLengths();
#'@export
showGeneLengths<-function(db.path=system.file("extdata", package="frLib")
		)
{
	cat("** **The following gene length data are precompiled in the package:\n");
	cat("\t\n")
	lfile<-file.path(db.path, "geneLengthData.RData")
	if(file.exists(lfile))
	{
	 load(lfile);
	 print(names(len.data));
	} else {
		cat("N/A\n");
	}
}
#'@title add the gene length data (administrative usage only)
#'@description  add the length data to the length data base file.
#'@details this is for the package maintainer usage only. we need to check the 
#'	source version of the data base and then update and install.
#'@param lns vector data 
#'@param ver version information should be unique as the key for the length data 
#'@param db.path !!!!should be the one in the source code!!!!
#'
#'@examples
#'#DON'T RUN
#'####first add salmon calculated gene lengths 
#'#GTFfile<-"/home/feng/Windows/D/feng/LAB/MSI/mouse_genome_MM9/Homo_sapiens.GRCh38.91.gtf"
#'# format="gtf"
#'#tx2gene	= getTx2Gene(GTFfile, format=format);
#'#lens<-getGeneLengthSalmon(qfile, tx2gen=tx2gene);
#'#db.path<-file.path("/home/feng/Feng/hg/frLib/inst/extdata/")
#'#addGeneLengths(lns=lens, ver="GRCh38.91.salmon",db.path=db.path);
#'#showGeneLengths()
#'#showGeneLengths(db.path);
#'#
#'#add GTF direct genelength data
#'#GTFfile<-"/home/feng/Windows/D/feng/LAB/MSI/mouse_genome_MM9/Homo_sapiens.GRCh38.91.gtf"
#'#lens<-getGeneLengthGff(GTFfile, format="gtf", genome="GRCh38.91") 
#'#addGeneLengths(lns=lens, ver="GRCh38.19.GTF",db.path=db.path);
#'#'#showGeneLengths()
#'#showGeneLengths(db.path);
#'
#'@export 
addGeneLengths<-function(lns, ver, db.path)
{
	if(missing(lns)||missing(ver)||missing(db.path))
	{
		stop ("please specify the length data ver b.path data");
	}
	len.data<-list();
	#check if there already a saved data in the folder
	#file_dir<-system.file("extdata", package="frLib")
	lfile<-file.path(db.path, "geneLengthData.RData");
	if(file.exists(lfile))
	{
		load(lfile);
	} 
	
	#now check whether the file has the ver, if yes,
	if(is.null(len.data[[ver]]))
	{
		cat("a new entry of lens data has been added to the data base\n");
	}
	else {
		cat("an existing entry has been udpated");
	}
	len.data[[ver]]<-lns;
	save(file=lfile, len.data);
}
#'@title load the pre-compiled gene length data 
#'@description load the specific version of precompiled gene length data from the package
#'@param ver either index or version of the gene length data.
#'@param db.path to the gene length data base 
#'@return gene length data
#'@examples  
#'	loadGeneLengths(1);
#'@export
loadGeneLengths<-function(ver=1, db.path=system.file("extdata", package="frLib"))
{
	
	dbfile<-file.path(db.path, "geneLengthData.RData")
	if(file.exists(dbfile))
	{
		load(dbfile);
		return(len.data[[ver]])
	} else {
		cat("ERROR****: data base file doesn't exist. NULL returned")
		return (NULL)
	}
	
}