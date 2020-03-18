#R code module for dealing with RNA Seq gene counts

# TPMK
# FPMK
# TPM

#check out the definition at 

#    These three metrics attempt to normalize for sequencing depth and gene length. Here’s how you do it for RPKM:
#
#    Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
#    Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
#    Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.

#    FPKM is very similar to RPKM. RPKM was made for single-end RNA-seq, where every read corresponded to a single fragment that was sequenced. FPKM was made for paired-end RNA-seq
#    TPM is very similar to RPKM and FPKM. The only difference is the order of operations. Here’s how you calculate TPM:

#    Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).#
#    Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#    Divide the RPK values by the “per million” scaling factor. This gives you TPM.

#    So you see, when calculating TPM, the only difference is that you normalize for gene length first, and then normalize for sequencing depth second. However, the effects of this difference are quite profound.
#    When you use TPM, the sum of all TPMs in each sample are the same. This makes it easier to compare the proportion of reads that mapped to a gene in each sample. In contrast, with RPKM and FPKM, the sum of the normalized reads in each sample may be different, and this makes it harder to compare samples directly.

#https://haroldpimentel.wordpress.com/2014/12/08/in-rna-seq-2-2-between-sample-normalization/
#https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/ 
#+++++++++++++++++++++

#'@title calculate the TPMs based on the raw counts
#'@description calculate the TPMs based on the raw counts and the effective lengths.
#'@details Take the inputs of counts and effective lengths and then calculate the TPMs. 
#' In this case, the effective lengths are the known input.
#'@param counts the data frame or matrix with a format of genes (rows) by samples (column)
#'@param effLen the data frame or matrix with a format identical to the counts data.
#'@return data frame or matrix
#'@seealso \code{\link{countsToTpm2}}  \code{\link{countsToFpkm}}
countToTpm <- function(counts, effLen)
{
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}

#'@title calculate the RPMKs/FPMKs based on the raw counts
#'@description calculate the RPMKs/FPMKs based on the raw counts and the effective lengths.
#'@details Take the inputs of counts and effective lengths and then calculate the TPMs. 
#' In this case, the effective lengths are the known input.
#'@param counts the data frame or matrix with a format of genes (rows) by samples (column)
#'@param effLen the data frame or matrix with a format identical to the counts data.
#'@return data frame or matrix
#'@seealso \code{\link{countsToTpm}} \code{\link{countsToTpm2}}
countsToFpkm <- function(counts, effLen)
{
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

#'@title convert RPMKs/FPMKs to TPMs
#'@description convert the RPMKs/FPMKs to the TPMs. The calculation is simple. 
#'@param fpmk the data frame or matrix with a format of genes (rows) by samples (column)
#'@return data frame or matrix
#'@seealso \code{\link{countsToTpm}} \code{\link{countsToTpm2}} \code{\link{countsToFpkm}}
#'@examples
#'################################################################################
#'# An example
#'################################################################################
#'cnts <- c(4250, 3300, 200, 1750, 50, 0)
#'lens <- c(900, 1020, 2000, 770, 3000, 1777)
#'countDf <- data.frame(count = cnts, length = lens)

#'# assume a mean(FLD) = 203.7
#'countDf$effLength <- countDf$length - 203.7 + 1
#'countDf$tpm <- with(countDf, countToTpm(count, effLength))
#'countDf$fpkm <- with(countDf, countToFpkm(count, effLength))
#'with(countDf, all.equal(tpm, fpkmToTpm(fpkm)))
#'countDf$effCounts <- with(countDf, countToEffCounts(count, length, effLength))
#'
fpkmToTpm <- function(fpkm)
{
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

#'@title calculate the effective counts 
#'@description the calculation is simple, counts* len/effLen, but not sure about what/why
#'		we need this. 
#'@param counts  the data frame or matrix with a format of genes (rows) by samples (column)
#'@param len  the data frame or matrix with a format of genes (rows) by samples (column) 
#'	for gene/transcript length. this should be by genes. identical for the same gene across samples 
#'@param effLen  the data frame or matrix with a format of genes (rows) by samples (column) 
#'	for effective gene/transcript length. this should be by genes. they could be different for the 
#'		same gene across different samples, since the samples might have different fragment lengths.
#'@return data frame or matrix
#'@seealso \code{\link{countsToTpm}} \code{\link{countsToTpm2}} \code{\link{countsToFpkm}}

countsToEffCounts <- function(counts, len, effLen)
{
    counts * (len / effLen)
}

#'@title calculate TPM based on the Raw counts (equal gene length for all samples)
#'@description Calculate the TPM based on the raw counts, the feature length and mean fragment length. 
#'	The feature length and mean fragment length are used to compute the effecture length. 
#'@details      TPM is very similar to RPKM and FPKM. The only difference is the order of operations. 
#'	Here’s how you calculate TPM:\cr
#'
#'    Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).\cr
#'    Count up all the RPK values in a sample and divide this number by 1,000,000. 
#'		This is your “per million” scaling factor.\cr
#'    Divide the RPK values by the “per million” scaling factor. This gives you TPM.\cr
#'
#'    So you see, when calculating TPM, the only difference is that you normalize for gene length first, 
#'	  and then normalize for sequencing depth second. However, the effects of this difference are quite profound.
#'    When you use TPM, the sum of all TPMs in each sample are the same. 
#'	  This makes it easier to compare the proportion of reads that mapped to a gene in each sample. 
#'	  In contrast, with RPKM and FPKM, the sum of the normalized reads in each sample may be different, 
#'	  and this makes it harder to compare samples directly.\cr
#'For effective length\cr
#'		(http://crazyhottommy.blogspot.com/2016/07/comparing-salmon-kalliso-and-star-htseq.html)\cr
#' I traced back to this paper by Lior pachter group Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation.
#'It is quite mathematical, but the general idea is:\cr
#'
#'    If we take the fragment length to be fixed, then the effective length is how many fragments can 
#'		occur in the transcript. This turns out to be length - frag_len +1. 
#'		The number of fragments coming from a transcript will be proportional to this number, 
#'		regardless of whether you sequenced one or both ends of the fragment. 
#'		In turn, when it comes to probabilistically assigning reads to transcripts the effective length 
#'		plays a similar role again. Thus for short transcripts, there can be quite a difference 
#'		between two fragment lengths. To go back to your example if you have transcript of length 310, 
#'		your effective length is 10 (if fragment length is 300) or 160 (if fragment length is 150) 
#'		in either case, which explains the discrepancy you see.
#'
#' From @Rob
#'
#'    The effective length is computed by using the fragment length distribution to 
#'		determine the effective number of positions that can be sampled on each transcript. 
#'		You can think of this as convolving the fragment length distribution with 
#'		the characteristic function (the function that simply takes a value of 1) 
#'		over the transcript. For example if we observe fragments of length 50 --- 1000, 
#'		a position more than 1000 bases from the end of the transcript will contribute 
#'		a value of 1 to the effective length, while a position 150 bases will contribute a value of F(150), 
#'      where F is the cumulative distribution function of the fragment length distribution.
#'		For single end data, where we can't learn an empirical FLD, we use a gaussian whose mean 
#'		and standard deviation can be set with --fldMean and --fldSD respectively.
#'
#'	From Harold Pimentel's post above. He is in Lior Pachter's group.\cr
#'
#'    Effective length refers to the number of possible start sites a feature 
#'	could have generated a fragment of that particular length. In practice, 
#'	the effective length is usually computed as: \cr
#'		l_e=l_i - U_FDL + 1
#'    where uFDL is the mean of the fragment length distribution 
#'      which was learned from the aligned read. If the abundance estimation method you’re using 
#'     incorporates sequence bias modeling (such as eXpress or Cufflinks), 
#'		the bias is often incorporated into the effective length by making the feature 
#'   shorter or longer depending on the effect of the bias.
#'@param counts maxtrix or data frame in a format of genes/transcripts as rows and samples as columns
#'@param featureLength vector of length of genes/transcripts. it should have the same length as transcripts/genes 
#'@param meanFragmentLength vector of mean fragment length. it should have the same length as samples. 
#'		By that, we could have different mean fragment length for different samples. But for each sample, 
#'		we assume all fragment has the same mean length.
#'
#'@return a data frame/matrix holding the TPM for genes and samples. 
 
countsToTpm2 <- function(counts, featureLength, meanFragmentLength) {

  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))

  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))

  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]

  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  return(tpm)
}

#'@title calculate TPM based on the Raw counts (unequal gene length across samples)
#'@description Calculate the TPM based on the raw counts, the feature length and mean fragment length. 
#'	The feature length and mean fragment length are used to compute the effecture length. 
#'@details check \code{\link{countsToTpm2}}
#'@param counts matrix of the gene counts 
#'@param featureLength matrix of the gene lengths
#'@param meanFragmentLength vector holding the fragment length 
#'
#'@export
#'
countsToTpm3 <- function(counts, featureLength, meanFragmentLength) {

	if(missing(counts)||missing(featureLength)||missing(meanFragmentLength))
	{
		stop("****RUN TIME ERROR***\n\t one of the inputs is not specified, please check\n\b")
	}
  # Ensure valid arguments.
  stopifnot(dim(featureLength) == dim(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))

  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength[,i] - meanFragmentLength[i] + 1
  }))

  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  #featureLength <- featureLength[idx,]

  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  return(tpm)
}
