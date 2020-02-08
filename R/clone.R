#R code to do clone related operation
#clone.R
#
#==defining 1)the functions to run sampling from clone assignment file
#==    2)run summary of clones.
#Feng @bu updated 2/27/2019

#'@title sampling to generate clone size distribution
#'@description generate the clone size distribution by bootstrap from input obtained 
#'	by the clone partation software.
#'@details be careful not to input the distribution size larger than the input files. 
#'	A warning will be issued when it happens.
#'@param filename string for the input file. Note this input file is the result by 
#'	the clonayst clone partition software. Usually it is named like "cloneAssignment.txt"
#'@param size integer for the size of the total number of members/Igs. This is like a 
#'	normalization which is identical to all sample so that we can compare between samples
#'@param times integer the number sampling to generate the distribution
#'@param sample.name string the sample name for the clone assignment.
#'
#'@return a data frame for the clone size distribution 
#'@examples
#' dir_file<-system.file("extdata", package="frLib")
#' ca<-read.table(file.path(dir_file,"CloneAssignments.txt"), header=T,sep="\t")
#' print(ca);
#'x<-generateCloneSizeBySampling(filename=file.path(dir_file,"CloneAssignments.txt"), size=2000, times=10, sample.name="sample13")
#'@export
#reading a file, named cloneAssignment, which is an output of 
#clonyast. each row of the file is for one sequence. This row shows
#that for this sequence which clones it belong to, etc.
#assuming a header and tab-delimited txt file.
#we basically will draw samples (size) and then 
#we will return a data.frame showing the clone size distribution
#like
#clone.size   #freq

#==========================
#'@import plyr

generateCloneSizeBySampling<-function(filename, size, times, sample.name)
{
	#clone.size<-data
	#first read the files 
	if(missing(filename))
	{
		stop("missing input file name, please specify....");
	}

	if(missing(size))
	{
		stop("missing input size (for sample size), please specify....");
	}
	if(missing(times))
	{
		stop("missing input times (repeat sampling times), please specify....");
	}
	if(missing(sample.name))
	{
		stop("missing sample name, please specify.......");
	}
	#filename="CloneAssignments.txt"
	clone.assignment<-read.table(file=filename, header=T, sep="\t");
	if(dim(clone.assignment)[1]<size)
	{
		stop("the specify the size is bigger than the total number of records, please check and try again")
	}
	
	num.records<-dim(clone.assignment)[1];
	
	clone.size<-NULL;
	for(i in 1:times)
	{
		index<-sample(c(1:num.records), size);
		clone.picked<-clone.assignment[index,];
		#now we get a sample clones, run stats
		
		clone.temp<-aggregate(clone.picked$CDR3Length, by=list(CloneID=clone.picked$CloneID), FUN=length)
		colnames(clone.temp)[2]<-"Members";
		
		clone.temp<-count(clone.temp, vars="Members")
		names(clone.temp)<-c("Members", paste0("freq_",i))
		if(i==1)
		{
			clone.size<-clone.temp;
		}else {
			clone.size<-merge(clone.size, clone.temp, all.x=T, all.y=T)
		}
	}
	#first get NA to be replaced by 0
	for(i in 1:dim(clone.size)[1])
	{
		clone.size[i,is.na(clone.size[i,])]<-0;
	}

	#now ran state for the samples
	clone.size<-cbind(clone.size[,1],apply(clone.size[,c(-1)], 1, mean))
	colnames(clone.size)<-c("Members", sample.name); 
	#done in here.
	clone.size;
}
