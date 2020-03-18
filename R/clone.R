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
#' ca<-read.table(file.path(dir_file,"CloneAssignments.txt"), header=T,sep="\t") ; #note this CloneAssignments file has been truncated 
#' print(ca);
#'x<-generateCloneSizeBySampling(filename=file.path(dir_file,"CloneAssignments.txt"), size=20, times=10, sample.name="sample13")
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

#convert the numeric data frame or data matrix into groups/categories
#'@title categorize the data on the basis of the input criteria
#'@description We count the frequencies of a numeric data set based on the input category intervals 
#'	. 
#'@param dat the numeric data that will be grouped.
#'@param criteria numeric vector describing the grouping criteria or the intervals. The grouping is done
#'	in such way: 1)group less than the minimum value of criteria; 2) group larger than minimum one and 
#'		less than 2nd smallest value; ..... n-1) group larger than 2nd largest value and less the largest
#'		value; n) group with values larger than the maximum x value.
#'@param by the numeric data that we used to check against the criteria for categorizing
#'@return a data matrix holding the newly categorized data.
#'@examples
#'
#create random data
#'set.seed(0)

#'dat<-data.frame(g1=runif(100, 0,1), g2=runif(100, 0,1))
#'criteria<-c(0.1, 0.5,0.7)
#'
#'cdat<-categorize_matrix(as.matrix(dat), criteria,FUN="sum");
#'print(cdat)
#'cdat.2<-categorize_matrix(as.matrix(dat), by=dat[,1], criteria,FUN="sum");
#'print(cdat.2)
#'
#'@export 
categorize_matrix<-function(dat, criteria, by,FUN)
{
	if(missing(dat))
	{
		stop("please specify the input data \"dat\"!!")
	}
	if(missing(criteria))
	{
		stop("please specify the criteria variable!!")
	}
	criteria<-sort(criteria)

	if(length(criteria)<1)
	{
		stop("the input criteria variable is zero-lengthed!!")
	}
	
	if(class(dat)!="matrix")
	{
		stop("Expecting a data matrix as input")
	}
	#dat<-as.matrix(dat);
	re.dat<-matrix(0, nrow=length(criteria)+1, ncol=dim(dat)[2]);
	compare.against.criteria<-dat
	if(!missing(by))
	{
		compare.against.criteria<-by
	}
	for(j in 1:dim(dat)[2])
	{
		#for the first case
		#compare.values<-dat[,j]
		if(missing(by))
		{
			compare.against.criteria<-dat[,j]
		}
		index<-which(compare.against.criteria<=criteria[1]);
		re.dat[1,j]<-(eval(call(FUN,dat[index,j])));#dat[index,j],2,sum)
		if(length(criteria)>1)
		{
			#now, stop doing the job
			for(i in 2:length(criteria))
			{
				index<-which(compare.against.criteria>criteria[i-1]&compare.against.criteria<=criteria[i]);
				re.dat[i,j]<-eval(call(FUN,dat[index,j]));#(dat[index,j],2,sum)
			}
		}
		index<-which(compare.against.criteria>criteria[length(criteria)]);
		re.dat[length(criteria)+1,j]<-eval(call(FUN,dat[index,j]));#apply(dat[index,j],2,sum);
	}
	#cat("set the names!!")
	colnames(re.dat)<-colnames(dat);
	return (re.dat);
}

#convert the numeric data frame or data matrix into groups/categories
#'@title categorize the data on the basis of the input criteria
#'@description We count the frequencies of a numeric data set based on the input category intervals 
#'	. 
#'@param dat the data frame that will be grouped. There is only on column of data we will categorize.
#'@param criteria numeric vector describing the grouping criteria or the intervals. The grouping is done
#'	in such way: 1)group less than the minimum value of criteria; 2) group larger than minimum one and 
#'		less than 2nd smallest value; ..... n-1) group larger than 2nd largest value and less the largest
#'		value; n) group with values larger than the maximum x value.
#'@param index.data the index of the data column to be categorized
#'@param index.group the index of the grouping variable.
#'@param FUN string to indicate the function for summarize the data
#'@return a data matrix holding the newly categorized data and the grouping variable.
#'@examples
#'
#create random data
#'set.seed(0)

#'dat<-data.frame(g1=runif(200, 0,1), g2=rep(c("group1","group2"),100 ))
#'criteria<-c(0.1, 0.5,0.7)

#'cdat<-categorize_dataframe(dat, criteria,1,2, "length");

#'@export 
categorize_dataframe<-function(dat, criteria,index.data, index.group, FUN)
{
	if(missing(dat))
	{
		stop("please specify the input data \"dat\"!!")
	}
	if(missing(criteria))
	{
		stop("please specify the criteria variable!!")
	}
	criteria<-sort(criteria)
	
	if(missing(index.data))
	{
		stop("please specify the input index.data")
	}
	if(missing(index.group))
	{
		stop("please specify the input index.group")
	}

	if(length(criteria)<1)
	{
		stop("the input criteria variable is zero-lengthed!!")
	}
	if(missing(FUN))
	{
		stop("please specify the FUN")
	}
	groupNames<-unique(dat[,index.group])
	#dat<-as.matrix(dat);
	re.dat<-data.frame(rep(0,(length(criteria)+1)*length(groupNames)),
						rep(groupNames,each=(length(criteria)+1))
						)
	names(re.dat)<-names(dat)[c(index.data,index.group)]
	for(j in 1:length(groupNames))
	{	
			#for the first case
		dat.group=dat[dat[,index.group]==groupNames[j],index.data]
		index<-which(dat.group<=criteria[1]);
		re.dat[1+(length(criteria)+1)*(j-1),1]<-eval(call(FUN,dat.group[index]));#dat[index,j],2,sum)
		if(length(criteria)>1)
		{
			#now, stop doing the job
			for(i in 2:length(criteria))
			{
				index<-which(dat.group>criteria[i-1]&dat.group<=criteria[i]);
				re.dat[i+(length(criteria)+1)*(j-1),1]<-eval(call(FUN,dat.group[index]));#(dat[index,j],2,sum)
			}
		}
		index<-which(dat.group>criteria[length(criteria)]);
		re.dat[length(criteria)+1+(length(criteria)+1)*(j-1),1]<-eval(call(FUN,dat.group[index]));#apply(dat[index,j],2,sum);
		
	}
	return (re.dat);
}