#IgSeqModule.R  
#	by Feng @ 6/9/2019

#this module contains the functions for doing IgSeq data analysis
#'@title get sorted groups
#'@description return a sorted groups based on the grouping fields.
#'@details with a input data frame, we sort each group by the field \"ord\".
#'		Only the sorted first number of records will be returned.
#'
#'@return the data frame containing only the top/bottom records sorted.
#'@param x the input data frame
#'@param groups the array of grouping fields in the data frame
#'@param ord the ordering field.
#'@param num number of records to be returned.
#'@param decreasing bool to specify the order of the sorting.
#'@examples
#'#making a data frame
#'  set.seed(1)
#'	group<-rep(c("A","B"),20)
#'	sample<-rep(c(1,2,3,4),10) 
#' x<-rnorm(40, 1,2)
#' x<-data.frame(group=group, sample=sample, x=x)
#' y<-getSortedGroups(x=x, groups=c("group", "sample"),ord="x", num=2)
#'
#'@export
#now sort the groups based on the clones size.
	getSortedGroups<-function(x, #dataframe input
						groups,  #the variable used to
						ord, #the ordering variable, only scalar expected
						num,  #num of first elements to return 
						decreasing=T  #sorting decreasing or increasing
						)
		{
			if(missing(x))
			{
				stop("please specify the input data frame x!!!")
			}
			if(missing(groups))
			{
				stop("please specify the input grouping fields--groups!!!")
			}
			if(missing(ord))
			{
				stop("please specify the input sorting field ord!!!")
			}
			if(missing(num))
			{
				stop("please specify the num of records to returned for each group of data num!!!")
			}
			if(class(x)!="data.frame"&&class(x)!="matrix")
			{
				stop("please specify the input x as data frame or matrix!!")
			}
			if(length(ord)!=1)
			{
				warning("more than one elements were specified for \"ord\", only the first one is used")
				ord<-ord[1];
			}
			if(class(num)!="numeric"&class(num)!="integer")
			{
				stop("please specify the input data num as an integer!!!")
			}
			#dtf<-df.clone; groups<-c("tissue","sampleName" );num<-10; decreasing=T; ord="X.Members"
			ret<-NULL;
			#build a sorting criterion data frame
			num.groups<-length(groups)
			#need to know the length of each groups
			#len.array<-rep(0, num.groups)
			#totalNum.permute<-1;
			#for(i in 1:num.groups)
			#{
			#	len.array[i]<-length(unique(dtf[,groups[i]]))
			#	if(len.array[i]!=0)
			#	{
			#		totalNum.permute<-totalNum.permute*len.array[i];
			#	}
			#}
			
			#build the grouping variable permutation
			permute.group<-unique(x[,groups])
			totalNum.permute<-dim(permute.group)[1];
			#for(i in 1:num.groups)
			#{
			#	#len.array[i]<-length(unique(dtf[,groups[i]]))
			#	if(len.array[i]!=0)
			#	{
			#		tmp<-unique(as.character(dtf[,groups[i]]));
			#		tmp<-rep(tmp,totalNum.permute/len.array[i])
			#		permute.group<-cbind(permute.group, tmp)
			#	}
			#}
			#colnames(permute.group)<-groups;
			
			#now we can try to get the groups, sorted
			#i<-1 ; j<-1
			for(j in 1:totalNum.permute)
			{
				tmp<-x
				for(i in 1:num.groups)
				{
					tmp<-tmp[tmp[,groups[i]]==permute.group[j,i],]
				}
				#now do sorting on
				tmp<-tmp[order(tmp[,ord],decreasing=decreasing)[1:num],]
				ret<-rbind(ret, tmp);
			}
			return(ret);
		} #end of 

#'@title accumulate sum 
#'@description obtaining the accumulate sum of a vector
#'@param x vector of numbers
#'@return a numeric sum
#'@examples
#'	x<-c(1:10)
#' accumulateSum(x)
#'@export 

###function definition
###used to obtain the accumulative sum of a col or row.
###input assuming a vector only
accumulateSum<-function(x)
{
	y<-x;
	if(length(x)<2)
	{
		return (y)
	}
	
	for(i in 2:length(x))
	{
		y[i]<-y[i]+y[i-1]
	}
	return(y);
}