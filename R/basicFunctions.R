#R modules for functions to do some basic manipulations

#Feng @BU 2020-05

#'@title to reformat the data table to stack columns
#'@description reformat the data frame or matrix, so that the cognate columns are 
#'	stacked into one column and the other attribute columns are repeated accordingly
#'@param data the input data table, either a data frame or a matrix
#'@param index the vector containing all the cognate columns to be stacked. 
#'	any other columns will be the attibute/factors defining the data columns and repeated.
#'
#'@return a new data frame 
#'@examples
#'	x<-data.frame(v1 =rnorm(10), v2=rnorm(10,10), v3=rnorm(10,100), c1=rep(c("a","b"), 5), c2=rep(c("c","d"),each=5), c3=c(1:10))
#'	y<-stackTableData(x, index=c(1:3))
#'   print(x)
#'	print(y)
#'@export

stackTableData<-function( data, index=c(1:dim(data)[2])) 
{
	if(missing(data))
	{
		stop("please specify the data table...............")
	}
	if(class(data)!="data.frame" &&class(data)!="matrix")
	{
		stop("can not take data input other than data.frame or matrix")
	}
	if(length(index)<0||length(index)>dim(data))
	{
		stop("index out of range. please check")
	}
	index.repeat<-setdiff(c(1:dim(data)[2]), index)
	re<-data.frame(x=data[,index[1]])
	#now start doing the reformatting.
	if(length(index)>1)
		{
			for(i in 2:length(index))
			{
					re<-rbind(re, data.frame(x=data[,index[i]]))
			}
		}
	names(re)[1]=names(data)[index[1]]
		if(length(index.repeat)==0)
		return (re);
		#now take care the factors/attributes
		for(j in 1:length(index.repeat))
		{
			re<-cbind(re, rep(data[,index.repeat[j]], length(index)))
		}
		names(re)[-1]=names(data)[index.repeat]
		return (re)
}