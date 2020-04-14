###In this module we define the composition-related functions. 
#-created 4/11/2020 Feng @BU

#'@import compositions 

#'@title build ilr/balance tree bases
#'@description obtain the ilr/balance tree bases/coefficients from the 
#'	sequential binary tree. The coefficients then can be used to 
#'	carry out the ilr transformation.
#'@details The input is the sequential binary tree specified as a numeric 
#'	data matrix. The tree was construct to specify how to subcomp the node.
#'	For each row we will only accept three numbers, -1, 0 and 1. -1, means components
#'	are included the trees that will be on the denominator; 1 means ones are  on the 
#'	numerator; 0 means no inclusion in either tree. The formular to calculate  
#'	the balace is \cr
#'	b_i=sqrt(n_p*n_n/(n_p+n_n))log((cp_1*cp_2*...*cp_i)^(1/n_p)/(cn_1*cn_2....*cn_j)^(1/n_n))
#'	\cr 
#'	n_p and n_n are the number of elements with positive sign and negative sign, respectively.
#'	cp_i and cn_j are the compositional proprotions of the entries. 
#'	The bases/coefficients are calcuated as\cr 
#'		a_+=1/n_p *sqrt(n_p*n_n/(n_p+n_n))
#'		\cr
#'		a_-=-1/n_n*sqrt(n_p*n_n/(n_p+n_n))
#'		\cr
#'   To build the sb tree, we have an example as follows \cr
#'	\itemize{
#'		\item {"name"}{ "A","B", "C", "D", "E", "F"}
#'		\item {"1"}                     {1,   -1,     -1,    -1,    -1,    -1}
#'		\item {"2"}                     {0,    1,       1,     -1,    -1,    -1}
#'		\item {"3"}                     {0,    1,      -1,      0,     0,      0}
#'		\item {"4"}                     {0,    0,        0,     1,    -1,    -1}
#'		\item {"5"}                     {0,    0,        0,     0,      1,    -1} 
#'	}
#'  In row one, we want compare A with (B, C, D, E, F);\cr
#'  In row two, we compare B with C;  \cr
#'  In row three, we compare D with E and F;\cr
#'  In row four, we compare E with F
#' See examples below.
#'@param sbtree data matrix holding the sequential binary tree. Remember to use only -1, 1 and 1
#'	to represent the tree.
#'@return a data matrix has the columns for the transformation. Each column holds the bases for one 
#'	row in the input tree. (transposed)
#'@examples
#'	#build the sb tree, describe in the description above
#'bmatrix<-matrix(0, nrow=5, ncol=6)
#'bmatrix[1,]=c(1, rep(-1, 5))
#'bmatrix[2,]=c(0, rep(1, 2), rep(-1,3))
#'bmatrix[3,]=c(0, 1,-1, rep(0,3))
#'bmatrix[4,]=c(rep(0,3), 1,-1,-1)
#'bmatrix[5,]=c(rep(0,4), 1,-1)
#'print(bmatrix)
#' x<-buildBalanceBase(bmatrix)
#'ilr(c(1,2,3,4,5,6), V=x)
#'
#'#build a D=3 basic tree for demo
#'	ilrBase(D=3)
#'	x<-buildBalanceBase(bmatrix[c(4,5), c(4,5,6)])
#'	ilr(c(1,2,3, V=ilrBase(D=3)))
#'	ilr(c(1,2,3), V=x)
#'@export
buildBalanceBase<-function(sbtree)
{
	if(missing(sbtree))
	{
		stop("Please specify the sequential binary tree in matrix format as input\n");
	}
	if(class(sbtree)!="matrix")
	{
		stop("Please specify the sequential binary tree in matrix format as input\n");
	}
	
	nr<-dim(sbtree)[1]
	#nc-dim(sbtree)[2]
	bb<-sbtree
	#for each row, we calculate the coefficient 
	for(i in 1:nr)
	{
#	cat("loop ", i, "\n");
		#get number of positive and negative
		p <- sum(sbtree[i,]>0);
		n <- sum(sbtree[i,]<0);
#		cat("p:", p, "; n:",n, "\n")
		a_p<-sqrt(p*n/(p+n))/p
		a_n<-sqrt(p*n/(p+n))/n
		bb[i, bb[i,]>0]<-bb[i, bb[i,]>0]*a_p
		bb[i, bb[i,]<0]<-bb[i, bb[i,]<0]*a_n
	}
	return(t(bb))
}
