#R wrapper for c code

#'@export
#'@useDynLib frLib add_rc
#
add_r_w<-function(x,y)
{
#here in this example, we using .C call, of which C function normally does not return values.
# and also pointer only as input. 
#but instead we using the input pointer to function as output. Inputs are passed as 
# list of pointers. so we get the output as the 3rd input.
	.C(add_rc, x, y, numeric(1))[[3]];
}

#
#====================
#'@title Rarefection based on incidence data (c code version, R wrapper)
#'@details rarefection_Inc is built on well known formula
#'              check chao and caldwell's work.
#' This one is the R wrapper to call c version for speed run. 
#'@param Q, incidence frequency count array, e.g. 
#'     c(1,250, 2,25,.....). Each pair is for one frequency count 
#'      in this case, there are 250 species with frequency of 1 in the data
#'@param T, is the total number of sampling unit/sample size. 
#'@param t, is the subsampling sizes. t>0 and t<T
#'@examples
#'require(iNEXT) #to use their data for testing.
#'data(ant)
#'iNEXT(ant$h2000m, q=0, datatype="incidence_freq", size=seq(1,300,20))
#'io<-iNEXT(ant$h500m, q=0, datatype="incidence_freq", size=seq(1,300,20))
#'
#'a<-ant$h2000m[-1]
#'aq<-aggregate(a, list(freq=a), length)
#'Qa<-rep(0, 2*dim(aq)[1])
#'Qa[seq(1, length(Qa), 2)]<-aq[,1]
#'Qa[seq(2, length(Qa), 2)]<-aq[,2]
#'getIndividuals(Q=Qa, T=200, seq(1,200,20))
#'x<-rarefection_Inc(Qa, T=200, seq(1,200,20))
#'y<-extrapolation_Inc(Qa, T=200,t=1:200, k=10)
#'plot(seq(1,200,20),x )
#'
#'b<-ant$h500m[-1]
#'bq<-aggregate(b, list(freq=b), length)
#'Qb<-rep(0, 2*dim(bq)[1])
#'Qb[seq(1, length(Qb), 2)]<-bq[,1]
#'Qb[seq(2, length(Qb), 2)]<-bq[,2]
#'getIndividuals(Q=Qb, T=230, seq(1,200,20))
#'x<-rarefection_Inc_w(Qb, T=230, seq(1,300,20))
#'y<-extrapolation_Inc(Qb, T=230,t=seq(1, 100,10), k=10)
#'@export
#'@useDynLib frLib rarefection_Inc_c
rarefection_Inc_w<-function(Q, T, t=c(1:T))
{
    #calling the 
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T);
    .Call(rarefection_Inc_c, as.integer(Qm[,1]), as.integer(Qm[,2]), 
                as.integer(dim(Qm)[1]), as.integer(T), as.integer(t), length(t));
}
#'@title Rarefection based on abundance data (R wrapper of c function)
#'@details rarefection_Inc is built on well known formula
#'              check chao and caldwell's work.
#'          NOTE: this is identical to rarefection_Inc, we only switched the parameter letters
#'          and still call the same functions inside. 
#' this is the c implemented code called by R. 
#'@param f, abunance frequency count array, e.g. 
#'     c(1,250, 2,25,.....). Each pair is for one frequency count 
#'      in this case, there are 250 species with frequency of 1 in the data
#'@param n, is the total number of individuals. 
#'@param m, is the subsampling size. m>0 and m<n
#'@examples
#'require(iNEXT) #to use their data for testing.
#'data(spider)
#'io<-iNEXT(spider$Girdled, q=0, datatype="abundance", size=seq(1,300,20))
#'
#'a<-spider$Girdled
#'aq<-aggregate(a, list(freq=a), length)
#'Qa<-rep(0, 2*dim(aq)[1])
#'Qa[seq(1, length(Qa), 2)]<-aq[,1]
#'Qa[seq(2, length(Qa), 2)]<-aq[,2]
#'x<-rarefection_Abu_w(Qa, n=sum(a), seq(1,200,20))
#'y<-extrapolation_Abu(Qa, n=sum(a),m=1:200, k=10)
#'plot(seq(1,200,20),x )
#'@export
rarefection_Abu_w<-function(f, n, m=c(1:n))
{
    #calling the 
    Qm<-matrix(f, nrow=length(f)/2, ncol=2, byrow=T)
    .Call(rarefection_Inc_c, as.integer(Qm[,1]), as.integer(Qm[,2]), 
                as.integer(dim(Qm)[1]), as.integer(n), as.integer(m), length(m));
}