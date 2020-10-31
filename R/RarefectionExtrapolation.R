#R code for doing Rarefection and Extrapolation
#based on Chao's work
#
#Feng 10/24/2020
#====================

#for doing the fitting.
#'@import minpack.lm

#'@title Rarefection based on incidence data
#'@details rarefection_Inc is built on well known formula
#'              check chao and caldwell's work. 
#'@param Qm, matrix with frequency and count as the columns
#'@param T, is the total number of sampling unit/sample size. 
#'@param t, is the subsampling size. t>0 and t<T
rarefection_Inc_each<-function(Qm, T, t)
{
        #build matrix
        Y<-Qm[,1] #<- all the possible frequency 
        numerator<-seq(0,t-1,1)
        denominator<-seq(0,t-1,1)
        
        #start computing.
        out<-rep(0, length(Y))
        for(i in 1:length(Y))
        {
            if((T-Y[i])<t){
                next;}
            n<-(T-Y[i])-numerator
            d<-T-denominator
            out[i]<-exp(sum(log(n/d)))
        }
        return(sum(Qm[,2])-sum(out*Qm[,2]))
}
#====================
#'@title Rarefection based on incidence data
#'@details rarefection_Inc is built on well known formula
#'              check chao and caldwell's work. 
#'@param Q, incidence frequency count array, e.g. 
#'     c(1,250, 2,25,.....). Each pair is for one frequency count 
#'      in this case, there are 250 species with frequency of 1 in the data
#'@param T, is the total number of sampling unit/sample size. 
#'@param t, is the subsampling size. t>0 and t<T
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
#'x<-rarefection_Inc(Qb, T=230, seq(1,300,20))
#'y<-extrapolation_Inc(Qb, T=230,t=seq(1, 100,10), k=10)
#'@export
rarefection_Inc<-function(Q,T,t=c(1:T))
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    out<-rep(0, length(t))
    for(i in 1:length(t))
    {
        out[i]=rarefection_Inc_each(Qm, T, t[i])
    }
    return (out);
}
#'@title Rarefection based on abundance data
#'@details rarefection_Inc is built on well known formula
#'              check chao and caldwell's work.
#'          NOTE: this is identical to rarefection_Inc, we only switched the parameter letters
#'          and still call the same functions inside. 
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
#'x<-rarefection_Abu(Qa, n=sum(a), seq(1,200,20))
#'y<-extrapolation_Abu(Qa, n=sum(a),m=1:200, k=10)
#'plot(seq(1,200,20),x )
#'@export
rarefection_Abu<-function(f,n,m=c(1:n))
{
    Qm<-matrix(f, nrow=length(f)/2, ncol=2, byrow=T)
    out<-rep(0, length(m))
    for(i in 1:length(m))
    {
        out[i]=rarefection_Inc_each(Qm, n, m[i])
    }
    return (out)
}
#'@title to calculate the individuals for doing rarefection
#'@param Q, incidence frequency count array, e.g. 
#'     c(1,250, 2,25,.....). Each pair is for one frequency count 
#'      in this case, there are 250 species with frequency of 1 in the data
#'@param T, is the total number of sampling unit/sample size. 
#'@param t, is the subsampling size. t>0 and t<T
#
#'@export
getIndividuals<-function(Q, T, t=c(1:T))
{
    #Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    out<-rep(0, length(t))
    N<-sum(Q[seq(1,length(Q),2)]*Q[seq(2,length(Q),2)])
    out<-t/T*N
    return(out)
}
#'@title to do the extrapolation_Inc on the basis of the incidence frequency count
#'@description  in this current version, we use Chao2 as the estimator of Q0_hat
#'@param Q, incidence frequency count array, e.g. 
#'     c(1,250, 2,25,.....). Each pair is for one frequency count 
#'      in this case, there are 250 species with frequency of 1 in the data
#'@param T, is the total number of sampling unit/sample size. 
#'@param t, is the subsampling size. t>0 and t<T
#'@examples 
#'#see the rarefection_Inc for examples 
#'@export
extrapolation_Inc<-function(Q, T, t,k=10)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    Sobs<-sum(Qm[,2]);
    Q0_hat<-as.numeric(Chao2(Q, T)[1]-Sobs)
    Q1<-as.numeric(getElementI(Qm, 1))
    #S_sample<-Sobs+Q0_hat*(1-exp(-1*t*Q1/(Q1+T*Q0_hat)))
    S_sample<-Sobs+Q0_hat*(1-(1-Q1/(Q1+T*Q0_hat))^t)
    return (S_sample);
}
#'@title to do the extrapolation_abundance on the basis of the abundance frequency count
#'@description  in this current version, we use Chao1 as the estimator of f0_hat.
#'          Note: this is almost identical to extrapolation_Inc except using Chao1
#'@param f, incidence frequency count array, e.g. 
#'     c(1,250, 2,25,.....). Each pair is for one abundance count 
#'      in this case, there are 250 species with frequency of 1 in the data
#'@param n, is the total number of individuals. 
#'@param m, is the subsampling size. m>0 and m<n
#'@examples 
#'#see the rarefection_Inc for examples 
#'@export
extrapolation_Abu<-function(f, n, m,k=10)
{
    Qm<-matrix(f, nrow=length(f)/2, ncol=2, byrow=T)
    Sobs<-sum(Qm[,2]);
    Q0_hat<-as.numeric(Chao1(f, n)[1]-Sobs)
    Q1<-as.numeric(getElementI(Qm, 1))
    #S_sample<-Sobs+Q0_hat*(1-exp(-1*t*Q1/(Q1+T*Q0_hat)))
    S_sample<-Sobs+Q0_hat*(1-(1-Q1/(Q1+n*Q0_hat))^m)
    return (S_sample);
}


##-==========the following is for doing collection model estimation of 
#'@title prepare input data for collection function modeling.
#'@description prepare input based on clone CDR3 data, clones identity and CDR3 length
#'                  check the ref Sobreon 1993 Accumulation function. 
#'                  note in ref Greiff 2017 they are not correct in two things: # of unique clones as independent 
#'                  variable and the log transformation of the independent variable.
#'@param clones a list of samples clone ids, each of  which is made up of v gene, j gene and CDR3
#'@param freqs a list of sample clone frequencies/abundances. This will be used to calculate 
#'          # of sequences for each sample. 
#'@param orders an array of indexes of the clones in the previous two variables.
#'@examples 
#'      # check the wl03R2/newPipeline/Repertoire_3_v1.0.r for an example.
#'@return a list of two inputs for modelling: total number of sequences and # of clones 
#'@export       
collection_data<-function(clones, freqs, orders=c(1:length(clones))) 
{
    #cr<-rep(0, length(orders))
    seq.num<-rep(0, length(orders))
    uClones.num<-rep(0, length(orders))
    uClones<-clones[[orders[1]]]
    ifreqs<-freqs[[orders[1]]]
    uClones.num[1]<-length(uClones)
    seq.num[1]<-sum(ifreqs)
    for(i in 2:length(orders))
    {
        seq.num[i]<-seq.num[i-1]+sum(freqs[[orders[i]]])
            uClones<-union(uClones,clones[[orders[i]]])
            uClones.num[i]<-length(uClones)
    }
    return (list(total=seq.num, unique=uClones.num));
}
#'@title run collection model fitting
#'@description fit the collection model
#'              ncs~ a(1-exp(-1*b*N))
#'@param ncs input array as the number of clones in the sample
#'@param N input array the number of total sequences in the sample
#'@param a initial value of a 
#'@param b intital value of b
#'@return a nls fitting object
#'@export 
collection_runFit<-function(ncs, N, init_a=3.9E9, init_b=1.7e-10)
{
    nl<-nlsLM(ncs~a*(1-exp(-1*b*N)), start=list(a=init_a,b=init_b), control=nls.lm.control(maxiter = 100))
    #m<-nls(cv~a*(1-exp(-1*b*N)), start=list(a=4e9,b=1.7e-10))
    return (nl)
}
#'@title collection model_fn3 
#'@description to calculate the expected clones based on # of samples
#'@param a the asymptotic total number of clones in the repertoire
#'@param b the slope for decreasing of detection probablity for a new clone.
#'@return the expected clones
#'@export
#'
model_fn3<-function (a, b, N)
{
    return (a*(1-exp(-1*b*N)))
}
#plot(getIndividuals(Q, T, 1:T),rarefection_Inc(Q, T, 1:T) )
#plot(1:T,rarefection_Inc(Q, T, 1:T) )
#t<-1:20000
#plot(t,extrapolation_Inc(Q, T,t, k=10))

####will do abundance
