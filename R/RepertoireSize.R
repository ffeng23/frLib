#R module for doing repertoire size operations. We mainly follow 
# Chao's work, SpadeR and iNEXT
# done by Feng 10/21/2020
#--------------------------------
#'@title print the information/example about the module Repertoire operations
#'@description most of the information is borrowed from Chao's work in two packages
#'      , SpadeR and iNEXT. check the original work at github or chao's website at
#'      http://chao.stat.nthu.edu.tw/wordpress/
#'@examples 
#'  printRepertoireInfo()
#'@export
#'
printRepertoireInfo<-function()
{
    cat("
    This work is rewrite the code according to Chao's work, since their work runs slow and 
    somehow their warnings for doing too many things that we don't need. The information here is
    borrowed from their work at github and Chao's website.
    We do two things here: 
    1) ChaoSpecies (estimating species richness for one community).
    
    2) Diversity (estimating a continuous diversity profile and various diversity indices in 
    one community including species richness, Shannon diversity and Simpson diversity). 
    This function also features plots of empirical and estimated continuous diversity profiles. 
    
    Data Types

    It is very important to prepare your data in correct format. Data are generally classified 
    as abundance data and incidence data and there are five types of data input formats options 
    (datatype=\"abundance\", \"abundance_freq_count\", \"incidence_freq\", \"incidence_freq_count\",
    \"incidence_raw\"). In this module, we only do abundance_freq_count and incidence_freq_count.
    Type (1A) abundance-frequency counts data only for a single community 
    (datatype = \"abundance_freq_count\"): input data are arranged as (1 f1 2 f2 ... r fr)
    (each number needs to be separated by at least one blank space or separated by rows), 
    where r denotes the maximum frequency and fk denotes the number of species represented 
    by exactly k individuals/times in the sample. Here the data (f1, f2,..., fr) are referred to as 
    \"abundance-frequency counts\".
    Type (2A) incidence-frequency counts data only for a single community 
    (datatype=\"incidence_ freq_count\"): input data are arranged as (T 1 Q1 2 Q2 ... r Qr) 
    (each number needs to be separated by at least one blank space or separated by rows), 
    where Qk denotes the number of species that were detected in exactly k sampling units, 
    while r denotes the number of sampling units in which the most frequent species were found. 
    The first entry must be the total number of sampling units, T. The data (Q1,Q2,...,Qr) are referred to 
    as \"incidence frequency counts\".
    Note: for our data, we only use abundance_freq_count and incidence_freq_count, and importantly
    we don't have T in the numerical array. We have only c(f1, f2,....) or c(Q1, Q2.....) and input T as a sparate
    input. see examples in other function with \"?\"
    also we need to install SpadeR and iNEXT to run some examples and debugging.
    #---------------
    (3) DESCRIPTION OF ESTIMATORS/MODELS:\n

Chao2 (Chao, 1987): This approach uses the frequencies of uniques and duplicates to estimate the number of undetected species; see Chao (1987) or Eq. (11a) in Chao and Chiu (2016b).
     
Chao2-bc: A bias-corrected form for the Chao2 estimator; see Chao (2005).
  
iChao2: An improved Chao2 estimator; see Chiu et al. (2014).

ICE (Incidence-based Coverage Estimator): A non-parametric estimator originally proposed by Lee and Chao (1994) in the context of capture-recapture data analysis. The observed species are separated as frequent and infrequent species groups; only data in the infrequent group are used to estimate the number of undetected species. The estimated CV for species in the infrequent group characterizes the degree of heterogeneity among species incidence probabilities. See Eq. (12b) of Chao and Chiu (2016b), which is an improved version of Eq. (3.18) in Lee and Chao (1994). This model is also called Model(h) in capture-recapture literature where h denotes \"heterogeneity\".

ICE-1: A modified ICE for highly-heterogeneous cases.

    \n")
}
#'@title incidence-based Coverage Estimator
#'@description call printRepertoireInfo() to see more information and also check the SpadeR
#'      user guide for formulas. For incidence data only.
#'@details 
#'              the formulus is 
#'@param Q the array which is 
#'Q incidence freq count array. assuming it has the format as
#'      c(1, n1, 2, n2, 3, n3, 4, n4......., r, nr). it has a pair like (freq, count)
#'      in 1 (freq) unit, n1 (count) cases/species. n2 species in 2 sampling unit, etc
#'      assume it is sorted. note this is different from Chao's data, which has T as the 
#'      first element of the array. This could be either abundance_freq_count or
#'      incidence_freq_count.
#'@param k integer cutoff for rare species group.
#'@param T integer the total number of groups which at least contains one from species group.
#'          most of the time this is the total number of the groups.
#'@examples 
#'#benchmark for debugging
#' require(SpadeR)
#'data(ChaoSpeciesData)
#'k<-10
#'ice<-ICE(Q=ChaoSpeciesData$Inci_count[-1,1],k=12,T=ChaoSpeciesData$Inci_count[1,1])
#'ice_1<-ICE_1(Q=ChaoSpeciesData$Inci_count[-1,1],k=12,T=ChaoSpeciesData$Inci_count[1,1])
#'chao2<-Chao2(Q=ChaoSpeciesData$Inci_count[-1,1],T=ChaoSpeciesData$Inci_count[1,1])
#'ichao2<-iChao2(Q=ChaoSpeciesData$Inci_count[-1,1],T=ChaoSpeciesData$Inci_count[1,1])
#'ch<-ChaoSpecies(ChaoSpeciesData$Inci_count, "incidence_freq_count", k=10, conf=0.95)
#'
#'@export
ICE<-function(Q, k, T)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    #Qm<-Qm[Qm[,2]>0,]
    #names(Qm)<-Qm[,1]
    #T<-dim(Qm)[1]
    
    C<-coverage(Q, k, T)
    D_freq<-sum(Qm[(Qm[,1]>k),2])
    D_infreq<-sum(Qm[(Qm[,1]<=k),2])
    
    gm<-CV_rare(Q,k,T)
    Q1<-getElementI(Qm, 1)
    ice<-D_freq + D_infreq/C+Q1/C*gm
    
    #=====variance
    CV_infreq_h <- sqrt(gm)
    u <- c(1:k)  
x <-(Qm[Qm[,2]>0,])
     x<-rep(x[,1], x[,2])   
t<-T 
A<-0;
if (Qf(1, x) > 0 & Qf(2, x) > 0){
      A <- 2*Qf(2, x)/((t-1)*Qf(1, x) + 2*Qf(2, x))
    } else if (Qf(1, x) > 0 & Qf(2, x) == 0){
      A <- 2/((t-1)*(Qf(1, x) - 1) + 2)
    } else {
      A <- 1
    }
    C_infreq <- 1 - Qf(1, x)/sum(x[which(x <= k)])*(1-A)
    
D_infreq <- as.numeric(length(x[which(x <= k)]))
   CV_infreq <- CV_infreq_h
  diff <- function(q){
    if (CV_infreq_h != 0){
      n_infreq <- sum(x[which(x <= k)])
      si <- sum(sapply(u, function(u)u*(u-1)*Qf(u, x)))
      if ( q == 1){
        dc_infreq <-  - (n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x))*2*Qf(1, x)*(t - 1) - 
                           (t - 1)*Qf(1, x)^2*((t - 1)*(Qf(1, x) + n_infreq) + 2*Qf(2, x)))/(n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          t/(t - 1)*(C_infreq^2*n_infreq*(n_infreq - 1)*(D_infreq*si + Qf(1, x)*si) - #g3
                       Qf(1, x)*D_infreq*si*(2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*(n_infreq - 1) + C_infreq^2*n_infreq)             
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 - 
          (C_infreq - Qf(1, x)*dc_infreq)/C_infreq^2 #g4
      } else if (q == 2){
        dc_infreq <-  - ( - (t - 1)*Qf(1, x)^2*(2*(t - 1)*Qf(1, x) + 2*(n_infreq + 2*Qf(2, x))))/(n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          t/(t - 1)*(C_infreq^2*n_infreq*(n_infreq - 1)*Qf(1, x)*(si + 2*D_infreq) - Qf(1, x)*D_infreq*si*( #g3
            2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*2*(n_infreq - 1) + C_infreq^2*n_infreq*2)
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 - 
          ( - Qf(1, x)*dc_infreq)/C_infreq^2 #g4
      }else if(q > k){
        d <- 1
      } else {
        dc_infreq <-  - ( - (t - 1)*Qf(1, x)^2*((t - 1)*Qf(1, x)*q + 2*Qf(2, x)*q))/(n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          t/(t - 1)*(C_infreq^2*n_infreq*(n_infreq - 1)*Qf(1, x)*(si + q*(q - 1)*D_infreq) - Qf(1, x)*D_infreq*si*( #g3
            2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*q*(n_infreq - 1) + C_infreq^2*n_infreq*q)
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 - 
          ( - Qf(1, x)*dc_infreq)/C_infreq^2 #g4
      }
      return(d)
    }else{
      n_infreq <- sum(x[which(x <= k)])
      si <- sum(sapply(u, function(u)u*(u-1)*Qf(u, x)))
      if ( q == 1){
        dc_infreq <-  - (n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x))*2*Qf(1, x)*(t - 1) - 
                           (t - 1)*Qf(1, x)^2*((t - 1)*(Qf(1, x) + n_infreq) + 2*Qf(2, x)))/(n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      } else if (q == 2){
        dc_infreq <-  - ( - (t - 1)*Qf(1, x)^2*(2*(t - 1)*Qf(1, x) + 2*(n_infreq + 2*Qf(2, x))))/(n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      }else if(q > k){
        d <- 1
      } else {
        dc_infreq <-  - ( - (t - 1)*Qf(1, x)^2*((t - 1)*Qf(1, x)*q + 2*Qf(2, x)*q))/(n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      }
      return(d)
    }
  }
  
   i <- rep(sort(unique(x)),each = length(unique(x)))
  j <- rep(sort(unique(x)),length(unique(x)))       # all combination
  
  var_ice <- sum(mapply(function(i, j)diff(i)*diff(j)*COV_f(i, j, x, ice), i, j))
  if (var_ice > 0){
    var_ice <- var_ice
  } else {
    var_ice <- NA
    cat("Warning: In this case, it can't estimate the variance of Model(h) estimation", "\n\n")
  }
    return (c(ice=ice, var=var_ice))
}
#'@title incidence-based Coverage Estimator for highly highly-heterogeneous data
#'@description call printRepertoireInfo() to see more information and also check the SpadeR
#'      user guide for formulas. For incidence data only. 
#'@details 
#'              the formulus is 
#'@param Q the array which is 
#'Q incidence freq count array. assuming it has the format as
#'      c(1, n1, 2, n2, 3, n3, 4, n4......., r, nr). it has a pair like (freq, count)
#'      in 1 (freq) unit, n1 (count) cases/species. n2 species in 2 sampling unit, etc
#'      assume it is sorted. note this is different from Chao's data, which has T as the 
#'      first element of the array. This could be either abundance_freq_count or
#'      incidence_freq_count.
#'@param k integer cutoff for rare species group.
#'@param T integer the total number of groups which at least contains one from species group.
#'          most of the time this is the total number of the groups.
#'@examples 
#'#benchmark for debugging
#' require(SpadeR)
#'data(ChaoSpeciesData)
#'k<-10
#'ice<-ICE(Q=ChaoSpeciesData$Inci_count[-1,1],k=12,T=ChaoSpeciesData$Inci_count[1,1])
#'ice_1<-ICE_1(Q=ChaoSpeciesData$Inci_count[-1,1],k=12,T=ChaoSpeciesData$Inci_count[1,1])
#'chao2<-Chao2(Q=ChaoSpeciesData$Inci_count[-1,1],T=ChaoSpeciesData$Inci_count[1,1])
#'ichao2<-iChao2(Q=ChaoSpeciesData$Inci_count[-1,1],T=ChaoSpeciesData$Inci_count[1,1])
#'ch<-ChaoSpecies(ChaoSpeciesData$Inci_count, "incidence_freq_count", k=10, conf=0.95)
#'
#'@export
ICE_1<-function(Q, k, T)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    #Qm<-Qm[Qm[,2]>0,]
    #names(Qm)<-Qm[,1]
    #T<-dim(Qm)[1]
    
    C<-coverage(Q, k, T)
    D_freq<-sum(Qm[(Qm[,1]>k),2])
    D_infreq<-sum(Qm[(Qm[,1]<=k),2])
    
    gm<-CV_rare(Q,k,T)
    Q1<-getElementI(Qm, 1)
    ice<-ICE(Q,k,T)
    gm_1<-CV_rare_1(Q, k,T, ice[1])
    ice_1<-D_freq+D_infreq/C+Q1/C*gm_1
    
    ###=====variance
    t<-T 
    CV_infreq_h1 <- sqrt(gm_1)
    x <-(Qm[Qm[,2]>0,])
     x<-rep(x[,1], x[,2])   
     A<-0
     if (Qf(1, x) > 0 & Qf(2, x) > 0){
      A <- 2*Qf(2, x)/((t-1)*Qf(1, x) + 2*Qf(2, x))
    } else if (Qf(1, x) > 0 & Qf(2, x) == 0){
      A <- 2/((t-1)*(Qf(1, x) - 1) + 2)
    } else {
      A <- 1
    }
    C_infreq <- 1 - Qf(1, x)/sum(x[which(x <= k)])*(1-A)
    
    u <- c(1:k)   
D_freq <- as.numeric(length(x[which(x > k)]))
CV_infreq<-CV_infreq_h1    
  diff <- function(q){
    if (CV_infreq_h1 != 0){
      n_infreq <- sum(x[which(x <= k)])
      si <- sum(sapply(u, function(u)u*(u-1)*Qf(u, x)))
      if ( q == 1){
        dc_infreq <-  - (n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x))*2*Qf(1, x)*(t - 1) - 
                           (t - 1)*Qf(1, x)^2*((t - 1)*(Qf(1, x) + n_infreq) + 2*Qf(2, x)))/(n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          t/(t - 1)*(C_infreq^2*n_infreq*(n_infreq - 1)*(D_infreq*si + Qf(1, x)*si) - #g3
                       Qf(1, x)*D_infreq*si*(2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*(n_infreq - 1) + C_infreq^2*n_infreq)             
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 - 
          (C_infreq - Qf(1, x)*dc_infreq)/C_infreq^2 + #g4
          (t/(t - 1))^2*(C_infreq^3*n_infreq^2*(n_infreq - 1)^2*(2*Qf(1, x)*D_infreq*si^2 + Qf(1, x)^2*si^2) - #g5
                           Qf(1, x)^2*D_infreq*si^2*(3*C_infreq^2*dc_infreq*n_infreq^2*(n_infreq - 1)^2 + C_infreq^3*2*n_infreq*(n_infreq - 1)^2 + C_infreq^3*n_infreq^2*2*(n_infreq - 1))
          )/C_infreq^6/n_infreq^4/(n_infreq - 1)^4 - 
          (t/(t - 1))*si*(C_infreq^2*n_infreq*(n_infreq - 1)*2*Qf(1, x) - Qf(1, x)^2*(2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*(n_infreq - 1) + C_infreq^2*n_infreq) #g6
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 
      } else if (q == 2){
        dc_infreq <-  - ( - (t - 1)*Qf(1, x)^2*(2*(t - 1)*Qf(1, x) + 2*(n_infreq + 2*Qf(2, x))))/(n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          t/(t - 1)*(C_infreq^2*n_infreq*(n_infreq - 1)*Qf(1, x)*(si + 2*D_infreq) - Qf(1, x)*D_infreq*si*( #g3
            2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*2*(n_infreq - 1) + C_infreq^2*n_infreq*2)
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 - 
          ( - Qf(1, x)*dc_infreq)/C_infreq^2 + #g4
          (t/(t - 1))^2*Qf(1, x)^2*(C_infreq^3*n_infreq^2*(n_infreq - 1)^2*(si^2 + D_infreq*2*si*2) - #g5
                                     D_infreq*si^2*(3*C_infreq^2*dc_infreq*n_infreq^2*(n_infreq - 1)^2 + C_infreq^3*2*n_infreq*2*(n_infreq - 1)^2 + C_infreq^3*n_infreq^2*2*(n_infreq - 1)*2)
          )/C_infreq^6/n_infreq^4/(n_infreq - 1)^4 - 
          t/(t - 1)*Qf(1, x)^2*(C_infreq^2*n_infreq*(n_infreq - 1)*2 - si*(2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*2*(n_infreq - 1) + C_infreq^2*2*n_infreq)
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2
      }else if(q > k){
        d <- 1
      } else {
        dc_infreq <-  - ( - (t - 1)*Qf(1, x)^2*((t - 1)*Qf(1, x)*q + 2*Qf(2, x)*q))/(n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          t/(t - 1)*(C_infreq^2*n_infreq*(n_infreq - 1)*Qf(1, x)*(si + q*(q - 1)*D_infreq) - Qf(1, x)*D_infreq*si*( #g3
            2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*q*(n_infreq - 1) + C_infreq^2*n_infreq*q)
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 - 
          ( - Qf(1, x)*dc_infreq)/C_infreq^2 + #g4
          (t/(t - 1))^2*Qf(1, x)^2*(C_infreq^3*n_infreq^2*(n_infreq - 1)^2*(si^2 + D_infreq*2*si*q*(q - 1)) - #g5
                                     D_infreq*si^2*(3*C_infreq^2*dc_infreq*n_infreq^2*(n_infreq - 1)^2 + C_infreq^3*2*n_infreq*q*(n_infreq - 1)^2 + C_infreq^3*n_infreq^2*2*(n_infreq - 1)*q)
          )/C_infreq^6/n_infreq^4/(n_infreq - 1)^4 -
          t/(t - 1)*Qf(1, x)^2*(C_infreq^2*n_infreq*(n_infreq - 1)*q*(q - 1) - #g6
                                 si*(2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*q*(n_infreq - 1) + C_infreq^2*n_infreq*q)
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2                                 
      }
      return(d)
    }else{
      n_infreq <- sum(x[which(x <= k)])
      si <- sum(sapply(u, function(u)u*(u-1)*Qf(u, x)))
      if ( q == 1){
        dc_infreq <-  - (n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x))*2*Qf(1, x)*(t - 1) - 
                           (t - 1)*Qf(1, x)^2*((t - 1)*(Qf(1, x) + n_infreq) + 2*Qf(2, x)))/(n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      } else if (q == 2){
        dc_infreq <-  - ( - (t - 1)*Qf(1, x)^2*(2*(t - 1)*Qf(1, x) + 2*(n_infreq + 2*Qf(2, x))))/(n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      }else if(q > k){
        d <- 1
      } else {
        dc_infreq <-  - ( - (t - 1)*Qf(1, x)^2*((t - 1)*Qf(1, x)*q + 2*Qf(2, x)*q))/(n_infreq*((t - 1)*Qf(1, x) + 2*Qf(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      }
      return(d)
    }
  } 
  i <- rep(sort(unique(x)),each = length(unique(x)))
  j <- rep(sort(unique(x)),length(unique(x)))       # all combination
  
  var_ice1 <- sum(mapply(function(i, j)diff(i)*diff(j)*COV_f(i, j,x, ice_1), i, j))
  if (var_ice1 > 0){
    var_ice1 <- var_ice1
  } else {
    var_ice1 <- NA
    cat("Warning: In this case, it can't estimate the variance of Model(h)-1 estimation", "\n\n")
  }
    return(c(ice_1=ice_1, var=var_ice1))
}
#'@title Chao2
#'@description call printRepertoireInfo() to see more information and also check the SpadeR
#'      user guide for formulus 
#'@details 
#'              the formulus is 
#'@param Q the array which is 
#'Q incidence freq count array. assuming it has the format as
#'      c(1, n1, 2, n2, 3, n3, 4, n4......., r, nr). it has a pair like (freq, count)
#'      in 1 (freq) unit, n1 (count) cases/species. n2 species in 2 sampling unit, etc
#'      assume it is sorted. note this is different from Chao's data, which has T as the 
#'      first element of the array. This could be either abundance_freq_count or
#'      incidence_freq_count.
#'@param T integer the total number of groups which at least contains one from species group.
#'          most of the time this is the total number of the groups.
#'@examples 
#'#benchmark for debugging
#' require(SpadeR)
#'data(ChaoSpeciesData)
#'ice<-ICE(Q=ChaoSpeciesData$Inci_count[-1,1],k=12,T=ChaoSpeciesData$Inci_count[1,1])
#'ice_1<-ICE_1(Q=ChaoSpeciesData$Inci_count[-1,1],k=12,T=ChaoSpeciesData$Inci_count[1,1])
#'chao2<-Chao2(Q=ChaoSpeciesData$Inci_count[-1,1],T=ChaoSpeciesData$Inci_count[1,1])
#'ichao2<-iChao2(Q=ChaoSpeciesData$Inci_count[-1,1],T=ChaoSpeciesData$Inci_count[1,1])
#'ch<-ChaoSpecies(ChaoSpeciesData$Inci_count, "incidence_freq_count", k=10, conf=0.95)
#'
#'@export
Chao2<-function(Q, T)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    
    Q1<-getElementI(Qm, 1)
    Q2<-getElementI(Qm, 2)
    D<-sum(Qm[,2])
    S<-D+(T-1)*Q1*(Q1-1)/(2*T)
    var_Chao2=(T-1)/T*Q1*(Q1-1)/2+((T-1)/T)^2*Q1*(2*Q1-1)^2/4-((T-1)/T)^2*Q1^4/4/S
    if(Q2!=0)
    {
        S<-D+(T-1)*Q1^2/(2*T*Q2)
        var_Chao2 <- Q2*((T - 1)/T*(Q1/Q2)^2/2 + ((T - 1)/T)^2*(Q1/Q2)^3 + ((T - 1)/T)^2*(Q1/Q2)^4/4)
    }
    #cat("s:", S)
        X<-c(S,var_Chao2)
        names(X)<-c("Chao2", "var")
    return (X)
}
#'@title Chao1 for abundance_freq_count data. 
#'@description call printRepertoireInfo() to see more information and also check the SpadeR
#'      user guide for formulus 
#'@details 
#'              the formulus is 
#'@param Q the array which is 
#'Q incidence freq count array. assuming it has the format as
#'      c(1, n1, 2, n2, 3, n3, 4, n4......., r, nr). it has a pair like (freq, count)
#'      in 1 (freq) unit, n1 (count) cases/species. n2 species in 2 sampling unit, etc
#'      assume it is sorted. note this is different from Chao's data, which has T as the 
#'      first element of the array. This could be either abundance_freq_count or
#'      incidence_freq_count.
#'@param T integer the total number of samples (individuals). This is sum(abundance * counts)
#'@examples 
#'#benchmark for debugging
#' require(SpadeR)
#'data(ChaoSpeciesData)
#'k<-10
#'ace<-ACE(Q=ChaoSpeciesData$Abu_count[,1],k=10)
#'ace_1<-ACE_1(Q=ChaoSpeciesData$Abu_count[,1],k=10)
#'chao1<-Chao1(Q=ChaoSpeciesData$Abu_count[,1],T=sum(ChaoSpeciesData$Abu_count[seq(2, 50,2),1]*ChaoSpeciesData$Abu_count[seq(1, 50,2),1]))
#'ichao1<-iChao1(Q=ChaoSpeciesData$Abu_count[,1],T=sum(ChaoSpeciesData$Abu_count[seq(2, 50,2),1]*ChaoSpeciesData$Abu_count[seq(1, 50,2),1]))
#'ch<-ChaoSpecies(ChaoSpeciesData$Abu_count, "abundance_freq_count", k=10, conf=0.95)
#'
#'@export
Chao1<-function(Q, T)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    
    Q1<-getElementI(Qm, 1)
    Q2<-getElementI(Qm, 2)
    D<-sum(Qm[,2])
    S<-D+(T-1)*Q1*(Q1-1)/(2*T)
    V<-  (T-1)/T*(Q1/(Q1-1)) /2 + ((T-1)/T)^2 *Q1*(2*Q1-1)^2/4 - ((T-1)/T)^2 *(Q1)^4/4/S
    
    if(Q2!=0)
    {
        S<-D+(T-1)*Q1^2/(2*T*Q2)
        V<-Q2* ( (T-1)/T*(Q1/Q2)^2/2 + ((T-1)/T)^2 *(Q1/Q2)^3 + ((T-1)/T)^2 *(Q1/Q2)^4/4)
    }
    
    z <- -qnorm((1 - 0.95)/2)
    t <- S - D
    K <- exp(z*sqrt(log(1 + V/t^2)))
    CI_Chao1 <- c(D + t/K, D + t*K)
    #cat("CI:", CI_Chao1)
    X<-c(Chao1=S, var=V, low=CI_Chao1[1], high=CI_Chao1[2])
    return (X)
}
#'@title improved chao2
#'@description call printRepertoireInfo() to see more information and also check the SpadeR
#'      user guide for formulus 
#'@details 
#'              the formulus is 
#'@param Q the array which is 
#'Q incidence freq count array. assuming it has the format as
#'      c(1, n1, 2, n2, 3, n3, 4, n4......., r, nr). it has a pair like (freq, count)
#'      in 1 (freq) unit, n1 (count) cases/species. n2 species in 2 sampling unit, etc
#'      assume it is sorted. note this is different from Chao's data, which has T as the 
#'      first element of the array. This could be either abundance_freq_count or
#'      incidence_freq_count.
#'@param k integer cutoff for rare species group.
#'@param T integer the total number of groups which at least contains one from species group.
#'          most of the time this is the total number of the groups.
#'@examples 
#'#benchmark for debugging
#' require(SpadeR)
#'data(ChaoSpeciesData)
#'k<-10
#'ice<-ICE(Q=ChaoSpeciesData$Inci_count[-1,1],k=12,T=ChaoSpeciesData$Inci_count[1,1])
#'ice_1<-ICE_1(Q=ChaoSpeciesData$Inci_count[-1,1],k=12,T=ChaoSpeciesData$Inci_count[1,1])
#'chao2<-Chao2(Q=ChaoSpeciesData$Inci_count[-1,1],T=ChaoSpeciesData$Inci_count[1,1])
#'ichao2<-iChao2(Q=ChaoSpeciesData$Inci_count[-1,1],T=ChaoSpeciesData$Inci_count[1,1])
#'ch<-ChaoSpecies(ChaoSpeciesData$Inci_count, "incidence_freq_count", k=10, conf=0.95)
#'
#'@export
iChao2<-function(Q, T)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    
    Q1<-getElementI(Qm, 1)
    Q2<-getElementI(Qm, 2)
    Q3<-getElementI(Qm, 3)
    Q4<-getElementI(Qm, 4)
    s_chao2<-Chao2(Q,T)
    if(Q4==0)
    {
        s_ichao2<-s_chao2[1] + (T - 3)/T*Q3/4/(Q4 + 1)*max(Q1 - (T - 3)/(T - 1)*Q2*Q3/2, 0)  
        return (c(s_ichao2, s_chao2[2]))
    }
    S<-s_chao2[1]+(T-3)/4/T*Q3/(Q4)*max(Q1-(T-3)/2/(T-1)*Q2*Q3/Q4,0)
    V<-s_chao2[2]
    #for variance
    x <-(Qm[Qm[,2]>0,])
     x<-rep(x[,1], x[,2])
     q1 <- f(1, x); q2 <- f(2, x); q3 <- f(3, x); q4 <- f(4, x)
    
    diff <- function(q, x){ # fq
    t<-T
    if(q4 == 0) q4 = 1
    if (q1 > 0 & q2 != 0){
      if (q == 1){
        d <- (t - 1)/t*q1/q2 - (t - 3)/t*q3/4/q4
      } else if (q == 2){
        d <- (t - 1)/t*q1^2/2/q2^2 - (t - 3)^2/t/(t - 1)*q3^2/8/q4^2
      } else if (q == 3){
        d <- (t - 3)/t*q1/4/q4
      } else {
        d <- -(t - 3)/t*q1*q3/4/q4^2 + (t - 3)^2/t/(t - 1)*q2*q3^2/4/q4^3
      }
    } else if (q1 > 1 & q2 == 0){
      if (q == 1){
        d <- (t - 1)/t*(2*q1 - 1)/2/(q2 + 1) + (t - 3)/t*q3/4/q4
      } else if (q == 2){
        d <- -(t - 1)/t*q1*(q1 - 1)/2/(q2 + 1)^2
      } else if (q == 3){
        d <- (t - 3)/t*q1/4/q4
      } else {
        d <- -(t - 3)/t*q1*q3/4/q4^2
      }    
    } else {
      d=0
    }
    return(d)
  }
  
  ind <- 1:4
  i <- rep(sort(unique(ind)),each = length(unique(ind)))
  j <- rep(sort(unique(ind)),length(unique(ind)))       # all combination
   var_iChao2<-0
  #    if (q1 - q2*q3/2/q4 > 0 & q3 != 0){
  if (q1 - (T - 3)/(T - 1)*q2*q3/2/q4 > 0 | q1 - (T - 3)/(T - 1)*q2*q3/2 > 0){
    var_iChao2 <- sum(mapply(function(i, j)diff(i, x)*diff(j, x)*COV_f(i, j, x , S), i, j))
  } else {
    var_iChao2 <- var_Chao2
  }
  
  if (var_iChao2 > 0){
    var_iChao2 <- var_iChao2
  } else {
    var_iChao2 <- NA
  }
    return (c(ichao2=S, var=var_iChao2))
}
#'@title improved chao2 for abundance data. 
#'@description call printRepertoireInfo() to see more information and also check the SpadeR
#'      user guide for formulus 
#'@details
#'              the formulus is 
#'@param Q the array which is 
#'Q abundance freq count array. assuming it has the format as
#'      c(1, n1, 2, n2, 3, n3, 4, n4......., r, nr). it has a pair like (freq, count)
#'      in 1 (freq) unit, n1 (count) cases/species. n2 species in 2 sampling unit, etc
#'      assume it is sorted. note this is abundance_freq_count .
#'@param T integer the total number of samples/individuals, sum(abundance*count)
#'@examples 
#'#benchmark for debugging
#' require(SpadeR)
#'data(ChaoSpeciesData)
#'k<-10
#'ace<-ACE(Q=ChaoSpeciesData$Abu_count[,1],k=10)
#'ace_1<-ACE_1(Q=ChaoSpeciesData$Abu_count[,1],k=10)
#'chao1<-Chao1(Q=ChaoSpeciesData$Abu_count[,1],T=sum(ChaoSpeciesData$Abu_count[seq(2, 50,2),1]*ChaoSpeciesData$Abu_count[seq(1, 50,2),1]))
#'ichao1<-iChao1(Q=ChaoSpeciesData$Abu_count[,1],T=sum(ChaoSpeciesData$Abu_count[seq(2, 50,2),1]*ChaoSpeciesData$Abu_count[seq(1, 50,2),1]))
#'ch<-ChaoSpecies(ChaoSpeciesData$Abu_count, "abundance_freq_count", k=10, conf=0.95)
#'
#'@export
iChao1<-function(Q, T)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    
    Q1<-getElementI(Qm, 1)
    Q2<-getElementI(Qm, 2)
    Q3<-getElementI(Qm, 3)
    Q4<-getElementI(Qm, 4)
    S<-Chao1(Q,T)
    V<-0
    if(Q4==0)
    {
        S<-S[1]+(T-3)/4/T*Q3/(Q4+1)*max(Q1-(T-3)/2/(T-1)*Q2*Q3/(Q4+1),0)
        V<- (T - 1)/T*Q1*(Q1 - 1)/2 + 
        ((T - 1)/T)^2*Q1*(2*Q1 - 1)^2/4 - ((T - 1)/T)^2*Q1^4/4/S[1]
    } else { #f4 is good 
        S<-S[1]+(T-3)/4/T*Q3/(Q4)*max(Q1-(T-3)/2/(T-1)*Q2*Q3/(Q4),0)
         V<-Q2*((T - 1)/T*(Q1/Q2)^2/2 + 
                         ((T - 1)/T)^2*(Q1/Q2)^3 + ((T - 1 )/T)^2*(Q1/Q2)^4/4)
    }
    #get variance
     x <-(Qm[Qm[,2]>0,])
     x<-rep(x[,1], x[,2])
     f1 <- f(1, x); f2 <- f(2, x); f3 <- f(3, x); f4 <- f(4, x)
        n<-T
    diff <- function(q, x){ # fq
      f1 <- f(1, x); f2 <- f(2, x); f3 <- f(3, x); f4 <- f(4, x) 
      if(f4 == 0) f4 = 1 
      if (f1 > 0 & f2 != 0){
        if (q == 1){
          d <- (n - 1)/n*f1/f2 - f3/4/f4
        } else if (q == 2){
          d <- (n - 1)/n*f1^2/2/f2^2 - f3^2/8/f4^2
        } else if (q == 3){
          d <- f1/4/f4
        } else {
          d <- -f1*f3/4/f4^2 + f2*f3^2/4/f4^3
        }
      } else if (f1 > 1 & f2 == 0){
        if (q == 1){
          d <- (n - 1)/n*(2*f1 - 1)/2/(f2 + 1) + f3/4/f4
        } else if (q == 2){
          d <- -(n - 1)/n*f1*(f1 - 1)/2/(f2 + 1)^2
        } else if (q == 3){
          d <- f1/4/f4
        } else {
          d <- -f1*f3/4/f4^2
        }    
      } else {
        d=0
      }
      return(d)
    }
    
    xx <- 1:4
    i <- rep(sort(unique(xx)),each = length(unique(xx)))
    j <- rep(sort(unique(xx)),length(unique(xx)))       # all combination
    
   
    if (f1 - f2*f3/2/f4 > 0 & f3 != 0){
        var_iChao1 <- sum(mapply(function(i, j)diff(i, x)*diff(j, x)*COV_f(i, j, x, S), i, j))
    } else {
        var_iChao1 <- V
    }
    
    if (var_iChao1 > 0){
      var_iChao1 <- var_iChao1
    } else {
      var_iChao1 <- NA
    }
        
    
    X<-c(iChao1=S, var=var_iChao1)
    return (X)
}
#'@title abundance-based coverage estimator
#'@description call printRepertoireInfo() to see more information and also check the SpadeR
#'      user guide for formulas. For incidence data only. 
#'@details 
#'              the formulus is 
#'@param Q the array which is 
#'abundance freq count array. assuming it has the format as
#'      c(1, n1, 2, n2, 3, n3, 4, n4......., r, nr). it has a pair like (freq, count)
#'      in 1 (freq) unit, n1 (count) cases/species. n2 species in 2 sampling unit, etc
#'      assume it is sorted. note this is abundance_freq_count .
#'@param k integer cutoff for rare species group.
#'
#'@examples 
#'#benchmark for debugging
#' require(SpadeR)
#'data(ChaoSpeciesData)
#'k<-10
#'ace<-ACE(Q=ChaoSpeciesData$Abu_count[,1],k=10)
#'ace_1<-ACE_1(Q=ChaoSpeciesData$Abu_count[,1],k=10)
#'chao1<-Chao1(Q=ChaoSpeciesData$Abu_count[,1],T=sum(ChaoSpeciesData$Abu_count[seq(2, 50,2),1]*ChaoSpeciesData$Abu_count[seq(1, 50,2),1]))
#'ichao1<-iChao1(Q=ChaoSpeciesData$Abu_count[,1],T=sum(ChaoSpeciesData$Abu_count[seq(2, 50,2),1]*ChaoSpeciesData$Abu_count[seq(1, 50,2),1]))
#'ch<-ChaoSpecies(ChaoSpeciesData$Abu_count, "abundance_freq_count", k=10, conf=0.95)
#'
#'@export
ACE<-function(Q, k)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    #Qm<-Qm[Qm[,2]>0,]
    #names(Qm)<-Qm[,1]
    #T<-dim(Qm)[1]
    
    C<-coverage_rare_abu(Q, k)
    D_freq<-sum(Qm[(Qm[,1]>k),2])
    D_infreq<-sum(Qm[(Qm[,1]<=k),2])
    
    gm<-CV_rare_abu(Q,k)
    Q1<-getElementI(Qm, 1)
    ace<-D_freq + D_infreq/C+Q1/C*gm
    
    #####get variance, first turn the array into a abundance array
    x <-(Qm[Qm[,2]>0,])
    x<-rep(x[,1], x[,2])
    
      n_rare <- as.numeric(sum(x[which(x <= k)]))
      D_rare <- as.numeric(length(x[which(x <= k)]))
      
      diff<- function(q){
        u<-c(1:k)
      if (gm != 0){
        si <- sum(sapply(u, function(u)u*(u - 1)*f(u, x)))
        if ( q == 1){
          d <- (1 - f(1, x)/n_rare + D_rare*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g1 
            ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*(D_rare*si + f(1, x)*si) - 
               f(1, x)*D_rare*si*(-2*(1 - f(1, x)/n_rare)*(n_rare - f(1, x))/n_rare^2*n_rare*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*(2*n_rare - 1))
            )/(1 - f(1, x)/n_rare)^4/n_rare^2/(n_rare - 1)^2 - #g2
            (1 - f(1, x)/n_rare + f(1, x)*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g3
            
        } else if(q > k){
          d <- 1
        } else {
          d <- (1 - f(1, x)/n_rare - D_rare*q*f(1, x)/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g1
            ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*f(1, x)*(si + D_rare*q*(q - 1)) - 
               f(1, x)*D_rare*si*(2*(1 - f(1, x)/n_rare)*f(1, x)*q/n_rare^2*n_rare*(n_rare - 1) + 
                                    (1 - f(1, x)/n_rare)^2*q*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*n_rare*q)
            )/(1 - f(1, x)/n_rare)^4/(n_rare)^2/(n_rare - 1)^2 + #g2
            (q*(f(1, x))^2/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g3
        }
        return(d)
      } else {
        if ( q == 1){
          d <- (1 - f(1, x)/n_rare + D_rare*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g1 
        } else if(q > k){
          d <- 1
        } else {
          d <- (1 - f(1, x)/n_rare - D_rare*q*f(1, x)/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g1
        }
        return(d)  
      }
    }
    i <- rep(sort(unique(x)),each = length(unique(x)))
    j <- rep(sort(unique(x)),length(unique(x)))       # all combination
    
    var_ace <- sum(mapply(function(i, j)diff(i)*diff(j)*COV_f(i, j,x, ace), i, j))
    if(var_ace<0)
    {
        var_ace<-NA
    }
    X<-c(ace=ace, var=var_ace)
    return (X)
}
#'@title incidence-based Coverage Estimator for highly highly-heterogeneous data
#'@description call printRepertoireInfo() to see more information and also check the SpadeR
#'      user guide for formulas. For incidence data only. 
#'@details 
#'              the formulus is 
#'@param Q the array which is 
#'abundance freq count array. assuming it has the format as
#'      c(1, n1, 2, n2, 3, n3, 4, n4......., r, nr). it has a pair like (freq, count)
#'      in 1 (freq) unit, n1 (count) cases/species. n2 species in 2 sampling unit, etc
#'      assume it is sorted. note this is abundance_freq_count .
#'@param k integer cutoff for rare species group.
#'
#'@examples 
#'#benchmark for debugging
#' require(SpadeR)
#'data(ChaoSpeciesData)
#'k<-10
#'ace<-ACE(Q=ChaoSpeciesData$Abu_count[,1],k=10)
#'ace_1<-ACE_1(Q=ChaoSpeciesData$Abu_count[-1,1],k=10)
#'chao1<-Chao1(Q=ChaoSpeciesData$Abu_count[,1],T=sum(ChaoSpeciesData$Abu_count[seq(2, 50,2),1]*ChaoSpeciesData$Abu_count[seq(1, 50,2),1]))
#'ichao1<-iChao1(Q=ChaoSpeciesData$Abu_count[,1],T=sum(ChaoSpeciesData$Abu_count[seq(2, 50,2),1]*ChaoSpeciesData$Abu_count[seq(1, 50,2),1]))
#'ch<-ChaoSpecies(ChaoSpeciesData$Abu_count, "abundance_freq_count", k=10, conf=0.95)
#'
#'@export
ACE_1<-function(Q, k)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    #Qm<-Qm[Qm[,2]>0,]
    #names(Qm)<-Qm[,1]
    #T<-dim(Qm)[1]
    
    C<-coverage_rare_abu(Q, k) #
    D_freq<-sum(Qm[(Qm[,1]>k),2])
    D_infreq<-sum(Qm[(Qm[,1]<=k),2])
    
    gm<-CV_rare1_abu(Q,k)
    Q1<-getElementI(Qm, 1)
    ace<-D_freq + D_infreq/C+Q1/C*gm
    
    ####
     x <-(Qm[Qm[,2]>0,])
    x<-rep(x[,1], x[,2])
    
      n_rare <-as.numeric(sum(x[which(x <= k)]))
      D_rare <-as.numeric( length(x[which(x <= k)]))
      #u <- c(1:k)   
      #f <- function(i, data){length(data[which(data == i)])}  
      
    diff <- function(q){
       u <- c(1:k)
      if (gm != 0){
        si <- sum(sapply(u, function(u)u*(u-1)*f(u, x)))
        if ( q == 1){
          d <- (1 - f(1, x)/n_rare + D_rare*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g1 
            ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*(D_rare*si + f(1, x)*si) - 
               f(1, x)*D_rare*si*(-2*(1 - f(1, x)/n_rare)*(n_rare - f(1, x))/n_rare^2*n_rare*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*(2*n_rare - 1))
            )/(1 - f(1, x)/n_rare)^4/n_rare^2/(n_rare - 1)^2 - #g2
            (1 - f(1, x)/n_rare + f(1, x)*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g3
            ((1 - f(1, x)/n_rare)^3*(n_rare*(n_rare - 1))^2*(2*f(1, x)*D_rare*si^2 + f(1, x)^2*si^2) - #g4
               f(1, x)^2*D_rare*si^2*(3*(1 - f(1, x)/n_rare)^2*(f(1, x) - n_rare)/(n_rare)^2*(n_rare*(n_rare - 1))^2 + 
                                        (1 - f(1, x)/n_rare)^3*2*n_rare*(n_rare - 1)^2 + (1 - f(1, x)/n_rare)^3*n_rare^2*2*(n_rare - 1)) 
            )/(1 - f(1, x)/n_rare)^6/n_rare^4/(n_rare - 1)^4 - 
            ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*(2*f(1, x)*si) - #g5
               f(1, x)^2*si*(2*(1 - f(1, x)/n_rare)*(f(1, x) - n_rare)/n_rare^2*n_rare*(n_rare - 1) + 
                               (1 - f(1, x)/n_rare)^2*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*n_rare) 
            )/(1 - f(1, x)/n_rare)^4/n_rare^2/(n_rare - 1)^2
        } else if(q > k){
          d <- 1
        } else {
          d <- (1 - f(1, x)/n_rare - D_rare*q*f(1, x)/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g1
            ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*f(1, x)*(si + D_rare*q*(q - 1)) - 
               f(1, x)*D_rare*si*(2*(1 - f(1, x)/n_rare)*f(1, x)*q/n_rare^2*n_rare*(n_rare - 1) + 
                                    (1 - f(1, x)/n_rare)^2*q*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*n_rare*q)
            )/(1 - f(1, x)/n_rare)^4/(n_rare)^2/(n_rare - 1)^2 + #g2
            (q*(f(1, x))^2/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g3
            ((1 - f(1, x)/n_rare)^3*n_rare^2*(n_rare - 1)^2*f(1, x)^2*(si^2 + 2*D_rare*si*q*(q - 1)) - #g4
               f(1, x)^2*D_rare*si^2*(3*(1 - f(1, x)/n_rare)^2*(f(1, x)*q/n_rare^2)*(n_rare*(n_rare - 1))^2 + 
                                        2*(1 - f(1, x)/n_rare)^3*n_rare*q*(n_rare - 1)^2 + 2*(1 - f(1, x)/n_rare)^3*n_rare^2*(n_rare - 1)*q)   
            )/(1 - f(1, x)/n_rare)^6/(n_rare)^4/(n_rare - 1)^4 - 
            ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*f(1, x)^2*q*(q - 1) - #g5
               f(1, x)^2*si*(2*(1 - f(1, x)/n_rare)*f(1, x)*q/n_rare^2*n_rare*(n_rare - 1) + 
                               (1 - f(1, x)/n_rare)^2*q*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*n_rare*q)
            )/(1 - f(1, x)/n_rare)^4/(n_rare)^2/(n_rare - 1)^2
        }
        return(d)
      } else {
        si <- sum(sapply(u, function(u)u*(u-1)*f(u, x)))
        if ( q == 1){
          d <- (1 - f(1, x)/n_rare + D_rare*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g1 
        } else if(q > k){
          d <- 1
        } else {
          d <- (1 - f(1, x)/n_rare - D_rare*q*f(1, x)/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g1
        }
        return(d)
      }
    }
    
    i <- rep(sort(unique(x)),each = length(unique(x)))
    j <- rep(sort(unique(x)),length(unique(x)))       # all combination
   
    var_ace <- sum(mapply(function(i, j)diff(i)*diff(j)*COV_f(i, j, x,ace), i, j))
    if(var_ace<0)
    {
        var_ace=NA
    }
    X<-c(ace1=ace, var=var_ace)
    return (X)
}
#####helper functions##################
#helper function, copied from SpadeR with modification, since
#in SpadeR, the do the abundance data inside their functions 
f <- function(i, data){
    return(length(data[which(data == i)]))
    }  
Qf <- function(i, data){
    return(length(data[which(data == i)]))
    }  
######helper functions 

#helper function, copied from SpadeR. copy right belongs to them.
COV_f <- function(i,j, x, ace){
      if (i == j){
        cov.f <- f(i, x)*(1 - f(i, x)/ace)
      } else {
        cov.f <- -f(i, x)*f(j, x)/ace
      }     
      return(cov.f)
    }

#helper function used by other functions.
coverage<-function(Q, k, T)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    #Qm<-Qm[Qm[,2]>0,]
    #names(Qm)<-Qm[,1]
    #T<-dim(Qm)[1]
    Q1<-0
    Q2<-0
    Q1<-getElementI(Qm, 1)
    Q2<-getElementI(Qm, 2)
    x<-Q1*(T-1)*Q1
    y<-0
    #temp<-Qm[Qm[,1]<=k,]
    for(i in 1:k)
    {
        y<-y+i*getElementI(Qm, i)*((T-1)*Q1+2*Q2)
    }
    C<-x/y
    return(1-C)
}
#here function for abundance data 
coverage_rare_abu<-function(Q, k)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    #Qm<-Qm[Qm[,2]>0,]
    #names(Qm)<-Qm[,1]
    #T<-dim(Qm)[1]
    Q1<-0
    #Q2<-0
    Q1<-getElementI(Qm, 1)
    #Q2<-getElementI(Qm, 2)
    x<-Q1
    y<-0
    #temp<-Qm[Qm[,1]<=k,]
    for(i in 1:k)
    {
        y<-y+getElementI(Qm, i)*i
    }
    C<-x/y
    return(1-C)
}

#helper function 
#Qm is the incidence freq count,
# the first column is the element freq (the species showing in how many sampling units) and second one is the count.
# the thing is that some of element is not in the array.
# for example, 
#   1 250
#   2 150
#   4 29
#...
# there is is no i=3, which is zero.
#
#'@export   
getElementI<-function(Qm,i)
{
    Qi<-0
    if(length(Qm[Qm[,1]==i,2])==1)
        {Qi<-Qm[Qm[,1]==i,2]}
        
    return(Qi)
}
#helper function 
CV_rare<-function(Q, k, T)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    #Qm<-Qm[Qm[,2]>0,]
    #names(Qm)<-Qm[,1]
    #T<-dim(Qm)[1]
    
    C<-coverage(Q, k, T)
    #D_freq<-sum(Qm[(Qm[,1]>k),2])
    D_infreq<-sum(Qm[(Qm[,1]<=k),2])
    x<-0
    y<-0
    temp<-Qm[(Qm[,1]<=k),]
    for(i in 1:k)
    {
        x<-x+i*(i-1)*getElementI(Qm, i)
        y<-y+(i*getElementI(Qm, i))
    }
    gma<-D_infreq/C * T /(T-1) * x/y/(y-1) -1
    
    return(max(gma, 0))
}
#for abandance data
CV_rare_abu<-function(Q, k)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    #Qm<-Qm[Qm[,2]>0,]
    #names(Qm)<-Qm[,1]
    #T<-dim(Qm)[1]
    
    C<-coverage_rare_abu(Q, k)
    #D_freq<-sum(Qm[(Qm[,1]>k),2])
    D_infreq<-sum(Qm[(Qm[,1]<=k),2])
    x<-0
    y<-0
    #temp<-Qm[(Qm[,1]<=k),]
    for(i in 1:k)
    {
        x<-x+i*(i-1)*getElementI(Qm, i)
        y<-y+(i*getElementI(Qm, i))
    }
    gma<-D_infreq/C  * x/y/(y-1) -1
    
    return(max(gma, 0))
}
#helper function 
CV_rare_1<-function(Q, k, T, S_ice)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    #Qm<-Qm[Qm[,2]>0,]
    #names(Qm)<-Qm[,1]
    #T<-dim(Qm)[1]
    
    C<-coverage(Q, k, T)
    D_freq<-sum(Qm[(Qm[,1]>k),2])
    #D_infreq<-sum(Qm[(Qm[,1]<=k),2])
    x<-0
    y<-0
    temp<-Qm[(Qm[,1]<=k),]
    for(i in 1:k)
    {
        x<-x+i*(i-1)*getElementI(Qm, i)
        y<-y+(i*getElementI(Qm, i))
    }
    gma<-(S_ice-D_freq)* T /(T-1) * x/y/(y-1) -1
    
    return(max(gma, 0))
}
#helper function 
CV_rare1_abu<-function(Q, k)
{
    Qm<-matrix(Q, nrow=length(Q)/2, ncol=2, byrow=T)
    #Qm<-Qm[Qm[,2]>0,]
    #names(Qm)<-Qm[,1]
    #T<-dim(Qm)[1]
    
    C<-coverage_rare_abu(Q, k)
    #D_freq<-sum(Qm[(Qm[,1]>k),2])
    #D_infreq<-sum(Qm[(Qm[,1]<=k),2])
    x<-0
    y<-0
    #temp<-Qm[(Qm[,1]<=k),]
    for(i in 1:k)
    {
        x<-x+i*(i-1)*getElementI(Qm, i)
        y<-y+(i*getElementI(Qm, i))
    }
    gm<-CV_rare_abu(Q,k)
    gma<-(1+(1-C)/C* x/(y-1) )*gm
    
    return(max(gma, 0))
}

#####



#ch2_Ig<-Chao2(Q,T=14)
#ich2_Ig<-iChao2(Q, T=14)
#ice_Ig<-ICE(Q, k=12, T=14)
#ice1_Ig<-ICE_1(Q, k=12, T=14)