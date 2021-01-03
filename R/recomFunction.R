#R modules to take care of recombination event analysis 
#
#======12/19/2020
#So far we only have function to do deletion profiles (v genes and j gene.)
#more to come 
#---------------------\\
#'@title residual function for fitting the deletion profile (for V genes)
#'@description 
#'for v genes 
  #'this is the residual function used by by nls.lm fitting the deletion profile.
#'@param pm is the log transformed pwm matrix, is the one to be fitted 
#'@param observed is the probability for each observed CDR3 (6 nt long sub sequence) .
#'           vector of numberics
#'@param xx is the subsequence (6nt long) as the input (must be characters)
#'           xx is a data.frame with nrow=length(observed) and ncol=6
#'----not used for now vvvvvv
#' #@param cfs is the correction factor, because when we dealing with #del>2, we normalize the prob in that group. here we need to rescale 
#'   to the original space (# del is [-4, 15])
#'@param JGenes  gene segment names for each input. this will be used to cacluated overall probability to do normalization.
#'@param deleltionProfile precomputed deletion profile .
#'@return a vector of residuals.
#'@examples
#' #see the example code in the /wl03R2/newPipeline/recomb/recomb_Deletion_v2.0.r
#'@export 

#NOTE: observed need to be scaled/normalized to be a prob over (ndel >2) in the function 
#

residFun_V<-function(pm, observed, xx, JGenes, deletionProfile,  indexes.CDR3 =c(2:7))#,cfs)
{
    if(length(observed)!=dim(xx)[1] || dim(xx)[2]!=6)
    {
        stop("the data input dimensions/lengths are not correct, quit!!!")
    }
    if(length(observed)!=length(JGenes))
    {
        stop("the data input dimensions/lengths are not correct, quit2!!!")
    }
#    if(missing(cfs))
#    {
#        stop("please specify the input correction factors")
#    }
    #first doing
    y<-log(observed)
    pred<-rep(0,length(observed))
    
    for(i in 1:length(observed))
    {
        cdr3<-as.character(t(xx[i,]))
        pred[i]<-contribution_model(pm, cdr3, JGenes[i],deletionProfile,
                            segNames="VGene", ndelNames="RV", indexes.CDR3) #)*cfs[cfs$VGene==JGenes[i], "cf"]
    }
    return (observed-exp(pred)); #exp(pred)
}

#this is the residual function used by 
#pm is the log transformed pwm matrix, is the one to be fitted 
#observed is the probability for each observed CDR3 (6 nt long sub sequence) .
#           vector of numberics
#xx is the subsequence (6nt long) as the input (must be characters)
#           xx is a data.frame with nrow=length(observed) and ncol=6

#'@title residual function for fitting the deletion profile (for J genes)
#'@description 
#'for v genes 
  #'this is the residual function used by by nls.lm fitting the deletion profile.
#'@param pm is the log transformed pwm matrix, is the one to be fitted 
#'@param observed is the probability for each observed CDR3 (6 nt long sub sequence) .
#'           vector of numberics
#'@param xx is the subsequence (6nt long) as the input (must be characters)
#'           xx is a data.frame with nrow=length(observed) and ncol=6
#'----not used for now vvvvvv
#' #@param cfs is the correction factor, because when we dealing with #del>2, we normalize the prob in that group. here we need to rescale 
#'   to the original space (# del is [-4, 15])
#'@param JGenes  gene segment names for each input. this will be used to cacluated overall probability to do normalization.
#'@param deleltionProfile precomputed deletion profile .
#'@return a vector of residuals.
#'@examples
#' #see the example code in the /wl03R2/newPipeline/recomb/recomb_Deletion_v2.0.r
#'@export 

residFun_J<-function(pm, observed, xx, JGenes, deletionProfile,  indexes.CDR3 =c(2:7))
{
    if(length(observed)!=dim(xx)[1] || dim(xx)[2]!=6)
    {
        stop("the data input dimensions/lengths are not correct, quit!!!")
    }
    if(length(observed)!=length(JGenes))
    {
        stop("the data input dimensions/lengths are not correct, quit2!!!")
    }
#    if(missing(cfs))
#    {
#        stop("the data input cfs is mssing, quit........")
#    }
    #first doing
    y<-log(observed)
    pred<-rep(0,length(observed))
    
    for(i in 1:length(observed))
    {
        cdr3<-as.character(t(xx[i,]))
        pred[i]<-contribution_model(pm, cdr3, JGenes[i],deletionProfile, 
                        segNames="JGene", ndelNames="RJ", indexes.CDR3 )
    }
    return (observed-exp(pred))
}

#function to get the 
#'@title cacluated correction factor 
#'@description the correction factor for each sgement is basically the portiton of 
#'  the number of deletion selected (eg n>=2) over all possible deletions.
#'  Note: this one is not used for fitting the parameters, but instead for summarizing and scaling when 
#'          do plotting with all number of delections.
#'@param dat data frame the recombination summary for data analyzed by cloanalyst
#'@param cutoff the low number of deletion to be included in the calculation.
#'@return the portion of n>cutoff of deletion over all possible number of deletions. This portion is used 
#'  to scale the probabilities to plot together with all the deletion numbers.
#'@examples
#'  #see the example code in the /wl03R2/newPipeline/recomb/recomb_Deletion_v2.0.r
#'@export 
getCorrectionFactor<-function(dat, segNames=c("VGene", "JGene"), ndelNames=c("RV", "RJ"),cutoff=2)
{
    segNames<-match.arg(segNames)
    ndelNames<-match.arg(ndelNames)
    #
    cfm<-NULL
    #determine the proportion of n>2
    for(v in unique(dat[,segNames]))
    {
        dat_x<-aggregate(dat[dat[,segNames]==v, ndelNames], 
                            by=list(dat[dat[,segNames]==v, ndelNames]), 
                        FUN=length)
        cf<-sum(dat_x[dat_x[,1]>=cutoff,2])/sum(dat_x[,2])
        cfm<-rbind(cfm, data.frame(cf=cf, VGene=v))
    }
    rownames(cfm)<-cfm$VGene
    colnames(cfm)[2]<-segNames
    return (cfm);
}

#'@title calculate the prob for a CDR3 sequence based on contributions from position weight matrix
#'@description this is different from the icm on that this one is normalized for n>2 deletions.
#'@param segName the name of the gene segment to summary for, "VGene" pr "JGene"
#'@param ndelName the column name for the n of deleltion in the recombination summary  data frame 
#'@param pm is flatten matrix, c(pwm[1,], pwm[3,], pwm[4,]). We don't include 'C' probabilities because 
#'      we want to force a sum of 1 for A/C/T/G. In addition, pm is log'ed matrix.
#'@param CDR3 char array the CDR3 sequence.
#'@param JGenes  gene segment names for each input. this will be used to cacluated overall probability to do normalization.
#'@param deleltionProfile precomputed deletion profile .
#'@param pm is flatten pwm, c(pwm[1,],  pwm[3,], pwm[4,]) # ---excluding pwm[2], pwm[2]=1-rest -----
#'@param indexes.CDR3 numeric array indicating the indexes for CDR3 sequence in the probability distribution of deletion profile
#'@param indexes.CDR3 numeric array indicating the indexes for CDR3 sequence in the probability distribution of deletion profile
#'@return a log'ed total probability for all delections based on the input position weighted matrix.
#'@examples
#' #see the example code in the /wl03R2/newPipeline/recomb/recomb_Deletion_v2.0.r
#'@export
#the return is the log prob with base of euler number
contribution_model<-function(pm=c(pwm[1,],  pwm[3,], pwm[4,]), CDR3, JGene, deletionProfile,
                    segNames=c("VGene", "JGene"), ndelNames=c("RV", "RJ"), indexes.CDR3 =c(2:7)
                    )
{
    segNames<-match.arg(segNames)
    ndelNames<-match.arg(ndelNames)
    
    if(class(CDR3)!="character")
    {
        stop("the input CDR3 is not in the correct format!!! Quit!!!!")
    }
    if(length(CDR3)!=6)
    {
        stop("the input CDR3 doesn't have correct # of elements, quit!!!");
    }
    prob<-icm(pm, CDR3)
    #go through each deletion # in the
    z<-z_sigma(pm, JGene, deletionProfile=deletionProfile, segNames=segNames, ndelNames=ndelNames, indexes.CDR3=indexes.CDR3)
    
    return (prob-z)
}

#'@title sum of the contributions from the deletions profile 
#'@description this is the scaling/normalizing factor for the contribution model.
#'@param segName the name of the gene segment to summary for, "VGene" pr "JGene"
#'@param ndelName the column name for the n of deleltion in the recombination summary  data frame 
#'@param pm is flatten matrix, c(pwm[1,], pwm[3,], pwm[4,]). We don't include 'C' probabilities because 
#'      we want to force a sum of 1 for A/C/T/G. In addition, pm is log'ed matrix.
#'@param CDR3 char array the CDR3 sequence.
#'@param JGenes  gene segment names for each input. this will be used to cacluated overall probability to do normalization.
#'@param deleltionProfile precomputed deletion profile .
#'@param pm is flatten pwm, c(pwm[1,],  pwm[3,], pwm[4,]) # ---excluding pwm[2], pwm[2]=1-rest -----
#input pm is flatten pwm, c(pwm[1,],  pwm[3,], pwm[4,]), ---exclusing pwm[2,]----
#'@return a log'ed probability for the CDR3 based on the input position weighted matrix.

#J Gene is the J 
z_sigma<-function(pm=c(pwm[1,], pwm[3,], pwm[4,]), JGene, deletionProfile, 
                segNames=c("VGene", "JGene"), ndelNames=c("RV", "RJ"), indexes.CDR3=c(2:7))
{
    segNames<-match.arg(segNames)
    ndelNames<-match.arg(ndelNames)
    j_profile<-deletionProfile[deletionProfile[,segNames]==JGene,]
    z<-0
    for(n in j_profile[,ndelNames])
    {
        cdr3<-j_profile[j_profile[,ndelNames]==n, indexes.CDR3]
        cdr3<-as.character(t(cdr3))
        z<-z+exp(icm(pm, cdr3))
    }
    return(log(z))
}
#'@title individual contribution model 
#'@description individual contribution model based on pwm, this is unnormalized version.
#'just add up the "contribution" based on the 6nts surranding the deletion position
#'assuming pwm is log transformed.
#'the returning value is not prob yet, is the "contribution" without normalization and log'ed.
#'@param pm is flatten matrix, c(pwm[1,], pwm[3,], pwm[4,]). We don't include 'C' probabilities because 
#'      we want to force a sum of 1 for A/C/T/G. In addition, pm is log'ed matrix.
#'@param CDR3 char array the CDR3 sequence.
#'@return log'ed probality for the observed CDR3 (6 nts around the deleted nts).
#'@examples
#' #see the example code in the /wl03R2/newPipeline/recomb/recomb_Deletion_v2.0.r
#'
#'@export
icm<-function(pm=c(pwm[1,], pwm[3,], pwm[4,]), CDR3)
{
    if(class(CDR3)!="character")
    {
        stop("the input CDR3 is not in the correct format!!! Quit!!!!")
    }
    if(length(CDR3)!=6)
    {
        stop("the input CDR3 doesn't have correct # of elements, quit!!!");
    }
    CDR3<-factor(CDR3,levels=c("A", "C", "G","T", "-"))
    CDR3<-as.integer(CDR3)
    prob=0
    #i is the position of CDR3 1-6
    pm_full<-rep(0, 5*6)
    pm_full[1:6]<-pm[1:6]
    pm_full[(1+6*2):(6+6*2)]<-pm[(7:12)]
    pm_full[(1+6*3):(6+6*3)]<-pm[13:18]
    pm_full[7:12]<-log(rep(1,6)-exp(pm[1:6])-exp(pm[7:12])-exp(pm[13:18]))
    pm_full[25:30]<-rep(0,6)
    for(i in 1:6){
        prob=pm_full[(CDR3[i]-1)*6+i]+prob
    }
    return (prob)
}

#'@title extract deletion profile prob distribution 
#'@description calculate the deletion profile proble distribution for each gene segment based on the 
#'  recombination summary data frame.
#'@param rdat data frame the recombination summary data for each tissue+isotype 
#'@param low the lowest # of deleted nt for the distribution 
#'@param high the highest # of deleted nt for the distribution
#'@param segName the name of the gene segment to summary for, "VGene" pr "JGene"
#'@param ndelName the column name for the n of deleltion in the recombination summary  data frame 
#'@examples
#' #see the example code in the /wl03R2/newPipeline/recomb/recomb_Deletion_v2.0.r
#'@export
probDistDel<-function(rdat, low=2, high=18, segName=c("VGene", "JGene"), ndelName=c("RV", "RJ"  ))
{
        segName<-match.arg(segName)
        ndelName<-match.arg(ndelName)
        
        gs.index<-grep(segName, colnames(rdat))
        if(length(gs.index)!=1)
        {
            stop("the specified geneSeg name is not correct, please check. Quit.....")
        }
        ndel.index<-grep(ndelName, colnames(rdat))
        if(length(ndel.index)!=1)
        {
            stop("the specified ndel name is not correct, please check. Quit.....")
        }
        if(missing(rdat))
        {
            stop("please specify rdat (recomb summary data frame), quit.......")
        }
        rdat<-rdat[rdat[,ndel.index]<=high&rdat[,ndel.index]>=low,]
        rdat[,gs.index]<-as.character(rdat[,gs.index])
        IgVs<-unique(rdat[,gs.index])
        p_vs<-NULL
        for(v in IgVs)
        {
            #Get p_j(), dependent variable as inputs
             x<-aggregate(rdat[rdat[,gs.index]==v, ndel.index], 
                            by=list(rdat[rdat[,gs.index]==v, ndel.index]), 
                        FUN=length)
            x[,2]<-x[,2]/sum(x[,2])
            x$VGene<-v
            colnames(x)<-c(ndelName, "Prob", segName)
            
            p_vs<-rbind(p_vs, x)
        }
        return (p_vs)
}


#generate deletion sequence based on data as inputs)

#'@title based on the input to create table of deletion sequence profile for each gene segment 
#'@description NOTE: this is for V genes. We will have a different function J genes, since the way to 
#'      to get CDR3 is different in direction 3'->5' for V genes and 5'->3' for J genes . 
#'  NOTE: this is now hardcoded for VGene, will need to revise later.
#  #' @param rdat data frame of the recom data 
#'@param probDist data frame the prob of each ndel in each sample for each gene segment 
#'@param anns the input data table cross referencing the gene segment id and sequences. 
#'@param segSeqs data frame providing the gene segment sequences.
#'@param segName the name of the gene segment to summary for, "VGene" pr "JGene"
#'@param ndelName the column name for the n of deleltion in the recombination summary  data frame 
#'@param ndel the lowest number of deletions. We will calculate all deletions larger than this ndel, for all n>ndel   
#'          NO!!! the above description for ndel is NOT right. ndel is now specifying the number of  nts before the 
#'              deletion site. The number of deletion is specified by the porbDist$RJ.  They are talking about two different 
#'              things.
#'@examples 
#' #see the example code in the /wl03R2/newPipeline/recomb/recomb_Deletion_v2.0.r
#'@export
deletionProfileData<-function(probDist, anns, segSeqs, segNames=c("VGene", "JGene"), 
                                                            ndelNames=c("RV", "RJ"), ndel=2)
{
    if(missing(probDist)||missing(anns)||missing(segSeqs))#||missing(rdat))
    {
        stop("please specify the input data.....")
    }
    
    segNames<-match.arg(segNames)
    #cat("1\n")
    ndelNames<-match.arg(ndelNames)
    #cat("2\n")
    if(ndel<0)
    {
        stop("the input ndel is smaller than zero. we haven't implement such condition. quit.")
    }
    v<-NULL
    IgVs<-unique(probDist[,segNames])
    x<-NULL
    #cat("start doing ...............loop\n")
    for(j in IgVs)#for each J gene sequences, we need to figure out the pwm. This is similiar to multiple  factor regression
                        # log(p)= x1+x2+x3+x4+x5+x6.    p is prob( n del | CDR3 & n>2) and x1...x6 are the nts at six positions around the deletion site. 
                        # NOTE: the 6 position is the 2 positions in the deleted portion and 4 positions in the undeleted portion.
                        #for each position p, we will figure out the input ( the nt based on the j gene sequences), A/C/T/G
    {
        cat("....")
        #get the sequence and get deletion profile (# of del)
        nums<-probDist[probDist[,segNames]==j, ndelNames]
        #cat("\t numbs :", nums,"\n")
        uid<-anns[anns$Name==j,"UID"][1]
        jseq<-(segSeqs[uid][[1]])
        ln_seq<-length(jseq)
        x_j<-NULL
        #cat("doing inner loop...\n")
        for(n in nums)
        {
            if(segNames=="JGene")
            {
                xstart<-(n-ndel+1)
                #xend<-(n-ndel+5)
                x_temp<-data.frame(RJ=n, x1=jseq[xstart],x2=jseq[xstart+1],x3=jseq[xstart+2], x4=jseq[xstart+3], x5=jseq[xstart+4], x6=jseq[xstart+5] )
                names(x_temp)[1]<-ndelNames
                x_j<-rbind(x_j, x_temp)
            }
            else 
            {
                xstart<-ln_seq-(n-ndel)
                xend<-ln_seq-(n-ndel+5)
                x_temp<-data.frame(RJ=n, x1=jseq[xstart],x2=jseq[xstart-1],x3=jseq[xstart-2], x4=jseq[xstart-3], x5=jseq[xstart-4], x6=jseq[xstart-5] )
                names(x_temp)[1]<-ndelNames
                x_j<-rbind(x_j, x_temp)
            }
        }
        #cat("Done inner loop...\n")
        x_j[,segNames]<-j
        x<-rbind(x,x_j)
    }
    p_x<-cbind(data.frame(Prob=probDist$Prob), x)
    p_x$log_p<-log(p_x$Prob)
    cat("Done\n")
    return(p_x)
}

#'@title build the deletion profile for v gene segments
#'@description this is the deletion profile (all number of deleted ) for v gene segments of the observed data.
#'
#'@param rdata data frame for observed data results by cloanalyst. 
#'@param anns annotations for all gene segments 
#'@param segSeqs gene sequences for the segments 
#'@param nlow lowest number of deleted nts.
#'@param nhigh highest number of deleted nts.
#'@param segNames name of gene segment. Can only be one of c("VGene", "JGene") 
#'@param ndelNames name of deletion. Can only be one of c("RV", "RJ")
#'@param ndel the number of nts before the delections to be used as the CDR3 sequences for checking 
#'          this is different from nlow, which is the number of lowest deletion to check as the deletion profile .
#'          nlow must be larger than ndel. 
#'@examples
#'   #see the example code in the /wl03R2/newPipeline/recomb/recomb_Deletion_v2.0.r
#'@export
deletionProfileAll<-function(rdat, anns, segSeqs, nlow=2, nhigh=15,segNames=c("VGene", "JGene"), 
                                                            ndelNames=c("RV", "RJ"), ndel=nlow)
{
    if(missing(rdat)||missing(anns)||missing(segSeqs))#||missing(rdat))
    {
        stop("please specify the input data.....")
    }
    
    segNames<-match.arg(segNames)
    #cat("1\n")
    ndelNames<-match.arg(ndelNames)
    #cat("2\n")
    if(nlow<0)
    {
        stop("the input ndel is smaller than zero. we haven't implement such condition. quit.")
    }
    if(nlow < ndel )
    {
        stop("the input ndel is smaller than ndel. we haven't implement such condition. quit.")
    }
    
    IgVs<-unique(rdat[,segNames])
    
    deletionProfile<-NULL
    for(j in IgVs)#for each J gene sequences, we need to figure out the pwm. This is similiar to multiple  factor regression
                        # log(p)= x1+x2+x3+x4+x5+x6.    p is prob( n del | CDR3 & n>2) and x1...x6 are the nts at six positions around the deletion site. 
                        # NOTE: the 6 position is the 2 positions in the deleted portion and 4 positions in the undeleted portion.
                        #for each position p, we will figure out the input ( the nt based on the j gene sequences), A/C/T/G
    {
        cat(".....")
        #get the sequence and get deletion profile (# of del)
        nums<-c(nlow:nhigh)
        uid<-anns[anns$Name==j,"UID"][1]
        jseq<-(segSeqs[uid][[1]])
        ln_seq<-length(jseq)
        x_j<-NULL
        for(n in nums)
        {
            if(segNames=="JGene")
            {
                xstart<-(n-ndel +1)
                #xend<-(n-ndel+5)
                x_temp<-data.frame(RJ=n, x1=jseq[xstart],x2=jseq[xstart+1],x3=jseq[xstart+2], x4=jseq[xstart+3], x5=jseq[xstart+4], x6=jseq[xstart+5] )
                names(x_temp)[1]<-ndelNames
                x_j<-rbind(x_j, x_temp)
            }
            else 
            {
                xstart<-ln_seq-(n-ndel)
                xend<-ln_seq-(n-ndel+5)
                x_temp<-data.frame(RJ=n, x1=jseq[xstart],x2=jseq[xstart-1],x3=jseq[xstart-2], x4=jseq[xstart-3], x5=jseq[xstart-4], x6=jseq[xstart-5] )
                names(x_temp)[1]<-ndelNames
                x_j<-rbind(x_j, x_temp)
            }
            #xstart<-lseq-(n-nlow)
            #xend<-lseq-(n-nlow+5)
            #x_temp<-data.frame(RJ=n, x1=jseq[xstart],x2=jseq[xstart-1],x3=jseq[xstart-2], x4=jseq[xstart-3], x5=jseq[xstart-4], x6=jseq[xstart-5] )
            #x_v<-rbind(x_v, x_temp)
        }
        names(x_j)[1]=ndelNames
        x_j[,segNames]<-j
        deletionProfile<-rbind(deletionProfile,x_j)
    }
    cat("\n")
    return (deletionProfile);
}
#now, we need to build the model "manual" and do the fitting.
#'@title  build the position weight matrix based on observed data
#'@description this is build on data and will be used as the starting points for fitting.
#'@param probDist this is the matrix generated by the function deletionProfileV with data.
#'@param nrow number of rows in the position weighted matrix. 4 by default
#'@param ncol number of columns in the position weighted matrix. 6 by default 
#'@export 
#'@examples 
#'  #see the example code in the /wl03R2/newPipeline/recomb/recomb_Deletion_v2.0.r
buildPWM<-function(probDist, nrow=4, ncol=6)
{ 
    pwm<-matrix(0,nrow=4, ncol=6)
    #get the start values for pwm
    for(i in 3:8)
    {
        tmp<-aggregate(probDist[,i], by=list(probDist[,i]), length)
        if(dim(tmp)[1]>4){
            tmp<-tmp[1:4,]
            }
        pwm[,i-2]<-tmp[,2]/sum(tmp[,2])
        
    }
    rownames(pwm)<-tmp[,1]
    pwm<-log(pwm)
    return (pwm)
}