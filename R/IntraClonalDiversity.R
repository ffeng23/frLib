#This is the module contains functions to do intraClonalDiversity of IgSeq data.
#     The R code here are more for testing and concept-proving purpose. We have the c code to the real job, since R code run too slow.
#  Note: here the intraClonalDiversities function is not "exported". Can not be called directly.
#               11/16/2020
#
#'@title r function to determine the sequence length with sequences containing discountinued section.
#'@description in this case the two discontinued sections are well included in
#' the sequences (see below)
#' seq1 |----<.........>-----|
#' seq2    |----<.........>---|
#' but the boundaries of the two disc section could be falling into different case, crossed, separated
#'       or one is ahead of the other, etc.
#'Also we 1)ASSUME seqLength is the SHORTER length between seq1 and seq2 and we don't have to know whether it is seq1 or seq2.
#'2) ASSUME the seq1 disc section (discPositionStart1) is before seq2 disc section.
#'3) both discontinued sections are present 
#'@param seqLength int length of the shorter sequence between the two. It doesn't matter which one. We don't 
#       have to know which one this is.
#'@param discPositionStart1 start position of seq1. assume to be the discountinued section located ahead of the other one (seq2)
#'@param discPositionEnd1   end position of seq1.
#'@param discPositionStart2 start position of seq2, assume to be the after the start of seq1 (above)
#'@param discPositionEnd2 end postion of seq2
#'@return sequence length which the discountinued section has been subtracted from.
#'@export
#'@examples 
#'       #first test the first function 
#'  seqLength<-290
#'  discPositionStart1<-11
#'        discPositionEnd1<-20
#'        discPositionStart2<-28
#'       discPositionEnd2<-29
#'
#' getSeqLength_inclusive_seq1Disc_ahead(
#'        seqLength,
#'        discPositionStart1,
#'        discPositionEnd1,
#'        discPositionStart2,
#'        discPositionEnd2
#'    )
#'
#'  seqLength<-290
#'  discPositionStart1<-11
#'        discPositionEnd1<-20
#'        discPositionStart2<-21
#'       discPositionEnd2<-29
#'
#' getSeqLength_inclusive_seq1Disc_ahead(
#'        seqLength,
#'        discPositionStart1,
#'        discPositionEnd1,
#'        discPositionStart2,
#'        discPositionEnd2
#'    )
#'
getSeqLength_inclusive_seq1Disc_ahead <-function (
        seqLength,
        discPositionStart1,
        discPositionEnd1,
        discPositionStart2,
        discPositionEnd2
    )
{
        if(is.na(discPositionStart1)||is.na(discPositionStart2)||is.na(discPositionEnd1)||is.na(discPositionEnd2))
        {
            stop("one or both sequence discontinued section is not present. quit!!!")
        }
        if(discPositionStart1>discPositionStart2)
        {
            stop("the seq 1 discontinued section starts after seq2 one. quit!!")
        }
#NOTE: again assume discPositionStart1<discPositionStart2
        mg_discPositionStart<-discPositionStart1; #real start ; set default 
        mg_discPositionEnd<-discPositionEnd1;#read end ; set default
        num_nt<-seqLength;#mg_discPositionStart-mg_discPositionEnd;# the real final length; set default 
        if(discPositionEnd1<discPositionEnd2)
        {
            if(discPositionEnd1<discPositionStart2){
                    #two disc sections don't cross, we need to subtract both
                    num_nt<-seqLength-(mg_discPositionEnd-mg_discPositionStart+1);#default seq1 first
                    #seq2
                     mg_discPositionStart<-discPositionStart2;
                    mg_discPositionEnd<-discPositionEnd2;
                    num_nt<-num_nt-(mg_discPositionEnd-mg_discPositionStart+1);
            }else #both cross , since seq1 disc ends after seq2 starts.
            {
                mg_discPositionStart<-discPositionStart1;
                    mg_discPositionEnd<-discPositionEnd2;  
                    num_nt<-seqLength-(mg_discPositionEnd-mg_discPositionStart+1);
            }
        }
        else  #in this case, seq1 end after seq2 disc end, seq1 include seq2 
        {
            #using the default setting (of seq1 disc )
            num_nt<-seqLength-(mg_discPositionEnd-mg_discPositionStart+1);
        }
        #mg_discPositionStart<-discPositionStart.seq2
        return (num_nt);
}
#'@title caculate the actual length with seq1 discontinued section
#'@description 
#'this function deal with the condition where the seq1 discontinued section exists, but the seq2 disc section might or might not exists
#'most times seq1 total length and seq2 total length are the same.
#' we also assume seq1 is the shorter sequence between the two. 
#' do NOT assume discPosition 1 is ahead of discPosition2
#'@export
#'@examples
#'#now testing second function.
#‘
#' seqLength<-290
#' totalLength<-300
#' discPositionStart1<-11
#'        discPositionEnd1<-20
#'        discPositionStart2<-NA
#'        discPositionEnd2<-NA
#'  getSeqLength_disc1_existShort (
#'        seqLength,
#'        totalLength,
#'        discPositionStart1,
#'        discPositionEnd1,
#'        discPositionStart2,
#'        discPositionEnd2
#'    )
#'
#' seqLength<-290
#' totalLength<-300
#'  discPositionStart1<-11
#'        discPositionEnd1<-20
#'        discPositionStart2<-30
#'        discPositionEnd2<-100
#' getSeqLength_disc1_existShort (
#'        seqLength,
#'        totalLength,
#'        discPositionStart1,
#'        discPositionEnd1,
#'        discPositionStart2,
#'        discPositionEnd2
#'    )
getSeqLength_disc1_existShort<-function (
        seqLength,
        totalLength,
        discPositionStart1,
        discPositionEnd1,
        discPositionStart2,
        discPositionEnd2
    )
    {
        if(is.na(discPositionStart1)||is.na(discPositionStart1))
        {
            stop("seq1 discontinued section is not present!!! quit!!!")
        }
        #ASSUME disc seq1 exist and it is short than seq2
        #so we don't have to check the whether seq1 disc exist
        #set up the default 
        
        mg_discPositionStart<-discPositionStart1;
        mg_discPositionEnd<-discPositionEnd1;
        num_nt<-seqLength - (mg_discPositionEnd-mg_discPositionStart+1) 
                        
        #check the longer one disc positions, in case the short sequence start within the 
        #longer sequence discontinue section
        if(  !is.na(discPositionStart2) )  #case where seq2 disc exists 
        {
            #first check the location of seq1 start and seq2 disc section starts
             #the first case is seq1 starts after seq2 disc section starts
            if( (totalLength - seqLength) > discPositionStart2){  # note the seq1 and seq2 are all right aligned and they must have identical CDR3
                        if((totalLength - seqLength)<discPositionEnd2)
                        { #staring of seq1 is falling into disc section of seq2. so we need to reset seq2 disc section starting position 
                                discPositionStart2<- totalLength - seqLength;
                                #after resetting the starting point, we can not call the function do the job
                                #prepare the input 
                                
                                num_nt<-getSeqLength_inclusive_seq1Disc_ahead (#note: in this case seq2 disc is ahead of seq1 disc, so we swap the position and feed int the params 
                                            seqLength, discPositionStart2, discPositionEnd2, 
                                            discPositionStart1, discPositionEnd1
                                        );
                        }else {
                            # starting of seq1 is not include disc section , this is same as if there is no disc section on seq2 and we only need to look at seq1 
                            # so we do noghting in here. since originally above it is check seq1 alone.
                            #use the default setting the seq1 disc section boundaries
                            num_nt<-seqLength-(mg_discPositionEnd-mg_discPositionStart+1);
                        }
                
            } else #meaning (totalLength - num.seq1)<=discPositionStart.seq2)  in this case we have disc seq1 and disc seq2 exist and they all well are included in the shorter seq 
            {
                #now we need to check to see which one is ahead of the other to determine the order of 
                #calling 
                if(discPositionStart1<discPositionStart2){
                    num_nt<-getSeqLength_inclusive_seq1Disc_ahead ( #<=== function call need to filling detail
                            seqLength, discPositionStart1, discPositionEnd1,
                                            discPositionStart2, discPositionEnd2
                        )
                }else { # this is the one disc seq2 is ahead
                     num_nt<-getSeqLength_inclusive_seq1Disc_ahead (#note: in this case seq2 disc is ahead of seq1 disc, so we swap the position and feed int the params 
                                            seqLength, discPositionStart2, discPositionEnd2, 
                                            discPositionStart1, discPositionEnd1
                                        );                    
                }
                    
            }
        }   else   { #in this condition, no disc section in seq2
            # we do nothing, since the originally above it is set up to consider only the seq1
                    #jusing the default setting seq1 boundaries
                    #use the default setting the seq1 disc section boundaries
                    num_nt<-seqLength-(mg_discPositionEnd-mg_discPositionStart+1);
             
        }#end of if else first level
        
        return (num_nt);
} #end of function 
#'@title to calculate the actual sequence length with seq1 discontinued section not present
#'@description in this function, 
#' we process the case where seq 1 is shorter, 
#' and no discountinue section. 
#'@export
#'@examples 
#' #test the last one 
#'  seqLength<-290
#'totalLength<-300
#' 
#'        discPositionStart2<-NA
#'        discPositionEnd2<-NA 
#'        getSeqLength_disc1_NAShort(
#'        seqLength,
#'        totalLength,
#'        discPositionStart2,
#'        discPositionEnd2
#'    )
#'  seqLength<-290
#'  totalLength<-300
#' 
#'        discPositionStart2<-30
#'        discPositionEnd2<-39
#'        getSeqLength_disc1_NAShort(
#'        seqLength,
#'        totalLength,
#'        discPositionStart2,
#'        discPositionEnd2
#'    )
getSeqLength_disc1_NAShort<-function (
        seqLength,
        totalLength,
        discPositionStart2,
        discPositionEnd2
    )
{
    num_nt<-seqLength
    if(!is.na(discPositionStart2)){
        if((totalLength-seqLength)<discPositionStart2){
            #here the shorter seq1 is before seq2 discountinue section starts
                num_nt<-seqLength-(discPositionEnd2-discPositionStart2+1);
        } else { # here the shoter seq1 is after seq2 discontinue section starts 
            if((totalLength-seqLength)<discPositionEnd2){ # shorter seq1 starts in the middle of discontinued section 
                    discPositionStart2<-totalLength-seqLength
                    num_nt<-seqLength-(discPositionEnd2-discPositionStart2+1);
            } else { # the shorter seq1 is after seq2 discontinued section ends 
                #so there is no discontinued section 
                num_nt<-seqLength;
            }
        }
    }
    else {#both seq2 discontinued and seq1 discontinued section are not present .
        num_nt<-seqLength; #do nothing.
    }
    return(num_nt); 
}
#'@title calculate intraclonal diversity for one clone
#'
#'@description    
#'now start doing parsing of the pairwise difference
#'return a list showing idi, mutations.ns (non-silent), mutations (total)
#'assumming in this function, we only do for one specific clone. 
#'it is the outside caller's responsibility to clean and feed in the 
#'correct input only relevant to this specific clone.
#clone.assignment=clone.S10.CID32;			mutations<-mutations.S10.CID32; ig.recSum<-ig.S10.CID32;indel.penalty=1;			vlengths<-vlength.S10.CID32;
#'@examples
#' #note: the following code example needs read the data from saved R data. 
#' # the original R data were saved in the wl03R2 project. 
#'  #make sure they can be load. see /home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline/DataAnalysis_v1.2_IntraClonalDiversity_nonSilentMu_cfun.R
#'  #setwd("/home/feng/Feng/LAB/Wetzler_PorinB/wl03R2/newPipeline")
#'
#'#try to get bm sample 10 cloneID 32
#'df.cloneAssign<-BM.df.cloneAssign[["IgG"]]
#'IG.RecSum<-IG.RecSum.bm[["IgG"]]
#'df.mutations<-BM.df.mutations[["IgG"]]
#'df.vlength<-BM.df.vlength[["IgG"]]
#'
#'clone.S10.CID32<-df.cloneAssign[df.cloneAssign$sampleName=="MB10"&df.cloneAssign$CloneID==1,]
#'ig.S10.CID32<-IG.RecSum[is.element(IG.RecSum$UID, clone.S10.CID32$ReadID),]
#'mutations.S10.CID32<-df.mutations[is.element(df.mutations$ReadID, clone.S10.CID32$ReadID),]
#'vlength.S10.CID32<-df.vlength[is.element(df.vlength$ReadID, clone.S10.CID32$ReadID),]
#'
#'clone.S10.CID12<-df.cloneAssign[df.cloneAssign$sampleName=="MB10"&df.cloneAssign$CloneID==12,]
#'ig.S10.CID12<-IG.RecSum[is.element(IG.RecSum$UID, clone.S10.CID12$ReadID),]
#'mutations.S10.CID12<-df.mutations[is.element(df.mutations$ReadID, clone.S10.CID12$ReadID),]
#'vlength.S10.CID12<-df.vlength[is.element(df.vlength$ReadID, clone.S10.CID12$ReadID),]
#'
#'clone.S1.CID26<-df.cloneAssign[df.cloneAssign$sampleName=="MB1"&df.cloneAssign$CloneID==26,]
#'ig.S1.CID26<-IG.RecSum[is.element(IG.RecSum$UID, clone.S1.CID26$ReadID),]
#'mutations.S1.CID26<-df.mutations[is.element(df.mutations$ReadID, clone.S1.CID26$ReadID),]
#'vlength.S1.CID26<-df.vlength[is.element(df.vlength$ReadID, clone.S1.CID26$ReadID),]
#'
#'  #### idi=0, no diversity.
#'	idi<-intraClonalDiversity(clone.assignment=clone.S1.CID26,
#'		mutations<-mutations.S1.CID26,
#'			ig.recSum<-ig.S1.CID26,
#'			indel.penalty=1,
#'			vlengths<-vlength.S1.CID26
#'			)
#'  ##### make it diversified
#'      mutations.S1.CID26[15,"Position"]<-259
#'	    idi<-intraClonalDiversity(clone.assignment=clone.S1.CID26,
#'		    mutations<-mutations.S1.CID26,
#'			ig.recSum<-ig.S1.CID26,
#'			indel.penalty=1,
#'			vlengths<-vlength.S1.CID26
#'			) #idi=0.005475725
#'  ###make one fallen into discontinued section
#'    mutations.S1.CID26[15,"Position"]<-150
#'	    idi<-intraClonalDiversity(clone.assignment=clone.S1.CID26,
#'		    mutations<-mutations.S1.CID26,
#'			ig.recSum<-ig.S1.CID26,
#'			indel.penalty=1,
#'			vlengths<-vlength.S1.CID26
#'			) #0.002737862
#'    ###modify the discontinued section 
#'  vlength.S1.CID26[3,"discPosStart"]<-161
#'idi<-intraClonalDiversity(clone.assignment=clone.S1.CID26,
#'		    mutations<-mutations.S1.CID26,
#'			ig.recSum<-ig.S1.CID26,
#'			indel.penalty=1,
#'			vlengths<-vlength.S1.CID26
#'			) #0.002737862
#'  vlength.S1.CID26[3,"discPosStart"]<-163
#' vlength.S1.CID26[3,"discPosEnd"]<-164
#' idi<-intraClonalDiversity(clone.assignment=clone.S1.CID26,
#'		    mutations<-mutations.S1.CID26,
#'			ig.recSum<-ig.S1.CID26,
#'			indel.penalty=1,
#'			vlengths<-vlength.S1.CID26
#'			) #0.002760536
#'
#' #set up to have NA discontinued section.
#'vlength.S1.CID26[1,c("discPosStart", "discPosEnd")]<-c(NA, NA)
##' idi<-intraClonalDiversity(clone.assignment=clone.S1.CID26,
#'		    mutations<-mutations.S1.CID26,
#'			ig.recSum<-ig.S1.CID26,
#'			indel.penalty=1,
#'			vlengths<-vlength.S1.CID26
#'			)  # idi=0.003714128
#' #set up to have NA on discontinued section 
#' #vlength.S1.CID26<-df.vlength[is.element(df.vlength$ReadID, clone.S1.CID26$ReadID),]
#' vlength.S1.CID26[2,c("discPosStart", "discPosEnd")]<-c(NA, NA)
#'  #counting the mutations 1 v 2, 0 mutations; 1 v 3, 1 mutations (288-(161-118+1+164-163+1))
#'  # 2 vs 3, 2 mutations (288-(164-163+1))
##' idi<-intraClonalDiversity(clone.assignment=clone.S1.CID26,
#'		    mutations<-mutations.S1.CID26,
#'			ig.recSum<-ig.S1.CID26,
#'			indel.penalty=1,
#'			vlengths<-vlength.S1.CID26
#'			)  # idi=0.003708413
#' # 2 have no disContinued section.  
#'  vlength.S1.CID26[1,c("discPosStart", "discPosEnd")]<-c(NA, NA)
#'idi<-intraClonalDiversity(clone.assignment=clone.S1.CID26,
#'                 mutations<-mutations.S1.CID26,
#'                     ig.recSum<-ig.S1.CID26,
#'                     indel.penalty=1,
#'                     vlengths<-vlength.S1.CID26
#'                     ) #idi=0.004662005 
#' #all 3 have no disc section.
#' vlength.S1.CID26[3,c("discPosStart", "discPosEnd")]<-c(NA, NA)
#'idi<-intraClonalDiversity(clone.assignment=clone.S1.CID26,
#'                 mutations<-mutations.S1.CID26,
#'                     ig.recSum<-ig.S1.CID26,
#'                     indel.penalty=1,
#'                     vlengths<-vlength.S1.CID26
#'                     ) #idi=0.00462963
#
#'@export
intraClonalDiversity<-function(  
					clone.assignment, #this is the sequences for this specific clone
					mutations,  #this is the mutations for the sequences for this clone 
					ig.recSum, #contains the VRG info for the sequences in this clone.
					vlengths, #contains the vlength for VRG matches and total v length
					indel.penalty=1  #how many substituations worth of one indel???, it could be zero
				)
	{
	   #input:
	   #clone.assignment<-clone.bm.S10.CID32;mutations<-mutations.bm.S10.CID32;ig.recSum<-ig.bm.S10.CID32;indel.penalty=1;vlengths<-vlength.bm.S10.CID32
	   #clone.assignment<-clone.bm.S10.CID12;mutations<-mutations.bm.S10.CID12;ig.recSum<-ig.bm.S10.CID12;indel.penalty=1;vlengths<-vlength.bm.S10.CID12
		#for each record/sequence,
		#we need to check the mutations table/data frame to count the pairwise difference
		num.pairs<-0;
		pair.diff<-0;
		ig.recSum$UID<-as.character(ig.recSum$UID)
		mutations$ReadID<-as.character(mutations$ReadID)
		clone.assignment$ReadID<-as.character(clone.assignment$ReadID)
		#here we assume that we have more than 1 sequences in the clones.
		#we don't consider that singleton clones.
		if(dim(clone.assignment)[1]<=1)
		{
			warning("singleton input to the function (intraClonalDiversity)!!!\n")
			return (0);
		}
		for(i in 1:dim(clone.assignment)[1])
		{
			#get the sequences
			id.seq1<-as.character(clone.assignment[i,"ReadID"]);
			
			num_nt.seq1<-ig.recSum[ig.recSum$UID==id.seq1,"X.VBases"];
			#cat("i:", i,"\n")
			if(i>=dim(clone.assignment)[1])
			{
				break; #we are done.
			}
			
			for(j in (i+1):dim(clone.assignment)[1])
			{
				#cat("\tj:", j,"\n")
				num.pairs<-num.pairs+1;
				mutation.seq1<-mutations[mutations$ReadID==id.seq1,c("Type", "Position")]
				totalV.seq1<-vlengths[vlengths$ReadID==id.seq1,"totalVBase"];
				#read the mutations from the mutations file
				#the logic is that the mutations is comparing with the UCA, 
				#list the point mutations at each
				#first total number of nts in the sequences
                
                #------updated 11/12/2020
                    #now we need to add one extra condition to get rid of mutations that is falling into the discontinued region of the other seq
                #for seq 1
                discPositionStart.seq1<-vlengths[vlengths$ReadID==id.seq1,"discPosStart"];#<----------
                discPositionEnd.seq1<-vlengths[vlengths$ReadID==id.seq1,"discPosEnd"];#<-----------
                
				id.seq2<-clone.assignment[j,"ReadID"];
				num_nt.seq2<-ig.recSum[ig.recSum$UID==id.seq2, "X.VBases"];
				totalV.seq2<-vlengths[vlengths$ReadID==id.seq2,"totalVBase"];
				
				mutation.seq2<-mutations[mutations$ReadID==id.seq2,c("Type", "Position")]
                discPositionStart.seq2<-vlengths[vlengths$ReadID==id.seq2,"discPosStart"];#<-------
                discPositionEnd.seq2<-vlengths[vlengths$ReadID==id.seq2,"discPosEnd"];#<--------
                
                #now determine the valid mutations<--------
                #seq1
                if((!is.na(discPositionStart.seq2))&&(!is.na(discPositionEnd.seq2)))
                {
                    mutation.seq1<-mutation.seq1[mutation.seq1$Position<discPositionStart.seq2|mutation.seq1$Position>discPositionEnd.seq2,]
                }
                #seq2
                if((!is.na(discPositionStart.seq1))&&(!is.na(discPositionEnd.seq1)))
                {
                    mutation.seq2<-mutation.seq2[mutation.seq2$Position<discPositionStart.seq1|mutation.seq2$Position>discPositionEnd.seq1,]
                }
				#==============start determining the valid sequence length take care of the discontinued sections 
                #id.long<-id.seq1;
				num_nt<-num_nt.seq2;
				if(num_nt.seq2>num_nt.seq1){ #seq1 is shorter.
                    #in here we assume the sequence is right aligned.!!!(because they have identical CDR3 region)
					num_nt<-num_nt.seq1; #seq 1 is shorter, we are on seq1 and checking against seq2
                    if(!is.na(discPositionStart.seq1)&&!is.na(discPositionEnd.seq1)) #<----
                    {
                        num_nt<-getSeqLength_disc1_existShort (
                                num_nt.seq1,
                                totalV.seq1,
                                discPositionStart.seq1,
                                discPositionEnd.seq1,
                                discPositionStart.seq2,
                                discPositionEnd.seq2
                            );
                    }else{ #seq1 doesn't has disc section 
                          num_nt<-getSeqLength_disc1_NAShort(
                                    num_nt.seq1,
                                    totalV.seq1,
                                    discPositionStart.seq2,
                                    discPositionEnd.seq2
                                );
                    }
                    
					mutation.seq2<-mutation.seq2[mutation.seq2$Position>(totalV.seq1-num_nt.seq1),];
					# the reason that we do this is because the totalV is the expect length, but sequencing might give us a shortened/truncated one.
					#  num_nt.seq1 or num_nt.seq2 are the observed length, could be shorter.  (totalV.Seq- num_nt.seq1) telling us the unsequenced nts.
					# the mutation has to be falling into the sequenced region. 
                    
                }else {  #in this case num_nt.seq2< num_nt.seq1, seq2 is shorter
					num_nt<-num_nt.seq2;
                    if(!is.na(discPositionStart.seq2)&&!is.na(discPositionEnd.seq2)) #<----
                    {
                        num_nt<-getSeqLength_disc1_existShort (
                                num_nt.seq2,
                                totalV.seq2,
                                discPositionStart.seq2,
                                discPositionEnd.seq2,
                                discPositionStart.seq1,
                                discPositionEnd.seq1
                            );
                        
                    } else { # seq2 disc is not present
                        num_nt<-getSeqLength_disc1_NAShort(
                                    num_nt.seq2,
                                    totalV.seq2,
                                    discPositionStart.seq1,
                                    discPositionEnd.seq1
                                );
                    }
					mutation.seq1<-mutation.seq1[mutation.seq1$Position>(totalV.seq2-num_nt.seq2),]
                }
             
				#now count the different
				num.diff.sub<-length(setdiff(mutation.seq1[mutation.seq1$Type=="Substitution","Position"],mutation.seq2[mutation.seq2$Type=="Substitution","Position"]))
				num.diff.sub<-num.diff.sub+length(setdiff(mutation.seq2[mutation.seq2$Type=="Substitution","Position"],mutation.seq1[mutation.seq1$Type=="Substitution","Position"]))
				
				num.diff.indel<-length(setdiff(mutation.seq1[mutation.seq1$Type!="Substitution","Position"],mutation.seq2[mutation.seq2$Type!="Substitution","Position"]))
				num.diff.indel<-num.diff.indel+length(setdiff(mutation.seq2[mutation.seq2$Type!="Substitution","Position"],mutation.seq1[mutation.seq1$Type!="Substitution","Position"]))
				
				num.diff.indel<- num.diff.indel*indel.penalty;
				
				pair.diff<-pair.diff+(num.diff.indel+num.diff.sub)/num_nt;
				#pair.diff<-
			}#inner for loop
		}#end of outer for loops
		
		return(pair.diff/num.pairs);
	}#end of function
    
 
	#function to go through
	intraClonalDiversities<-function(
			clones, #clone summary, 
			clone.assigns, #clone assignments for sequences
			mutations, #seq muations VRG
			Ig.RecSums, #Ig Recsum for VBases
			vlengths,  #total lengths of v 
			indel.penalty=1 #for penalty of indel.
		)
	{
		#clones<-df.clones.sp;clone.assigns<-df.cloneAssign.sp;mutations<-df.mutations.sp;Ig.RecSums<-IG.RecSum.sp;
		#vlengths<-df.vlength.sp;indel.penalty=1;
		
		#now go through each samples 
		samples<-unique(clones$sampleName)
		idis<-NULL;
		for(s in samples)
		{
			cat("doing sample ", s, "\n");
			flush.console();
			clones.sample<-clones[clones$sampleName==s,]
			clone.assigns.sample<-clone.assigns[clone.assigns$sampleName==s,]
			mutations.sample<-mutations[mutations$sampleName==s,]
			Ig.RecSums.sample<-Ig.RecSums[Ig.RecSums$sampleName==s,]
			vlengths.sample<-vlengths[vlengths$sampleName==s,]
			
			#get each clones in this sample
			clones.sample<-clones.sample[clones.sample$X.Members>1,]
			idis.sample<-data.frame(clones.sample$CloneID,sampleName=s, idi=rep(0, length(clones.sample$CloneID)),  XMembers=0, cloneSize=0);
			index<-0;
			for(j in clones.sample$CloneID)
			{
				index<-index+1;
				XMembers<-clones.sample[clones.sample$CloneID==j,"X.Members"];
				cloneSize<-clones.sample[clones.sample$CloneID==j,"cloneSize"];
				clone.assigns.clone<-clone.assigns.sample[clone.assigns.sample$CloneID==j,];
				mutations.clone<-mutations.sample[is.element(mutations.sample$ReadID, clone.assigns.clone$ReadID),]
				Ig.RecSums.clone<-Ig.RecSums.sample[is.element(Ig.RecSums.sample$UID, clone.assigns.clone$ReadID),]
				vlengths.clone<-vlengths.sample[is.element(vlengths.sample$ReadID, clone.assigns.clone$ReadID),]
				idi<-intraClonalDiversity(clone.assignment=clone.assigns.clone,
						mutations=mutations.clone, ig.recSum=Ig.RecSums.clone,
						vlengths=vlengths.clone, indel.penalty=1
						)
				idis.sample[index,c(3,4,5)]<-c(idi,XMembers,cloneSize )
			}#for loop each sample
			idis<-rbind(idis, idis.sample);
		}#end of outer for loops.
		
		return (idis)
	}#end of function.
#' 	idis<-intraClonalDiversities(clones=df.clones[df.clones$X.Members>5,], #clone summary, 
#'			clone.assigns=df.cloneAssign, #clone assignments for sequences
#'			mutations=df.mutations, #seq muations VRG
#'			Ig.RecSums=IG.RecSum, #Ig Recsum for VBases
#'			vlengths=df.vlength,  #total lengths of v 
#'			indel.penalty=1 #for penalty of indel.
#'			)   


#testing code
#vlength.S1.CID26[1,c("discPosStart", "discPosEnd")]<-c(NA, NA) 
