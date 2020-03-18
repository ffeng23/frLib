#Code Name: pickBalancedBarcode.R
#	R code to pick balancedBarcode, meaning each position we have diversed 
#		nt's. 1) show have A/C/G/T; 2) A/C/G/T are equaly distributed.
#		Point 1 is a must and Point 2 is used to select the best ones that
#		have passed the first point.
#########################

#Require library(Biostrings)
#for R library import(Biostrings)
library(Biostrings)

#Require library xlsx
#for R library import(xlsx)
library(xlsx)
library(rcPkg)

###############################function definition
#functions the flat for loops like below
#
# without replacement
# order doesn't matter
#	for(i in 1:num)
#		for(j in i+1:num)
#			.....
#				for(k in j+1:num)

# s is the total number of elements . 
# num is the num of elements of the subset of combination
#		say we need a subset of 4 sequences, for example
subsetFullCombination<-function(s,num)
{
	if(missing(s)||missing(num))
	{
		stop("Please specify the input .....")
	}
	if(class(num)!="integer" && class(num)!="numeric")
	{
		stop("please specify an interge for num of permutations")
	}
	if(num>s)
	{
		stop(" the num of elements of the permutatoins has to be smaller than the length of the input array.")
	}
	index.array<-c(1:num); #the starting index sequence
	#sset<-list(); #keep the sequences
	#sset.index<-data.frame(); #keep the indexes
	len<-s
	count<-0;
	flag_finished<-FALSE;
	total<-factorial(len)/factorial(num)/factorial(len-num);
	#inside the while loop we keep checking/updating index.array
	#   index.array keeps the index for each of the num of elements.
	#	we check for it and determine what to increment the current index order
	#	go back up to incrementing the upper(previous) index
	
	#pre-allocate the space 
	sset.index<-matrix(0, nrow=round(total), ncol=num) 
	while(TRUE)
	{
		#need to make the first combination can be run!!!
		
		#let's do the job to pick one element for one subset.
		count<-count+1;
		cat("***count:", count,"\n")
		flush.console();
		sset.index[count,]<-index.array
		#sset[[count]]<-s[c(index.array)]
		
		if(count%%5000==0)
		{
			cat("subsetFullCombination():", count,"/", total, "......\n");
			flush.console();
		}
		cat("\tbefore updating......\n")
		flush.console();
		# updating condition
	
		for(i in num:1) #checking all indexing positions in the output array
		{
			if(i==num)
			{
				index.array[i]<-index.array[i]+1;
				if(index.array[i]<=len-(num-i))
				{
					#in this case, we simply increment ignore and go next round of update, skipping out of the for loop 
					#index.array[i]<-index.array[i]+1;
					break; 
				} #else we will update the next array index in the for loop
			} else { #need to check the previous one (big number index)
				if(index.array[i+1]>len-(num-(i+1))) #the big number one is "overflowed", we will keep incrementing this current one
				{
					index.array[i]<-index.array[i]+1; #update the current one
					
					if(index.array[i]<=len-(num-i))
					{
						#update the previouse index/ reset them to it "original"
						index.array[num:(i+1)]<-c((num-i):1)+index.array[i]
						#in this case, we simply increment ignore and go next round of update, skipping out of the for loop 
						#index.array[i]<-index.array[i]+1;
						break; 
					} else {#we will update the next array index in the for loop
						if(i==1)
						{
							#we are done;
							flag_finished<-TRUE;
						}
					}
				} else { #in this case, it is impossible. we simply do nothing 
					#index.array[i]<-index.array[i]+1;
					#break;
					cat("imppossible case\n")
				}
			}
			
		}#end of for loop
		cat("\tdone for one run\n")
		if(flag_finished)
		{
			break; #done!!
		}
	}
	cat("DONE!!!!\n")
	return(sset.index)
}
#check the diversity of the dnastrings at each position
#s is a dnastringset
diversity<-function(s)
{
	if(class(s)!="DNAStringSet")
	{
		stop("the input is not of DNAStringSet!!please check")
	}
	freq<-rep(0,4)
	names(freq)<-c("A","C","T","G")
	
	len<-length(s[[1]])
	diverse<-rep(TRUE,len);
	diverge<-0
	x<-toupper(as.matrix(s))
	for(i in 1:len)
	{
		diverse[i]<-setequal(unique(x[,i]), names(freq))
		for(j in 1:length(freq))
		{
			freq[j]<-sum(x[,i]==names(freq)[j])
		}
		#now get the freq
		freq<-freq/sum(freq)
		freq<-freq-c(0.25,0.25,0.25,0.25)
		cat("i:", sum(freq^2), "\n");
		diverge<-diverge+sum(freq^2)
	}
	
	list(diverse=diverse, freq=diverge/len)
}
subsetCombination4<-function(s)
{
	ret<-list();
	if(length(s)<4)
	{
		stop("the number of elements in the DNAStringSet is less than 4. please add more sequences in the set")
	}
	#if(length(s)==4)
	#{
	#	cat("4")
	#	return (list(s))
	#}
	len<-length(s)
	index<-0;
	for(i in 1:(len-3))
	{
		cat("i:",i, "\n")
		for(j in (i+1):(len-2))
		{
			#cat("\tj:",j, "\n")
			for(k in (j+1):(len-1))
			{
				#cat("\t\tk:",k, "\n")
				for(m in (k+1):(len))
				{
					#cat("\t\t\tm:",m, "\n")
					index<-index+1;
					ret[[index]]<-s[c(i,j,k,m)]
				}
			}
		}
	}
	return (ret)
}

getDiversedSet4<-function(s, include, s.pair, include.pair )
{
	if(missing(s))
	{
		stop("please specify input....")
	}
	
	#if(!missing(include))
	#{
	#	s<-append(s, include);
	#}
	#if(!missing(include.pair)&&!missing(s.pair))
	#{
	#	s.pair<-append(s.pair,include.pair)
	#}
	s.com<-s
	if(!missing(s.pair))
	{#combine both Reads
		s.com<-DNAStringSet(paste0(as.character(s), as.character(s.pair)))
	}
	st<-subsetCombination4(s.com)
	st.diverse<-list();
	diverge.score<-c();
	
	s.com.include<-NULL;
	if(!missing(include)&&!missing(include.pair))
	{
		s.com.include<-DNAStringSet(paste0(as.character(include), as.character(include.pair)))
	}
	
	for(i in 1:length(st))
	{
		#get diversity
		#cat("i:", i, "\n")
		score<-NULL;
		if(!missing(s.com.include))
		{
			score<-diversity(append(st[[i]], s.com.include))
		} else {
			score<-diversity(st[[i]]);
		}
		if(all(score[[1]]))
		{
			cat("\tgot one, i:", i, "\n")
			diverge.score<-c(diverge.score, score[[2]])
			st.diverse[[length(st.diverse)+1]]<-st[[i]]
		}
	}
	st.diverse[order(diverge.score, decreasing=TRUE)];

	#save the data for using.
	#setwd("E:\\feng\\LAB\\hg\\frLib\\dev")
	#save(s.com.include, st.diverse, file="balanced_umi_10nt_4Plus4.RData")
	# 

	return (st.diverse)
}

getDiverseSet<-function(s, num, include, s.pair, include.pair )
{
	if(missing(s))
	{
		stop("please specify input \"s\"....")
	}

	if(missing(num))
	{
		stop("please specify input \"num\"....")
	}
	#if(!missing(include))
	#{
	#	s<-append(s, include);
	#}
	#if(!missing(include.pair)&&!missing(s.pair))
	#{
	#	s.pair<-append(s.pair,include.pair)
	#}
	s.com<-s
	if(!missing(s.pair))
	{#combine both Reads
		s.com<-DNAStringSet(paste0(as.character(s), as.character(s.pair)))
		names(s.com)<-names(s);
	}
	#
	cat("generating the permutated sequence set.......\n")
	st0<-subsetFullCombination(s.com, num)
	st<-st0[[1]]
	st.index<-st0[[2]]
	
	st.diverse<-list();
	diverge.score<-c();
	
	s.com.include<-NULL;
	if(!missing(include)&&!missing(include.pair))
	{
		s.com.include<-DNAStringSet(paste0(as.character(include), as.character(include.pair)))
		names(s.com.include)<-names(include);
	}
	
	index<-c();
	cat("Find the balanced sets..")
	for(i in 1:length(st))
	{
		#get diversity
		if((i/length(st)*100)%%10==0)
		{
			cat("getDiverseSet():", i/length(st)*100, "%\n")
			flush.console();
		}
		score<-NULL;
		if(!missing(s.com.include))
		{
			score<-diversity(append(st[[i]], s.com.include))
		} else {
			score<-diversity(st[[i]]);
		}
		if(all(score[[1]]))
		{
			cat("\tgetDiverseSet():find one diversed set\n");
			flush.console();
			diverge.score<-c(diverge.score, score[[2]])
			st.diverse[[length(st.diverse)+1]]<-st[[i]]
			index<-c(index,i)
		}
	}
	cat("Done!\n")
	if(length(diverge.score)==0){
		return(list())
	}
	st.diverse<-st.diverse[order(diverge.score, decreasing=F)];
	index<-index[order(diverge.score, decreasing=F)];
	
	#save the data for using.
	#setwd("E:\\feng\\LAB\\hg\\frLib\\dev")
	#save(s.com.include, st.diverse, file="balanced_umi_10nt_4Plus4.RData")
	# 

	return (list(diverse=st.diverse, index=st.index[index,]))
}
############################

#Require library(Biostrings)
#for R library import(Biostrings)
library(Biostrings)

#Require library xlsx
#for R library import(xlsx)
library(xlsx)
library(rcPkg)

#now read the data file
setwd("E:\\feng\\LAB\\order")

file.name<-"Feng Feng - SO 15060176_nexteraUniqueBarcodes96.xlsx"

uniBarcodes<-read.xlsx(file.name, 1)

#now we need to turn the barcodes into fasta
barcodes<-DNAStringSet(uniBarcodes[,6])
names(barcodes)<-uniBarcodes[,8]

#split into 2 set, R1 and R2
barcodes.R1<-barcodes[1:96]
barcodes.R2<-barcodes[97:192]

	#brutal force evaluating diversity
	num.subset<-8
	
	#We are doing subset without replacement and order doesn't matter
	#


#to include
index.include<-c("H9", "H10", "H11", "H12")

include<-barcodes.R1[index.include]

excluded<-setdiff(barcodes.R1, index.include)

include.pair<-barcodes.R2[index.include]

s<-excluded
s.pair<-setdiff(barcodes.R2, index.include)

#for testing 10 subsetCombination 
s<-barcodes.R1[1:20]
include<-barcodes.R1[93:96]
s.pair<-barcodes.R2[1:20]
include.pair<-barcodes.R2[93:96]
ms<-as.matrix(s);
ind<-subsetFullCombination_w(16,15);

getDiverseSet_w(ms, ind)

#start running..........
s20_4<-getDiverseSet(s=s, s.pair=s.pair, num=4)
s96_4<-getDiverseSet(s=barcodes.R1, s.pair=barcodes.R2, num=4)

s96_6<-getDiverseSet(s=barcodes.R1, s.pair=barcodes.R2, num=6)