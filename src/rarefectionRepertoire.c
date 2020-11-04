///////////////////////////////
// rarefctionRepertoire.c
// C Function File to do rarefction for repertoire size estiamtion
//
#include <R.h> // R memory io, in the standard list of system directories
#include <Rinternals.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

//used to be called by .C() 
void add_rc(double* x, double* y, double* o)
{
	Rprintf("feng in side c function\n");
	*o=*x+*y;
}
//this is c version of rarefection single run. this will be called by c only, and 
//not going to be called by R
//input: freq unsigned pointer, freq array, 
//              count unsigned pointer to count of the freq array
//          size_freq, unsigned size of freq 
//          T, unsigned total number of sampling units
//           t, unsigned subsampling size, t>0 and t<T
double rarefection_Inc_each_c(unsigned* freq, unsigned* count,  unsigned size_freq, unsigned T, unsigned t)
{
    //Rprintf("size_freq:%d\n", size_freq);
    //Rprintf("T:%d\n", T);
    //Rprintf("t:%d\n",t);
    double out=0;
    for(unsigned i =0; i<size_freq;i++)
    {
        //Rprintf("\tround i:%d\n", i );
        
        //we are here do the job.
        unsigned n=T-freq[i];
        unsigned d=T;
        double p=1;
        //Rprintf("\tn:%d, d:%d, p:%f\n", n,d,p);
        
        //skip if the freq[i] is too big
        if(T-freq[i]<t){
            //Rprintf("\tskipp.....\n");
            p=0;
        }
        else {
            for(unsigned j=0;j<t;j++)
            {
                //n-=j;
                //d-=j;
                p*=((double)(n-j))/(d-j); 
                //Rprintf("\t\tj:%d,==n:%d, d:%d, p:%f\n",j, n,d,p);
            }//end of inner loop: j
        }
        //Rprintf("\tp:%f\n", p);
        out+=(1-p)*count[i];
    }//end of out for loop: i
        return out;
}
//this is c version of rarefection run. this will be called by R wrapper 
//input: freq unsigned pointer, freq array, 
//              count unsigned pointer to count of the freq array
//          size_freq, unsigned size of freq 
//          T, unsigned total number of sampling units
//           t, unsigned pointer to subsampling size arrays, 
//
SEXP rarefection_Inc_c(SEXP freq /*unsigned* freq*/
                                                           , SEXP count /*unsigned* count*/
                                                           , SEXP size_freq /*unsigned size_freq*/
                                                           , SEXP T /*unsigned T*/
                                                           , SEXP ts /*unsigned t*/
                                                           , SEXP size_ts
                                                           )
{
    //now we need to parse the input 
    int* p_freq=INTEGER(freq);
    int* p_count=INTEGER(count);
    unsigned c_size_freq=asInteger(size_freq);
    unsigned c_T=asInteger(T);
    unsigned* p_t=INTEGER(ts);
    unsigned c_size_ts=asInteger(size_ts);
    
    //now let's do it.
    SEXP out=PROTECT(allocVector(REALSXP,  c_size_ts));
	double* c_out=REAL(out);//matrix(0, nrow=round(total), ncol=num)
		memset(c_out,0,c_size_ts*sizeof(double));
   for(int i=0;i<c_size_ts;i++)
    {
        c_out[i]=rarefection_Inc_each_c(p_freq, p_count, c_size_freq, c_T, p_t[i]);
//        c_out[1]=rarefection_Inc_each_c(p_freq, p_count, c_size_freq, c_T, 4);
        //Rprintf("the first one is :%f\n", c_out[1]);
    }
    UNPROTECT(1);
    return out;
    //return c_out[1];
}