#include "skda.h"
#include "fmin.h"
double skdaCD_obj(double c, void *arg)
{const skdaCDinfo *info = (const skdaCDinfo *)arg;
  double ans;
  int p,j,k;
  int n;
  double *phat;
  const int *y;
double prtest=0.0;
   phat=(info->phat);
   y=(info->y);
   p=*(info->p);
   n=*(info->n);
   for(j=0; j<p; j++) {
       *((info->lam)+j)=(info->tau)*( c*(*((info->lamtmp1)+j)- *((info->lamtmp0)+j)) +(*((info->lamtmp0)+j)));
   }

  condSMXprob(info->n, info->p, info->nclass, info->y,  info->lam, info->dsqx, info->nsqp, info->np, info->ny, info->prop, info->phat);

          ans=0.0;
          for(k=0; k<n; k++) {
	     ans+=log( phat[k+(*((info->ny)+y[k])) ]     );
          }
  ans=-ans;
  return ans;
}

/*********************************************/


/*********************************************/



void coordescent(const int *n, const int *p, const int *nclass, const double *x,
                 const int *y, const double *tau,
                 const int *maxct, const double *prop, double *lam, double *phat) {

int i, j, k;
int lp;
double *lamtmp1, *lamtmp0;
double *templam;
double *oldlam;
double *bestlam;
int ct=1;
double ferr0;
double pe, peold;
double f0;
double temp;
double normerror;
double test=1.0;
const double tol = 1e-7;
double cmin;
double pinv=1.0/((double) (*p));
double prtest=0.0;
int *ny; /* stores info of n*y[k] */
double *dsqx;
int mybase;
int *nsqp; /* to store 0, n^2, 2n^2, 3n^3, ..., (p-1)n^2*/
int *np; /* to store 0, n, 2n, 3n, ..., (n-1)n */
 skdaCDinfo info;
 double wutempn;


     
  ny=(int*)Calloc((*nclass), int);
   k=0;
   for (j=0; j<(*nclass); j++) {
     ny[j]=k;
     k+=(*n);
   }
  np=(int*)Calloc((*n), int);
   k=0;
   for (j=0; j<(*n); j++) {
   np[j]=k;
   k+=(*n);
   }
   
nsqp=(int*)Calloc((*p),int);
  k=0;
  i=(*n)*(*n);
  for (j=0; j<(*p); j++) {
    nsqp[j]=k;
    k+=i;
  }
 
lamtmp1=(double*)Calloc((*p), double);
lamtmp0=(double*)Calloc((*p), double);
templam=(double*)Calloc((*p), double);
dsqx=(double*)Calloc((*p)*(*n)*(*n), double);
lp=0;
mybase=0;
for (j=0; j<(*p); j++) {
   for (i=0; i<(*n); i++) {
      for (k=0; k<(*n); k++) {
       temp=x[mybase+i]-x[mybase+k];
       dsqx[lp]=temp*temp;
       lp+=1;
      }
   }
   mybase+=(*n);
}
       


info.p=p;
info.n=n;
info.nclass=nclass;
info.y=y;
info.dsqx=dsqx;
info.lamtmp0=lamtmp0;
info.lamtmp1=lamtmp1;
info.lam=templam;
info.phat=phat;
info.tau=*tau;
info.ny=ny;
info.nsqp=nsqp;
info.np=np;
info.prop=prop;

 oldlam=(double*)Calloc((*p), double);
bestlam=(double*)Calloc((*p), double);


for (i=0; i<(*p); i++) {
    bestlam[i]=0.0;
}

ferr0=0.0;
for(i=0; i<(*p); i++) {
  ferr0+=lam[i];
}

pe=0.0;
for(i=0; i<(*p); i++) {
 if(lam[i]<pe) {
   pe=lam[i];
   }
}

if(fabs(ferr0-(*tau))>1e-3 || pe < 0.0 ) {
  temp=(*tau)*pinv;
 for (i=0; i<(*p); i++) {
    lam[i]=temp;
 }
} 

peold=-1e8;
for(ct=1; ct<4; ct++)  {

   for(j=0; j<(*p); j++) {
      oldlam[j]=lam[j];
   }

   for(i=0; i<(*p); i++) {
      for(j=0; j<(*p); j++) {
             lamtmp1[j]=lam[j];
             lamtmp0[j]=0.0;
       }
      lamtmp1[i]=0.0;
      lamtmp0[i]=1.0;
     temp=0.0;
     for(j=0; j<(*p); j++) {
          temp+=lamtmp1[j];
     }

     if(temp>1e-4) {
       temp=1.0/temp;
       for(j=0; j<(*p); j++) {
          lamtmp1[j]=lamtmp1[j]*temp;
       }
      cmin=0.55;
      temp=skdaCD_obj(cmin,(void *)(&info));


       wutempn=0.0;
       for(j=0; j<=10; j++) {
         f0=skdaCD_obj(wutempn,(void *)(&info));
         if(temp> f0)
         {  cmin=wutempn;
         }
       wutempn+=0.1;
       }
           
       for(k=0; k<(*p);k++) 
       { 
	 lam[k]=(*tau)*(cmin*(lamtmp1[k]-lamtmp0[k])+lamtmp0[k]);
       }

       } /* end if */
    }  /* end for (i) */

    
    condSMXprob(n, p, nclass, y, lam, dsqx, nsqp, np, ny, prop, phat);

     ferr0=0.0;
     for(i=0;i<(*p); i++) {
       if(fabs(lam[i]-bestlam[i])>ferr0) {
          ferr0= fabs(lam[i]-bestlam[i]);
       }
     }
     if(ferr0<1e-3) {
       break;
     }
     pe=0.0;
     for(k=0; k<(*n); k++) {
       pe+=log(phat[k+ny[y[k]]]);
     }
    


    if(peold<pe) {
      peold=pe;
      for(i=0; i<(*p); i++) {
      bestlam[i]=lam[i];
      }
    }
   test=0.0;
   for(i=0; i<(*p); i++) {
      if(fabs(lam[i]-oldlam[i])>1e-3)  {
      test=1.0;
      }
   }

} /* end for loop in CD*/



test=0.0;
ct=1;
peold=-1e8;
while(test>0.5) {

   for(j=0; j<(*p); j++) {
      oldlam[j]=lam[j];
   }

   for(i=0; i<(*p); i++) {
      for(j=0; j<(*p); j++) {
             lamtmp1[j]=lam[j];
             lamtmp0[j]=0.0;
       }
      lamtmp1[i]=0.0;
      lamtmp0[i]=1.0;
     temp=0.0;
     for(j=0; j<(*p); j++) {
          temp+=lamtmp1[j];
     }

     if(temp>1e-4) {
       temp=1.0/temp;
       for(j=0; j<(*p); j++) {
          lamtmp1[j]=lamtmp1[j]*temp;
       }
      cmin = pz_fmin(0.0, 1.0, &skdaCD_obj, (void *)(&info), tol);
       temp=skdaCD_obj(cmin,(void *)(&info));
       f0=skdaCD_obj(0.0,(void *)(&info));
       if(temp> f0)
       {  cmin=0.0;
       }

       f0=skdaCD_obj(1.0,(void *)(&info));
       if(temp> f0)
       {  cmin=1.0;
       }
    
       for(k=0; k<(*p);k++) 
       { 
	 lam[k]=(*tau)*(cmin*(lamtmp1[k]-lamtmp0[k])+lamtmp0[k]);
       }

       } /* end if */
    }  /* end for (i) */

    
    condSMXprob(n, p, nclass, y, lam, dsqx, nsqp, np, ny, prop, phat);

     ferr0=0.0;
     for(i=0;i<(*p); i++) {
       if(fabs(lam[i]-bestlam[i])>ferr0) {
          ferr0= fabs(lam[i]-bestlam[i]);
       }
     }
     if(ferr0<1e-3) {
       break;
     }
     pe=0.0;
     for(k=0; k<(*n); k++) {
       pe+=log(phat[k+ny[y[k]]]);
     }
    


    if(peold<pe) {
      peold=pe;
      for(i=0; i<(*p); i++) {
      bestlam[i]=lam[i];
      }
    }
   test=0.0;
   for(i=0; i<(*p); i++) {
      if(fabs(lam[i]-oldlam[i])>1e-3)  {
      test=1.0;
      }
   }


   ct+=1;
   if(ct>(*maxct)) {
  

     for(i=0; i<(*p); i++) {
      lam[i]=bestlam[i];
     }
     condSMXprob(n, p, nclass, y, lam, dsqx, nsqp, np, ny,  prop, phat);
    break;
   }
} /* end while*/



Free(lamtmp0); Free(lamtmp1); Free(templam); Free(oldlam);
Free(bestlam); Free(ny); Free(dsqx); Free(nsqp); Free(np); 

}




/**************************/

void condprob(const int *n, const int *p, const int *nclass,
              const int *y, const double *x, const double *lam,
              const int *m, const double *newx, const double *prop,
              double *phat)

{ 
    int i, j, k, w;
    double tmp;
    double tmpw;
    double tempprod;
    int *poslamind;
    int poslamct;
    double *temp;
 
double prtest=0.0;

    int *mk, *nk;


 mk=(int*)Calloc((*p), int);
 nk=(int*)Calloc((*p), int);

   i=0;
   j=0;
   for (k=0; k<(*p); k++) {
     nk[k]=i;
     mk[k]=j;
     i+=*n;
     j+=*m;
   }






   temp=(double*)Calloc((*nclass), double);
   poslamind=(int*)Calloc((*p), int);
 
   poslamct=0;
   for (i=0; i<(*p); i++) {
      if(lam[i]>1e-3) {
        poslamind[poslamct]=i;
        poslamct+=1;
      }
   }


    for(i=0; i<(*m); i++) 
    {
     for(k=0; k<(*nclass); k++)
      {temp[k]=0.0;
      }
	for(j=0; j<(*n); j++) 
	{
          tmpw=0.0;
	    for(w=0; w<poslamct; w++) 
	    {  k=poslamind[w];
 
	       tmp=lam[k]*(newx[i+mk[k]]-x[j+nk[k]]);
              tmpw+=tmp*tmp;
            }
          tempprod=exp(-tmpw);
         temp[y[j]]+=tempprod;


     	}
      tmpw=0.0;
      for(k=0; k<(*nclass); k++) 
       {temp[k]=temp[k]*prop[k];
       tmpw+=temp[k];
       }
       tmpw=1.0/tmpw;
     for(k=0; k<(*nclass); k++)
      {
      phat[i+mk[k]]=temp[k]*tmpw;
      }
    }




Free(poslamind);Free(temp);   Free(mk); Free(nk);

}



/**************************/
void condSMXprob(const int *n, const int *p, const int *nclass,
              const int *y,  const double *lam,
              const double *dsqx, const int *nsqp, const int *np,
              const int *ny, const double *prop,
              double *phat)
/* this is another version of condprob with x=newx */
{ 
    int i, j, k, w;
    double tmp;
    double tmpw;
    double tempprod;
    int *poslamind;
    int poslamct;
    double *temp;
double prtest=0.0;
double *lamsq;


   temp=(double*)Calloc((*nclass), double);
   poslamind=(int*)Calloc((*p), int);
   lamsq=(double*)Calloc((*p), double);
    
   poslamct=0;
   for (i=0; i<(*p); i++) {
      if(lam[i]>1e-3) {
        poslamind[poslamct]=i;
        poslamct+=1;
      }
      lamsq[i]=lam[i]*lam[i];
   }


 for(i=0; i<(*n); i++) 
    { 
     for(k=0; k<(*nclass); k++)
      {temp[k]=0.0;
      }
    
     for(j=0; j<(*n); j++) 
	     {
    
       tmpw=0.0;
	     for(w=0; w<poslamct; w++) 
	       {k=poslamind[w];
    
	       tmp=lamsq[k]*dsqx[nsqp[k]+np[i]+j];
              tmpw+=tmp;
            }
          tempprod=exp(-tmpw);
         temp[y[j]]+=tempprod;


     	}
      
      tmpw=0.0;
      for(k=0; k<(*nclass); k++) 
       { temp[k]=temp[k]*prop[k];
        tmpw+=temp[k];
       }
       tmpw=1.0/tmpw;
     for(k=0; k<(*nclass); k++)
      {
      phat[i+ny[k]]=temp[k]*tmpw;
      }
    }

Free(poslamind);Free(temp); Free(lamsq);

i=10;
}
  






