
#include "gft_radiometric.h"


namespace gft{

  namespace Image32{

    
    sCurve *Histogram(sImage32 *img){
      int i,p,n,nbins;
      sCurve *hist=NULL;
      
      nbins = GetMaxVal(img)+1;
      hist  = gft::Curve::Create(nbins);
      n     = img->ncols*img->nrows;
      for (p=0; p < n; p++)
	hist->Y[img->data[p]]++;
      for (i=0; i < nbins; i++) 
	hist->X[i] = i;
      
      return(hist);
    }


    sCurve *NormHistogram(sImage32 *img){
      int i,sum;
      sCurve *hist=NULL,*nhist=NULL;
      
      hist  = Histogram(img);
      sum   = img->ncols*img->nrows;
      nhist = gft::Curve::Create(hist->n);
      for (i=0; i < nhist->n;i++){
	nhist->Y[i] = hist->Y[i]/sum;
	nhist->X[i] = hist->X[i];
      }
      gft::Curve::Destroy(&hist);
      return(nhist);
    }


    sCurve *NormalizeHistogram(sCurve *hist){
      sCurve *nhist;
      double sum;
      int i;
      
      nhist = gft::Curve::Clone(hist);
      sum = 0.0;
      for(i=0; i<nhist->n; i++)
	sum += nhist->Y[i];
      for(i=0; i<nhist->n; i++)
	nhist->Y[i] /= sum;
      
      return (nhist);
    }


    sCurve  *RemoveEmptyBins(sCurve *hist){
      sCurve *C;
      int n,i,j;
      n = 0;
      for(i=0; i<hist->n; i++){
	if(hist->Y[i]>0.0) n++;
      }
      j = 0;
      C = gft::Curve::Create(n);
      for(i=0; i<hist->n; i++){
	if(hist->Y[i]>0.0){
	  C->Y[j] = hist->Y[i];
	  C->X[j] = hist->X[i];      
	  j++;
	}
      }
      return C;
    }
   

    sImage32 *LinearStretch(sImage32 *img,
			    int omin, int omax,
			    int nmin, int nmax){
      sImage32 *simg=NULL;
      int p,n;
      float a;
      
      simg = Create(img->ncols,img->nrows);
      n    = img->ncols*img->nrows;
      if (omin != omax) 
	a = (float)(nmax-nmin)/(float)(omax-omin);
      else
	a = INT_MAX;
      
      for (p=0; p < n; p++){
	if (img->data[p] < omin)
	  simg->data[p] = nmin;
	else 
	  if (img->data[p] > omax)
	    simg->data[p] = nmax;
	  else {
	    if (a != INT_MAX)	  
	      simg->data[p] = (int)(a*(img->data[p]-omin)+nmin);
	    else{
	      simg->data[p] = nmax;
	    }   
	  }
      }
      return(simg);
    }


    void LinearStretchinplace(gft::sImage32 *img, 
			      int omin, int omax, 
			      int nmin, int nmax){
      int p,n;
      float a;
      n = img->n;
      if (omin != omax) 
	a = (float)(nmax-nmin)/(float)(omax-omin);
      else
	a = INT_MAX;
      
      for (p=0; p < n; p++){
	if (img->data[p] < omin)
	  img->data[p] = nmin;
	else 
	  if (img->data[p] > omax)
	    img->data[p] = nmax;
	  else {
	    if (a != INT_MAX)	  
	      img->data[p] = (int)(a*(img->data[p]-omin)+nmin);
	    else{
	      img->data[p] = nmax;
	    }   
	  }
      }
    }

    

  } //end Image32 namespace


} //end gft namespace


