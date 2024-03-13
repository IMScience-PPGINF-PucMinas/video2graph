
#include "gft_radiometric3.h"


namespace gft{

  namespace Scene32{

    sCurve *NormHistogram(sScene32 *scn){
      int i, sum;
      sCurve *nhist;
      
      nhist = Histogram(scn);
      sum = scn->n;
      for (i = 0; i < nhist->n; i++){
	nhist->Y[i] /= sum;
	nhist->X[i]=i;
      }
      return (nhist);
    }

    
    sCurve *Histogram(sScene32 *scn){
      int i, n, nbins;
      sCurve *hist = NULL;

      nbins = GetMaximumValue(scn)+1;
      hist  = gft::Curve::Create(nbins);
      n = scn->n;
      for (i = 0; i < n; i++)
	hist->Y[scn->data[i]]++;
      for (i = 0; i < nbins; i++)
	hist->X[i] = i;
     
      return (hist);
    }

    
    sCurve *Histogram(sScene32 *scn, int binwidth){
      int i,nbins,maxbins,b;
      sCurve *hist = NULL;

      maxbins = GetMaximumValue(scn)+1;

      nbins = maxbins/binwidth;
      if(maxbins%binwidth!=0) nbins+=1;

      hist  = Curve::Create(nbins);

      for(i=0; i<scn->n; i++){
	b = scn->data[i]/binwidth;
	hist->Y[b]++;
      }
      for(i=0; i<nbins; i++)
	hist->X[i] = i*binwidth + (float)(binwidth-1)/2.0;
      return (hist);
    }


    void LinearStretchinplace(sScene32 *scn, 
			      int omin,int omax,
			      int nmin,int nmax){
      float tmp,d,dn;
      int   p;

      d  = (float)(omax - omin);
      dn = (float)(nmax - nmin);
      for(p=0; p<scn->n; p++){
	if(scn->data[p] < omin)
	  scn->data[p] = nmin;
	else if(scn->data[p] > omax)
	  scn->data[p] = nmax;
	else if( d <= 0.0 ) scn->data[p] = nmin;
	else{
	  tmp = ((float)(scn->data[p] - omin))/d;
	  scn->data[p] = nmin + ROUND(tmp*dn);
	}
      }
      scn->maxval = nmax;
    }


    sScene32 *LinearStretch(sScene32 *scn, 
			    int omin,int omax,
			    int nmin,int nmax){
      sScene32 *sscn;
      sscn = Clone(scn);
      LinearStretchinplace(sscn, omin,omax, nmin,nmax);
      return sscn;
    }


  } //end Scene32 namespace


  namespace Scene16{

    sCurve *Histogram(sScene16 *scn, int binwidth){
      int i,nbins,maxbins,b;
      sCurve *hist = NULL;
      
      maxbins = GetMaximumValue(scn)+1;
      
      nbins = maxbins/binwidth;
      if(maxbins%binwidth!=0) nbins+=1;
      
      hist  = Curve::Create(nbins);
      
      for(i=0; i<scn->n; i++){
	b = scn->data[i]/binwidth;
	hist->Y[b]++;
      }
      for(i=0; i<nbins; i++)
	hist->X[i] = i*binwidth + (float)(binwidth-1)/2.0;
      return (hist);
    }


    void LinearStretchinplace(sScene16 *scn, 
			      int omin,int omax,
			      int nmin,int nmax){
      float tmp,d,dn;
      int   p;
      
      d  = (float)(omax - omin);
      dn = (float)(nmax - nmin);
      for(p=0; p<scn->n; p++){
	if(scn->data[p] < omin)
	  scn->data[p] = nmin;
	else if(scn->data[p] > omax)
	  scn->data[p] = nmax;
	else if( d <= 0.0 ) scn->data[p] = nmin;
	else{
	  tmp = ((float)(scn->data[p] - omin))/d;
	  scn->data[p] = nmin + ROUND(tmp*dn);
	}
      }
      scn->maxval = nmax;
    }


    sScene16 *LinearStretch(sScene16 *scn, 
			    int omin,int omax,
			    int nmin,int nmax){
      sScene16 *sscn;
      sscn = Clone(scn);
      LinearStretchinplace(sscn, omin,omax, nmin,nmax);
      return sscn;
    }

  } //end Scene16 namespace


  namespace Scene8{

    sCurve *Histogram(sScene8 *scn, int binwidth){
      int i,nbins,maxbins,b;
      sCurve *hist = NULL;

      maxbins = GetMaximumValue(scn)+1;

      nbins = maxbins/binwidth;
      if(maxbins%binwidth!=0) nbins+=1;

      hist  = Curve::Create(nbins);

      for(i=0; i<scn->n; i++){
	b = scn->data[i]/binwidth;
	hist->Y[b]++;
      }
      for(i=0; i<nbins; i++)
	hist->X[i] = i*binwidth + (float)(binwidth-1)/2.0;
      return (hist);
    }

    void LinearStretchinplace(sScene8 *scn, 
			      int omin,int omax,
			      int nmin,int nmax){
      float tmp,d,dn;
      int   p;
      
      d  = (float)(omax - omin);
      dn = (float)(nmax - nmin);
      for(p=0; p<scn->n; p++){
	if(scn->data[p] < omin)
	  scn->data[p] = nmin;
	else if(scn->data[p] > omax)
	  scn->data[p] = nmax;
	else if( d <= 0.0 ) scn->data[p] = nmin;
	else{
	  tmp = ((float)(scn->data[p] - omin))/d;
	  scn->data[p] = nmin + ROUND(tmp*dn);
	}
      }
      scn->maxval = nmax;
    }


    sScene8 *LinearStretch(sScene8 *scn, 
			   int omin,int omax,
			   int nmin,int nmax){
      sScene8 *sscn;
      sscn = Clone(scn);
      LinearStretchinplace(sscn, omin,omax, nmin,nmax);
      return sscn;
    }

  } //end Scene8 namespace


  namespace Scene{

    sCurve *Histogram(sScene *scn, int binwidth){
      switch(scn->nbits){
      case  8:
	return Scene8::Histogram(scn->ptr.scn8, binwidth);
      case 16:
	return Scene16::Histogram(scn->ptr.scn16, binwidth);
      case 32:
	return Scene32::Histogram(scn->ptr.scn32, binwidth);
      }
      return NULL;
    }


    void LinearStretchinplace(sScene *scn, 
			      int omin,int omax,
			      int nmin,int nmax){
      switch(scn->nbits){
      case  8:
	Scene8::LinearStretchinplace(scn->ptr.scn8,omin,omax,nmin,nmax);
	break;
      case 16:
	Scene16::LinearStretchinplace(scn->ptr.scn16,omin,omax,nmin,nmax);
	break;
      case 32:
	Scene32::LinearStretchinplace(scn->ptr.scn32,omin,omax,nmin,nmax);
	break;
      }
    }


    sScene *LinearStretch(sScene *scn, 
			  int omin,int omax,
			  int nmin,int nmax){
      sScene *sscn;
      sscn = Clone(scn);
      LinearStretchinplace(sscn, omin,omax, nmin,nmax);
      return sscn;
    }


  } //end Scene namespace


} //end gft namespace




