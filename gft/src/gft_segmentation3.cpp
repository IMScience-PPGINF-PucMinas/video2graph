
#include "gft_segmentation3.h"

namespace gft{

  namespace Scene8{

  sScene8   *Threshold(sScene8 *scn, int lower, int higher){
    sScene8 *bin=NULL;
    int p;
    bin = Create(scn);
    for(p=0; p<scn->n; p++){
      if((scn->data[p]>=lower)&&(scn->data[p]<=higher))
	bin->data[p]=1;
      else
	bin->data[p]=0;
    }
    bin->maxval = 1;
    return(bin);
  }


  sScene8 *GetBoundaries(sScene8 *scn, sAdjRel3 *A){
    sScene8 *hscn=NULL;
    int p,q,i;
    Voxel u,v;
    hscn = Create(scn);
    for(u.c.z=0; u.c.z<hscn->zsize; u.c.z++){
      for(u.c.y=0; u.c.y<hscn->ysize; u.c.y++){
	for(u.c.x=0; u.c.x<hscn->xsize; u.c.x++){
	  p = GetVoxelAddress(hscn,u);
	  if(scn->data[p] != 0){
	    for(i=1; i<A->n; i++){
	      v.v = u.v + A->d[i].v;
	      if(IsValidVoxel(hscn,v)){
		q = GetVoxelAddress(hscn,v);
		if(scn->data[p] != scn->data[q]){
		  hscn->data[p] = scn->data[p];
		  break;
		}
	      } 
	      else {
		hscn->data[p] = scn->data[p];
		break;
	      }
	    }
	  }
	}
      }
    }
    return(hscn);
  }


  sScene8 *GetTransitions(sScene8 *scn, sAdjRel3 *A){
    sScene8 *hscn=NULL;
    int p,q,i;
    Voxel u,v;
    hscn = Create(scn);
    for(u.c.z=0; u.c.z<hscn->zsize; u.c.z++){
      for(u.c.y=0; u.c.y<hscn->ysize; u.c.y++){
	for(u.c.x=0; u.c.x<hscn->xsize; u.c.x++){
	  p = GetVoxelAddress(hscn,u);
	  if(scn->data[p] != 0){
	    for(i=1; i<A->n; i++){
	      v.v = u.v + A->d[i].v;
	      if(IsValidVoxel(hscn,v)){
		q = GetVoxelAddress(hscn,v);
		if(scn->data[p] != scn->data[q]){
		  hscn->data[p] = scn->data[p];
		  break;
		}
	      } 
	    }
	  }
	}
      }
    }
    return(hscn);
  }

    

  } //end Scene8 namespace


  namespace Scene16{

  sScene8   *Threshold(sScene16 *scn, int lower, int higher){
    sScene8 *bin=NULL;
    int p;
    bin = Scene8::Create(scn->xsize,scn->ysize,scn->zsize);
    bin->dx = scn->dx;
    bin->dy = scn->dy;
    bin->dz = scn->dz;
    for(p=0; p<scn->n; p++){
      if((scn->data[p]>=lower)&&(scn->data[p]<=higher))
	bin->data[p]=1;
      else
	bin->data[p]=0;
    }
    bin->maxval = 1;
    return(bin);
  }

  } //end Scene16 namespace


  namespace Scene32{

  sScene8 *Threshold(sScene32 *scn, int lower, int higher){
    sScene8 *bin=NULL;
    int p;
    bin = Scene8::Create(scn->xsize,scn->ysize,scn->zsize);
    bin->dx = scn->dx;
    bin->dy = scn->dy;
    bin->dz = scn->dz;
    for(p=0; p<scn->n; p++){
      if((scn->data[p]>=lower)&&(scn->data[p]<=higher))
	bin->data[p]=1;
      else
	bin->data[p]=0;
    }
    bin->maxval = 1;
    return(bin);
  }


  int Otsu(gft::sScene32 *scn){
    gft::sCurve *hist=gft::Scene32::NormHistogram(scn);
    double p1,p2,m1,m2,s1,s2,J,Jmax=-1.0;
    int i,T,Topt=0,Imax=gft::Scene32::GetMaximumValue(scn);
    
    for (T=1; T < Imax; T++){
      p1 = 0.0;
      for (i=0; i <= T; i++) 
	p1 += hist->Y[i];
      p2 = 1.0 - p1;
      if ((p1 > 0.0)&&(p2 > 0.0)){
	m1 = 0.0;
	for (i=0; i <= T; i++) 
	  m1 += hist->Y[i]*i;
	m1 /= p1;
	m2 = 0.0;
	for (i=T+1; i <= Imax; i++) 
	  m2 += hist->Y[i]*i;
	m2 /= p2;
	s1 = 0.0;
	for (i=0; i <= T; i++) 
	  s1 += hist->Y[i]*(i-m1)*(i-m1);
	s1 /= p1;
	s2 = 0.0;
	for (i=T+1; i <= Imax; i++) 
	  s2 += hist->Y[i]*(i-m2)*(i-m2);
	s2 /= p2;
	J = (p1*p2*(m1-m2)*(m1-m2))/(p1*s1+p2*s2);
      }else{
	J = 0.0;      
      }
      if (J > Jmax){
	Jmax = J;
	Topt = T;
      }
    }
    //printf("Otsu: %d\n",Topt);
    gft::Curve::Destroy(&hist);
    return(Topt);
  }


    
  } //end Scene32 namespace



  namespace Scene{

  sScene8 *Threshold(sScene *scn, int lower, int higher){
    switch(scn->nbits){
    case  8:
      return Scene8::Threshold(scn->ptr.scn8, lower,higher);
    case 16:
      return Scene16::Threshold(scn->ptr.scn16, lower,higher);
    case 32:
      return Scene32::Threshold(scn->ptr.scn32, lower,higher);
    }
    return NULL;
  }


  } //end Scene namespace



} //end gft namespace

