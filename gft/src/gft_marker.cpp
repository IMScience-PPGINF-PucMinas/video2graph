
#include "gft_marker.h"


namespace gft{
  namespace Image32{

    float MaxRadiusByErosion(sImage32 *bin){
      sImage32 *bin1=NULL,*edt=NULL;
      sAdjRel *A=NULL;
      int Imax,Emax;

      A = gft::AdjRel::Circular(1.0);
      bin1 = Threshold(bin, 0, 0);
      edt = Mask2EDT(bin1, A, INTERIOR, INT_MAX, 0);
      Emax = GetMaxVal(edt);
      Destroy(&bin1);
      Destroy(&edt);
      
      bin1 = Threshold(bin, 1, INT_MAX);
      edt = Mask2EDT(bin1, A, INTERIOR, INT_MAX, 0);
      Imax = GetMaxVal(edt);
      Destroy(&bin1);
      Destroy(&edt);
      gft::AdjRel::Destroy(&A);
      
      return MIN(sqrtf(Imax),sqrtf(Emax));
    }


    float MaxObjRadiusByErosion(sImage32 *bin){
      sImage32 *bin1=NULL,*edt=NULL;
      sAdjRel *A=NULL;
      int Imax;
      
      A = gft::AdjRel::Circular(1.0);
      bin1 = Threshold(bin, 1, INT_MAX);
      edt = Mask2EDT(bin1, A, INTERIOR, INT_MAX, 0);
      Imax = GetMaxVal(edt);
      Destroy(&bin1);
      Destroy(&edt);
      gft::AdjRel::Destroy(&A);
      
      return (sqrtf(Imax));
    }


    float MaxBkgRadiusByErosion(sImage32 *bin){
      sImage32 *bin1=NULL,*edt=NULL;
      sAdjRel *A=NULL;
      int Emax;
      
      A = gft::AdjRel::Circular(1.0);
      bin1 = Threshold(bin, 0, 0);
      edt = Mask2EDT(bin1, A, INTERIOR, INT_MAX, 0);
      Emax = GetMaxVal(edt);
      Destroy(&bin1);
      Destroy(&edt);
      gft::AdjRel::Destroy(&A);
      
      return (sqrtf(Emax));
    }



    sImage32 *BkgMaskByErosion(sImage32 *bin, float radius){
      sSet *S=NULL;
      sImage32 *bin1=NULL,*bin2=NULL;
      
      bin1 = Threshold(bin, 0, 0);
      bin2 = ErodeBin(bin1, &S, radius);
      Destroy(&bin1);
      gft::Set::Destroy(&S);
      return bin2;
    }
    

    sImage32 *ObjMaskByErosion(sImage32 *bin, float radius){
      sSet *S=NULL;
      sImage32 *bin1=NULL,*bin2=NULL;
      
      bin1 = Threshold(bin, 1, INT_MAX);
      bin2 = ErodeBin(bin1, &S, radius);
      Destroy(&bin1);
      gft::Set::Destroy(&S);
      return bin2;
    }



    int *GetMarkers(sImage32 *label, sAdjRel *A){
      int *S = NULL;
      int p, q, i, k, nseeds;
      Pixel u,v;
      nseeds = 0;
      for(p=0; p<label->n; p++){
	if(label->data[p]!=NIL) nseeds++;
      }
      S = (int *)malloc((nseeds+1)*sizeof(int));
      if(S == NULL) return NULL;
      
      if(A == NULL){
	S[0] = nseeds;
	i = 1;
	for(p=0; p<label->n; p++){
	  if(label->data[p]!=NIL){
	    S[i] = p;
	    i++;
	  }
	}
      }
      else{
	i = 1;
	for(p=0; p<label->n; p++){
	  if(label->data[p]==NIL) continue;
	  u.x = p%label->ncols;
	  u.y = p/label->ncols;
	  for (k=1; k < A->n; k++){
	    v.x = u.x + A->dx[k];
	    v.y = u.y + A->dy[k];
	    if (gft::Image32::IsValidPixel(label, v.x, v.y)){
	      q = v.x + v.y*label->ncols;
	      if(label->data[q] != label->data[p]){
		S[i] = p;
		i++;
		break;
	      }
	    }
	  }
	}
	S = (int *)realloc(S, i*sizeof(int));
	if(S == NULL) return NULL;
	S[0] = i-1;
      }
      return S;
    }
    

    
  } /*end Image32 namespace*/
} /*end gft namespace*/
