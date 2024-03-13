
#include "gft_image32f.h"

namespace gft{
  namespace Image32f{

    sImage32f *Create(sImage32f *img){
      sImage32f *nimg=NULL;
      nimg = Create(img->ncols, img->nrows);
      nimg->dx = img->dx;
      nimg->dy = img->dy;
      return nimg;
    }

    sImage32f *Create(sImage32 *img){
      sImage32f *nimg=NULL;
      nimg = Create(img->ncols, img->nrows);
      nimg->dx = img->dx;
      nimg->dy = img->dy;
      return nimg;
    }
    
    sImage32f *Create(int ncols, int nrows){
      sImage32f *img=NULL;
      int i;
      img = (sImage32f *) calloc(1,sizeof(sImage32f));
      if (img == NULL){
	gft::Error((char *)MSG1,(char *)"Image32f::Create");
      }
      img->data  = gft::AllocFloatArray(nrows*ncols);
      img->ncols = ncols;
      img->nrows = nrows;
      img->n = nrows*ncols;
      img->dx = img->dy = 1.0;
      
      img->array = (float**)malloc(nrows*sizeof(float*));
      if(img->array == NULL){
	gft::Error((char *)MSG1,(char *)"Image32f::Create");
      }
      for(i = 0; i < nrows; i++){
	img->array[i] = (img->data + i*ncols);
      }
      return(img);
    }
    

    void Destroy(sImage32f **img){
      sImage32f *aux;
      if(img != NULL){
	aux = *img;
	if (aux != NULL){
	  if(aux->data !=  NULL) gft::FreeFloatArray(&aux->data);
	  if(aux->array != NULL) free(aux->array);
	  free(aux);
	  *img = NULL;
	}
      }
    }
    
    
    sImage32f *Clone(sImage32f *img){
      sImage32f *imgc;
      imgc = Create(img->ncols,img->nrows);
      memcpy(imgc->data,img->data,img->ncols*img->nrows*sizeof(float));
      imgc->dx = img->dx;
      imgc->dy = img->dy;
      return(imgc);
    }


    sImage32f *Clone(sImage32 *img){
      sImage32f *imgc;
      int i;
      imgc = Create(img->ncols,img->nrows);
      for (i=0; i < img->n; i++)
	imgc->data[i] = (float)img->data[i];
      imgc->dx = img->dx;
      imgc->dy = img->dy;
      return(imgc);
    }

    
    void Set(sImage32f *img, float value){
      int i;
      for (i=0; i < img->n; i++)
	img->data[i] = value;
    }


    float GetMinVal(sImage32f *img){
      int i,n;
      float min;
      n = img->n;
      min = img->data[0];
      for (i=1; i < n; i++)
        if (img->data[i] < min)
	  min = img->data[i];
      
      return(min);
    }
    
    float GetMaxVal(sImage32f *img){
      int i,n;
      float max;
      n = img->n;
      max = img->data[0];
      for (i=1; i < n; i++)
        if (img->data[i] > max)
	  max = img->data[i];
      
      return(max);
    }
    
    
    bool IsValidPixel(sImage32f *img, int x, int y){
      return ((x >= 0)&&(x < img->ncols)&&
	      (y >= 0)&&(y < img->nrows));
    }

    
    sImage32f *AddFrame(sImage32f *img, int sz, float value){
      sImage32f *fimg;
      int y,nbytes,offset;
      float *dst,*src;
      
      fimg = Create(img->ncols+(2*sz), img->nrows+(2*sz));
      Set(fimg, value);
      nbytes = sizeof(float)*img->ncols;
      offset = sz + fimg->ncols*sz;
      for (y=0,src=img->data,dst=fimg->data+offset; y < img->nrows;y++,src+=img->ncols,dst+=fimg->ncols){
	memcpy(dst, src, nbytes);
      }
      return(fimg);
    }


    sImage32f *RemFrame(sImage32f *fimg, int sz){
      sImage32f *img;
      int y,nbytes,offset;
      float *dst,*src;
      
      img = Create(fimg->ncols-(2*sz), fimg->nrows-(2*sz));
      nbytes = sizeof(float)*img->ncols;
      offset = sz + fimg->ncols*sz;
      for (y=0,src=fimg->data+offset,dst=img->data; y < img->nrows;y++,src+=fimg->ncols,dst+=img->ncols){
	memcpy(dst, src, nbytes);
      }
      return(img);
    }

  } /*end Image32f namespace*/


  
  namespace Image32{

    sImage32 *Clone(sImage32f *img, float multiplier){
      sImage32 *imgc;
      int p;
      imgc = gft::Image32::Create(img->ncols, img->nrows);
      for(p = 0; p < img->n; p++){
	imgc->data[p] = ROUND(img->data[p]*multiplier);
      }
      return imgc;
    }
    
  } /*end Image32 namespace*/  

  
} /*end gft namespace*/

