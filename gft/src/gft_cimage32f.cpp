
#include "gft_cimage32f.h"

namespace gft{
  namespace CImage32f{


    sCImage32f *Create(int ncols, int nrows){
      sCImage32f *cimg=NULL;
      int i;
     
      cimg = (sCImage32f *) calloc(1, sizeof(sCImage32f));
      for (i=0; i < 3; i++)
	cimg->C[i] = Image32f::Create(ncols,nrows);
      return(cimg);
    }


    sCImage32f *Create(sCImage32f *cimg){
      return Create(cimg->C[0]->ncols, cimg->C[0]->nrows);
    }


    sCImage32f *Create(sImage32f *img){
      return Create(img->ncols, img->nrows);
    }


    void    Destroy(sCImage32f **cimg){
      sCImage32f *tmp;
      int i;
      
      tmp = *cimg;
      if (tmp != NULL) {
	for (i=0; i < 3; i++)
	  Image32f::Destroy(&(tmp->C[i]));
	free(tmp);
	*cimg = NULL;
      }
    }

    
    sCImage32f *Clone(sCImage32f *cimg){
      sCImage32f *imgc;
      int i;
      
      imgc = (sCImage32f *) calloc(1,sizeof(sCImage32f));
      if (imgc == NULL){
	gft::Error((char *)MSG1,(char *)"CImage32f::Clone");
      }
      for (i=0; i<3; i++)
	imgc->C[i] = Image32f::Clone(cimg->C[i]);
      return imgc;
    }


    sCImage32f *Clone(sCImage *cimg){
      sCImage32f *imgc;
      int i;
      
      imgc = (sCImage32f *) calloc(1,sizeof(sCImage32f));
      if (imgc == NULL){
	gft::Error((char *)MSG1,(char *)"CImage32f::Clone");
      }
      for (i=0; i<3; i++)
	imgc->C[i] = Image32f::Clone(cimg->C[i]);
      return imgc;
    }


    void    Set(sCImage32f *cimg, float r, float g, float b){
      Image32f::Set(cimg->C[0], r);
      Image32f::Set(cimg->C[1], g);
      Image32f::Set(cimg->C[2], b);
    }



    sCImage32f *RGB2Lab(sCImage *cimg){
      sCImage32f *cimg_lab;
      double l,a,b;
      int p,n;
      n = cimg->C[0]->n;
      cimg_lab = Create(cimg->C[0]->ncols, cimg->C[0]->nrows);

      //pragma omp parallel for private(l,a,b)
      for(p = 0; p < n; p++){
	gft::Color::RGB2Lab(cimg->C[0]->data[p],
			    cimg->C[1]->data[p],
			    cimg->C[2]->data[p],
			    l, a, b);
	cimg_lab->C[0]->data[p] = l;
	cimg_lab->C[1]->data[p] = a;
	cimg_lab->C[2]->data[p] = b;
      }
      return cimg_lab;
    }



    sCImage32f *AddFrame(sCImage32f *cimg, int sz, float r, float g, float b){
      sCImage32f *fcimg;
      fcimg = (sCImage32f *) calloc(1,sizeof(sCImage32f));
      if (fcimg == NULL)
	gft::Error((char *)MSG1,(char *)"CImage32f::AddFrame");
      fcimg->C[0] = gft::Image32f::AddFrame(cimg->C[0], sz, r);
      fcimg->C[1] = gft::Image32f::AddFrame(cimg->C[1], sz, g);
      fcimg->C[2] = gft::Image32f::AddFrame(cimg->C[2], sz, b);
      return fcimg;
    }

    
    sCImage32f *RemFrame(sCImage32f *cimg, int sz){
      sCImage32f *fcimg;
      fcimg = (sCImage32f *) calloc(1,sizeof(sCImage32f));
      if (fcimg == NULL)
	gft::Error((char *)MSG1,(char *)"CImage32f::AddFrame");
      fcimg->C[0] = gft::Image32f::RemFrame(cimg->C[0], sz);
      fcimg->C[1] = gft::Image32f::RemFrame(cimg->C[1], sz);
      fcimg->C[2] = gft::Image32f::RemFrame(cimg->C[2], sz);
      return fcimg;
    }
  
    
  } /*end CImage32f namespace*/
} /*end gft namespace*/
