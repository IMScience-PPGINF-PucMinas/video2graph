#ifndef _GFT_IMAGE32F_H_
#define _GFT_IMAGE32F_H_

#include "gft_common.h"
#include "gft_image32.h"

namespace gft{

  namespace Image32f{

    /**
     * It supports both linear and two-dimensional access 
     * (i.e., img->data[p] or img->array[y][x] for a pixel
     * (x,y) at address p=x+y*xsize).
     */
    struct sImage32f {
      float *data;
      float **array;
      int nrows; /* numero de linhas (altura) */
      int ncols; /* numero de colunas (largura) */
      int n;     /* numero de pixels */
      float dx;
      float dy;
    };

    /**
     * \brief A constructor.
     */
    sImage32f *Create(int ncols,int nrows);
    sImage32f *Create(sImage32f *img);
    sImage32f *Create(sImage32 *img);
    
    /**
     * \brief A destructor.
     */
    void    Destroy(sImage32f **img);

    /**
     * \brief A copy constructor.
     */
    sImage32f *Clone(sImage32f *img);
    sImage32f *Clone(sImage32 *img);

    float   GetMinVal(sImage32f *img);
    float   GetMaxVal(sImage32f *img);
    
    void    Set(sImage32f *img, float value);

    bool    IsValidPixel(sImage32f *img, int x, int y);


    sImage32f *AddFrame(sImage32f *img, int sz, float value);
    sImage32f *RemFrame(sImage32f *fimg, int sz);

    
  } //end Image32f namespace

  typedef Image32f::sImage32f sImage32f;



  namespace Image32{
    sImage32 *Clone(sImage32f *img, float multiplier);
  } //end Image32 namespace



} //end gft namespace



#endif

