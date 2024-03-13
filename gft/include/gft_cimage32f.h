#ifndef _GFT_CIMAGE32F_H_
#define _GFT_CIMAGE32F_H_

#include "gft_image32f.h"
#include "gft_cimage.h"
#include "gft_color.h"

namespace gft{
  namespace CImage32f{

    struct sCImage32f {
      sImage32f *C[3];
    };

    sCImage32f *Create(int ncols, int nrows);
    sCImage32f *Create(sCImage32f *cimg);
    sCImage32f *Create(sImage32f *img);
    void       Destroy(sCImage32f **cimg);
    sCImage32f *Clone(sCImage32f *cimg);
    sCImage32f *Clone(sCImage *cimg);

    void    Set(sCImage32f *cimg, float r, float g, float b);
    
    sCImage32f *RGB2Lab(sCImage *cimg);
 
    sCImage32f *AddFrame(sCImage32f *cimg, int sz, float r, float g, float b);
    sCImage32f *RemFrame(sCImage32f *cimg, int sz);
    
  } //end CImage32f namespace

  typedef CImage32f::sCImage32f sCImage32f;

} //end gft namespace

#endif

