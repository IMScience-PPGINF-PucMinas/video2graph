
#ifndef _GFT_FILTERING_H_
#define _GFT_FILTERING_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_adjrel.h"
#include "gft_cimage.h"
#include "gft_cimage32f.h"

namespace gft{

  namespace Kernel{

    struct sKernel {
      float *val;
      sAdjRel *adj;
    };

    sKernel *Make(char *coefs);
    sKernel *Create(sAdjRel *A);
    sKernel *Clone(sKernel *K);
    sKernel *Normalize(sKernel *K);    
    void     Destroy(sKernel **K);

    sKernel *Gaussian(sAdjRel *A, float stddev);
    
  } //end Kernel namespace

  typedef Kernel::sKernel sKernel;
  
  namespace Image32{

    sImage32 *GaussianBlur(sImage32 *img, float stddev);
    sImage32 *GaussianBlur(sImage32 *img);
    
    sImage32 *SobelFilter(sImage32 *img);
    sImage32 *LinearFilter(sImage32 *img, sKernel *K);

    sImage32 *ImageMagnitude(sImage32 *imgx, sImage32 *imgy);

    sImage32 *MedianFilter(sImage32 *img, sAdjRel *A);

    void ModeFilterLabel(sImage32 *label, float r);
    
  } //end Image32 namespace


  namespace CImage{
    sImage32 *SobelFilter(sCImage *cimg);
    
  } //end CImage namespace    


  namespace CImage32f{
    sImage32f *SobelFilter(sCImage32f *cimg);
  } //end CImage32f namespace


  namespace Image32f{
    sImage32f *SobelFilter(sImage32f *cimg);
    sImage32f *LinearFilter(sImage32f *img, sKernel *K);
    sImage32f *ImageMagnitude(sImage32f *imgx, sImage32f *imgy);
  } //end Image32f namespace  

  
} //end gft namespace


#endif


