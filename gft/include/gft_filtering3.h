
#ifndef _GFT_FILTERING3_H_
#define _GFT_FILTERING3_H_

#include "gft_common.h"
#include "gft_scene.h"
#include "gft_adjrel3.h"
#include "gft_curve.h"
#include "gft_radiometric3.h"

namespace gft{

  namespace Kernel3{

    struct sKernel3 {
      float *val;
      sAdjRel3 *adj;
      int xsize,ysize,zsize;
    };

    sKernel3 *Create(sAdjRel3 *A);
    sKernel3 *Clone(sKernel3 *K);
    void      Destroy(sKernel3 **K);

    sKernel3 *Normalize(sKernel3 *K);

    sKernel3 *SphericalGaussian(float R, float s, float f);

  } //end Kernel3 namespace

  typedef Kernel3::sKernel3 sKernel3;


  
  namespace Scene8{

    void ModeFilterLabel(sScene8 *label, float r);

  } //end Scene8 namespace


  namespace Scene32{
    
    void ModeFilterLabel(sScene32 *label, float r);

    sScene32  *AccAbsDiff(sScene32 *scn, float r);

    void   SuppressHighIntensities(sScene32 *scn);

    sScene32 *Convolution(sScene32 *scn, sKernel3 *K);
    sScene32 *OptConvolution(sScene32 *scn, sKernel3 *K);

    sScene32 *GaussianBlur(sScene32 *scn);
    sScene32 *OptGaussianBlur(sScene32 *scn);
    //sScene32 *FastGaussianBlur(sScene32 *scn);
    //sScene32 *FastOptGaussianBlur(sScene32 *scn);
 
    sScene32 *Subsampling(sScene32 *scn);

    //-----------------

    sScene32  *LaplacianFilter(sScene32 *orig);
    sScene32  *SobelFilter(sScene32 *scn);
    sScene32  *SphericalGradient(sScene32 *scn, float r);    
    
  } //end Scene32 namespace


} //end gft namespace


#endif






