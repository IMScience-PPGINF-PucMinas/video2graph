
#ifndef _GFT_RADIOMETRIC_H_
#define _GFT_RADIOMETRIC_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_curve.h"

namespace gft{

  namespace Image32{

    sCurve *Histogram(sImage32 *img);
    sCurve *NormHistogram(sImage32 *img);
    sCurve *NormalizeHistogram(sCurve *hist);
    sCurve *RemoveEmptyBins(sCurve *hist);
    
    sImage32 *LinearStretch(sImage32 *img,
			    int omin, int omax,
			    int nmin, int nmax);
    void LinearStretchinplace(sImage32 *img, 
			      int omin, int omax, 
			      int nmin, int nmax);
    
  } //end Image32 namespace
    

} //end gft namespace

#endif


