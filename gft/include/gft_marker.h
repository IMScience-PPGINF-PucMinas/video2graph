
#ifndef _GFT_MARKER_H_
#define _GFT_MARKER_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_set.h"
#include "gft_adjrel.h"
#include "gft_analysis.h"
#include "gft_morphology.h"


namespace gft{
  namespace Image32{
    
    float  MaxRadiusByErosion(sImage32 *bin);
    float  MaxObjRadiusByErosion(sImage32 *bin);
    float  MaxBkgRadiusByErosion(sImage32 *bin);
    
    sImage32 *ObjMaskByErosion(sImage32 *bin, float radius);
    sImage32 *BkgMaskByErosion(sImage32 *bin, float radius);

    int *GetMarkers(sImage32 *label, sAdjRel *A);
    
 
   } //end Image32 namespace
} //end gft namespace

 
#endif
