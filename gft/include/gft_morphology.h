
#ifndef _GFT_MORPHOLOGY_H_
#define _GFT_MORPHOLOGY_H_

#include "gft_image32.h"
#include "gft_adjrel.h"
#include "gft_set.h"
#include "gft_pqueue32.h"


namespace gft{
  namespace Image32{

    sImage32 *Dilate(sImage32 *img, sAdjRel *A);
    sImage32 *Erode(sImage32 *img, sAdjRel *A);

    sImage32 *ErodeBin(sImage32 *bin, sSet **seed, float radius);
    
    sImage32 *MorphGrad(sImage32 *img, sAdjRel *A);

    void SupRec_Watershed(sAdjRel *A, 
			  sImage32 *I, sImage32 *J, 
			  sImage32 *L, sImage32 *V);

    /*Removes all background connected components from the stack
      of binary images of I whose area (number of pixels) is <=
      a threshold and outputs a simplified image.*/
    sImage32 *AreaClosing(sAdjRel *A, sImage32 *I, int T);

    sImage32 *VolumeClosing(sAdjRel *A, sImage32 *I, int T);
    
    sImage32 *CloseHoles(sImage32 *img);
    
  } //end Image32 namespace
} //end gft namespace

#endif
