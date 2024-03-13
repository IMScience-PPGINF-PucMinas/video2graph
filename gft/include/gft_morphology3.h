
#ifndef _GFT_MORPHOLOGY3_H_
#define _GFT_MORPHOLOGY3_H_

#include "gft_scene8.h"
#include "gft_adjrel3.h"
#include "gft_set.h"
#include "gft_pqueue32.h"


namespace gft{
  namespace Scene8{

    sScene8 *Dilate(sScene8 *scn, sAdjRel3 *A);
    sScene8 *Erode(sScene8 *scn, sAdjRel3 *A);
    
    sScene8 *ErodeBin(sScene8 *bin, sSet **seed, float radius);
    
  } //end Scene8 namespace
} //end gft namespace

#endif
