
#ifndef _GFT_SEGMENTATION3_H_
#define _GFT_SEGMENTATION3_H_

#include "gft_common.h"
#include "gft_scene.h"
#include "gft_adjrel3.h"
#include "gft_radiometric3.h"

namespace gft{

  namespace Scene8{

    sScene8 *Threshold(sScene8 *scn, int lower, int higher);
    sScene8 *GetBoundaries(sScene8 *scn, sAdjRel3 *A);
    sScene8 *GetTransitions(sScene8 *scn, sAdjRel3 *A);

  } //end Scene8 namespace



  namespace Scene16{

    sScene8 *Threshold(sScene16 *scn, int lower, int higher);

  } //end Scene16 namespace



  namespace Scene32{

    sScene8 *Threshold(sScene32 *scn, int lower, int higher);
    int Otsu(sScene32 *scn);
    
  } //end Scene32 namespace


  namespace Scene{

    sScene8 *Threshold(sScene *scn, int lower, int higher);

  } //end Scene namespace


} //end gft namespace


#endif

