
#ifndef _ADJREL3F_H_
#define _ADJREL3F_H_

#include "gft_common.h"

namespace gft{
  namespace AdjRel3f{

    struct sAdjRel3f {
      float *dx;
      float *dy;
      float *dz;
      int n;
    };
    
    sAdjRel3f *Create(int n);
    void       Destroy(sAdjRel3f **A);
    sAdjRel3f *Clone(sAdjRel3f *A);
    
    sAdjRel3f *ChangeOrientationToLPS(sAdjRel3f *A,
				      char *ori);

  } /*end AdjRel3f namespace*/

  typedef AdjRel3f::sAdjRel3f sAdjRel3f;

} /*end gft namespace*/

#endif

