#ifndef _GFT_GEOMETRY_H_
#define _GFT_GEOMETRY_H_

#include "gft_common.h"

namespace gft{

  namespace Point{
    
    struct sPoint{
      float x;
      float y;
      float z;
    };
    
  } //end Point namespace
  
  typedef Point::sPoint sPoint;

  
  namespace Vector{

    struct sVector{
      float x;
      float y;
      float z;
    };

    float    ScalarProd(sVector v1, sVector v2);
    sVector  VectorProd(sVector v1, sVector v2);
    sVector  Rotate(sVector v1, float R[4][4]);
    void     Normalize(sVector *v);
    float    Magnitude(sVector *v);

  } //end Vector namespace

  typedef Vector::sVector sVector;

} //end gft namespace
    
    
#endif

