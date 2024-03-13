
#include "gft_geometry.h"

namespace gft{
  namespace Vector{

    float ScalarProd(sVector v1, sVector v2){
      return(v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
    }

    sVector VectorProd(sVector v1, sVector v2){
      sVector v;
      
      v.x = v1.y*v2.z - v1.z*v2.y;
      v.y = v1.z*v2.x - v1.x*v2.z;
      v.z = v1.x*v2.y - v1.y*v2.x;
      
      return(v);
    }
    
    sVector Rotate(sVector v1, float R[4][4]) {
      sVector v;
      
      v.x = v1.x*R[0][0] + v1.y*R[0][1] + v1.z*R[0][2];
      v.y = v1.x*R[1][0] + v1.y*R[1][1] + v1.z*R[1][2];
      v.z = v1.x*R[2][0] + v1.y*R[2][1] + v1.z*R[2][2];
      
      return(v);
    }
    
    void Normalize(sVector *v) {
      float norm;
      
      norm= sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
      
      if (norm != 0.0) {
	v->x = v->x/norm;
	v->y = v->y/norm;
	v->z = v->z/norm;
      }
    }

    float Magnitude(sVector *v) {
      return (sqrt(v->x*v->x + v->y*v->y + v->z*v->z));
    }


  } //end Vector namespace
} //end gft namespace

