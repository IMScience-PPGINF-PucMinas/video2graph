
#ifndef _GFT_ADJREL3_H_
#define _GFT_ADJREL3_H_

#include "gft_common.h"

namespace gft{
  namespace AdjRel3{

    union sDisplacement3 {
      v4si v;
      int  data[4];
      struct{ int x,y,z; } axis;
    };

    struct sAdjRel3 {
      sDisplacement3 *d;
      int n;
    };
    
    struct sAdjVxl {
      int *dp;
      int n;
    };


    sAdjRel3 *Create(int n);
    void      Destroy(sAdjRel3 **A);
    sAdjRel3 *Clone(sAdjRel3 *A);

    sAdjRel3 *Spheric(float r);
    sAdjRel3 *Ellipsoid(float rx, float ry, float rz);
    sAdjRel3 *SphericalShell(float inner_radius,
			     float outer_radius);
    sAdjRel3 *Box(int xsize, int ysize, int zsize);


    void     Scale(sAdjRel3 *A,
		   float Sx, float Sy, float Sz);
    void     ClipX(sAdjRel3 *A, int lower, int higher);
    void     ClipY(sAdjRel3 *A, int lower, int higher);
    void     ClipZ(sAdjRel3 *A, int lower, int higher);

    int      GetFrameSize(sAdjRel3 *A);
    void     GetFrameSize(sAdjRel3 *A, 
			  int *sz_x, 
			  int *sz_y, 
			  int *sz_z);
    float   *GetDistanceArray(sAdjRel3 *A);

    void     DestroyAdjVxl(sAdjVxl **N);
    
  } //end AdjRel3 namespace

  typedef AdjRel3::sAdjRel3 sAdjRel3;
  typedef AdjRel3::sAdjVxl sAdjVxl;

} //end gft namespace


#include "gft_scene.h"

namespace gft{
  namespace AdjRel3{

    sAdjVxl  *AdjVoxels(sAdjRel3 *A, int xsize, int ysize);
    sAdjVxl  *AdjVoxels(sAdjRel3 *A, sScene32 *scn);
    sAdjVxl  *AdjVoxels(sAdjRel3 *A, sScene16 *scn);
    sAdjVxl  *AdjVoxels(sAdjRel3 *A, sScene8 *scn);
    sAdjVxl  *AdjVoxels(sAdjRel3 *A, sScene *scn);

    sAdjRel3 *Spherical_mm(sScene32 *scn, float r_mm);
    sAdjRel3 *SphericalGrid_mm(float dx, float dy, float dz,
			       float r_mm, float spacement_mm);
    sAdjRel3 *SphericalGrid_mm(sScene32 *scn, 
			       float r_mm, float spacement_mm);

    float   *GetDistanceArray_mm(sAdjRel3 *A, sScene32 *scn);
   
  } //end AdjRel3 namespace
} //end gft namespace

#endif

