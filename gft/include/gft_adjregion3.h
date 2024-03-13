
#ifndef _GFT_ADJREGION3_H_
#define _GFT_ADJREGION3_H_

#include "gft_common.h"
#include "gft_scene.h"
#include "gft_adjrel3.h"

namespace gft{
  namespace AdjRegion3{

    union sDisplacement3 {
      v4si v;
      int  data[4];
      struct{ int x,y,z; } axis;
    };

    struct sAdjRegion3 {
      sDisplacement3 *d;
      int n;
      int *dp;

      //private:
      int xsize;
      int xysize;
      int max[3];
      int min[3];
    };


    sAdjRegion3 *Create(int n);
    sAdjRegion3 *Create(sAdjRel3 *A);
    sAdjRegion3 *Create(sScene8 *mask, Voxel Ref);
    void        Destroy(sAdjRegion3 **adjreg);
    sAdjRegion3 *Clone(sAdjRegion3 *adjreg);
    sAdjRegion3 *Merge(sAdjRegion3 *r1, sAdjRegion3 *r2);

    sScene8 *Export2Mask(sAdjRegion3 *adjreg);

    void    Draw(sAdjRegion3 *adjreg,
		 sScene8 *scn,
		 Voxel u,
		 uchar val);
    void    Draw(sAdjRegion3 *adjreg,
		 sScene16 *scn,
		 Voxel u,
		 ushort val);
    void    Draw(sAdjRegion3 *adjreg,
		 sScene32 *scn,
		 Voxel u,
		 int val);
    /*
    void    MT_Draw(sAdjRegion3 *adjreg,
		    sScene8 *scn,
		    Voxel u,
		    uchar val);
    */
    void    DrawOpt(sAdjRegion3 *adjreg,
		    sScene8 *scn,
		    int p, uchar val);
    void    DrawOpt(sAdjRegion3 *adjreg,
		    sScene16 *scn,
		    int p, ushort val);
    void    DrawOpt(sAdjRegion3 *adjreg,
		    sScene32 *scn,
		    int p, int val);
    /*
    void    MT_DrawOpt(sAdjRegion3 *adjreg,
		       sScene8 *scn,
		       int p, uchar val);
    */
    void    Optimize(sAdjRegion3 *adjreg,
		     int xsize, int ysize);
    void    Optimize(sAdjRegion3 *adjreg,
		     sScene8 *scn);
    void    Optimize(sAdjRegion3 *adjreg,
		     sScene16 *scn);
    void    Optimize(sAdjRegion3 *adjreg,
		     sScene32 *scn);

    void    RefreshLimits(sAdjRegion3 *adjreg);
    void    GetLimits(sAdjRegion3 *adjreg,
		      int *dx_min, int *dy_min, int *dz_min,
		      int *dx_max, int *dy_max, int *dz_max);

    bool    FitInside(sAdjRegion3 *adjreg, Voxel vx,
		      int xsize, int ysize, int zsize, int sz);
    bool    FitInside(sAdjRegion3 *adjreg, Voxel vx,
		      sScene8 *scn, int sz);
    bool    FitInside(sAdjRegion3 *adjreg, Voxel vx,
		      sScene16 *scn, int sz);
    bool    FitInside(sAdjRegion3 *adjreg, Voxel vx,
		      sScene32 *scn, int sz);

    float   InnerMean(sAdjRegion3 *adjreg,
		      Voxel vx,
		      sScene16 *scn);

    float   InnerSum(sAdjRegion3 *adjreg,
		     Voxel vx,
		     sScene16 *scn);


  } //end AdjRegion3 namespace

  typedef AdjRegion3::sAdjRegion3 sAdjRegion3;

} //end gft namespace
    
#endif

