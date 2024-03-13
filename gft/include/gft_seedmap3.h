
#ifndef _GFT_SEEDMAP3_H_
#define _GFT_SEEDMAP3_H_

#include "gft_common.h"
#include "gft_scene.h"
#include "gft_adjregion3.h"
#include "gft_segmentation3.h"
#include "gft_fuzzycloud3.h"
#include "gft_gradient3.h"


namespace gft{
  /**
   * \brief Data structures for fast accessing seed pixels of a fuzzy model (e.g., CSM).
   */
  namespace AdjSeedmap3{


    struct sAdjSeedmap3 {
      sAdjRel3 *disp;  //displacement.
      sAdjRegion3 **uncertainty;
      sAdjRegion3 **object;
      sAdjRegion3 **obj_border;
      sAdjRegion3 **bkg_border;
      int nobjs;
    };


    sAdjSeedmap3 *Create(int nobjs);
    sAdjSeedmap3 *Create(sRegionCloud3 *rcloud);
    void          Destroy(sAdjSeedmap3 **asmap);

    
    void DrawObject(sScene8 *scn,
		    sAdjSeedmap3 *asmap,
		    Voxel u, int l, uchar val);
    void DrawObjBorder(sScene8 *scn,
		       sAdjSeedmap3 *asmap,
		       Voxel u, int l, uchar val);
    void DrawBkgBorder(sScene8 *scn,
		       sAdjSeedmap3 *asmap,
		       Voxel u, int l, uchar val);
    void DrawUncertainty(sScene8 *scn,
			 sAdjSeedmap3 *asmap,
			 Voxel u, int l, uchar val);
    
    void CopyUncertainty(sScene8 *dest, 
			 sScene8 *src,
			 sAdjSeedmap3 *asmap,
			 Voxel u, int l);

    void AddUncertainty(sScene8 *dest, 
			sScene8 *src,
			sAdjSeedmap3 *asmap,
			Voxel u, int l);

    void CloudArcWeight(sScene16 *arcw,
			sScene16 *wobj,
			sGradient3 *grad,
			Voxel u,
			sBorderCloud3 *bcloud,
			sAdjSeedmap3 *asmap,
			int l, float w);
    

  } //end AdjSeedmap3 namespace

  typedef AdjSeedmap3::sAdjSeedmap3 sAdjSeedmap3;
  
} //end gft namespace

#endif


