#ifndef _GFT_FUZZYCLOUD3_H_
#define _GFT_FUZZYCLOUD3_H_

#include "gft_common.h"
#include "gft_adjrel3.h"
#include "gft_scene32.h"
#include "gft_filelist.h"
#include "gft_adjrel3f.h"
#include "gft_gradient3.h"

/*#include "cloud3.h"*/
/*#include "gradient.h"*/
/*#include "scene_addons.h"*/

#define MAX_PROB 100000

namespace gft{
  namespace RegionCloud3{

    struct sRegionCloud3 {
      sAdjRel3f *disp;  //displacement.
      sScene32  **prob;
      int nobjs;
      int nimages;
      //private:
      sAdjRel3f *fdisp;  //float displacement.
    };

    sRegionCloud3 *ByLabelList(sFileList *L);
    void           Destroy(sRegionCloud3 **rcloud);
    void           GetVoxelSize(sRegionCloud3 *rcloud,
				float *dx, 
				float *dy, 
				float *dz);
    sRegionCloud3 *Subsampling(sRegionCloud3 *rcloud);
    sRegionCloud3 *LinearInterp(sRegionCloud3 *rcloud,
				float dx,float dy,float dz);
    sRegionCloud3 *GaussianBlur(sRegionCloud3 *rcloud);
    sRegionCloud3 *ChangeOrientationToLPS(sRegionCloud3 *rcloud,
					  char *ori);
    
    sRegionCloud3 *Read(char *filename);
    void           Write(sRegionCloud3 *rcloud, 
			 char *filename);
    
    void RemoveElem(sRegionCloud3 *rcloud,
		    sScene32 *label);

  } //end RegionCloud3 namespace

  typedef RegionCloud3::sRegionCloud3 sRegionCloud3;

} //end gft namespace
   
    //---------------------------------------------

namespace gft{
  namespace BorderCloud3{
    
    struct sBorderCloud3 {
      sAdjRel3 *disp;  //displacement.
      sGradient3 **prob;
      int nobjs;
      //private:
      sAdjRel3f *fdisp;  //float displacement.
    };

    sBorderCloud3 *ByRegionCloud(sRegionCloud3 *rcloud);
    
    void          Destroy(sBorderCloud3 **bcloud);
    void          GetVoxelSize(sBorderCloud3 *bcloud,
			       float *dx, 
			       float *dy, 
			       float *dz);
    sBorderCloud3 *Subsampling(sBorderCloud3 *bcloud);
    sBorderCloud3 *LinearInterp(sBorderCloud3 *bcloud,
			       float dx,float dy,float dz);
    sBorderCloud3 *ChangeOrientationToLPS(sBorderCloud3 *bcloud,
					 char *ori);
    void          Normalize(sBorderCloud3 *bcloud,
			    int omin,int omax,
			    int nmin,int nmax);
    
    sBorderCloud3 *Read(char *filename);
    void           Write(sBorderCloud3 *bcloud,
			 char *filename);
    
  } //end BorderCloud3 namespace

  typedef BorderCloud3::sBorderCloud3 sBorderCloud3;

} //end gft namespace
   
#endif


