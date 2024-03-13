
#ifndef _GFT_SCENE64_H_
#define _GFT_SCENE64_H_

#include "gft_common.h"

extern "C" {
#include "nifti1_io.h"
}

namespace gft{
  namespace Scene64{

  /**
   * It supports both linear and three-dimensional access 
   * (i.e., scn->data[p] or scn->array[z][y][x] for a voxel
   * (x,y,z) at address p=x+y*xsize+z*xsize*ysize).
   */
  struct sScene64 {
    long long *data;
    long long ***array;
    int xsize,ysize,zsize;
    float dx,dy,dz;
    long long maxval;
    int n;
    nifti_image *nii_hdr;
  };

  /**
   * \brief A constructor.
   */
  sScene64 *Create(int xsize,int ysize,int zsize);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  sScene64 *Create(sScene64 *scn);
  /**
   * \brief A destructor.
   */
  void     Destroy(sScene64 **scn);
  void     Copy(sScene64 *dest, sScene64 *src);
  void     Copy(sScene64 *dest, sScene64 *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  sScene64 *Clone(sScene64 *scn);
  sScene64 *SubScene(sScene64 *scn, Voxel l, Voxel h);
  sScene64 *SubScene(sScene64 *scn,
		     int xl, int yl, int zl,
		     int xh, int yh, int zh);
  void     Fill(sScene64 *scn, long long value);


  inline int GetAddressX(sScene64 *scn, int p);
  inline int GetAddressY(sScene64 *scn, int p);
  inline int GetAddressZ(sScene64 *scn, int p);
  inline int GetVoxelAddress(sScene64 *scn, Voxel v);
  inline int GetVoxelAddress(sScene64 *scn, int x, int y, int z);

  inline bool  IsValidVoxel(sScene64 *scn, int x, int y, int z);
  inline bool  IsValidVoxel(sScene64 *scn, Voxel v);
  
  long long    GetMaximumValue(sScene64 *scn);
  long long    GetMinimumValue(sScene64 *scn);


  //---------inline definitions------------------
  inline bool  IsValidVoxel(sScene64 *scn, int x, int y, int z){
    if((x >= 0)&&(x < scn->xsize)&&
       (y >= 0)&&(y < scn->ysize)&&
       (z >= 0)&&(z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline bool  IsValidVoxel(sScene64 *scn, Voxel v){
    if((v.c.x >= 0)&&(v.c.x < scn->xsize)&&
       (v.c.y >= 0)&&(v.c.y < scn->ysize)&&
       (v.c.z >= 0)&&(v.c.z < scn->zsize))
      return(true);
    else
      return(false);
  }

  
  inline int GetAddressX(sScene64 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))%scn->xsize);
  }
  inline int GetAddressY(sScene64 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))/scn->xsize);
  }
  inline int GetAddressZ(sScene64 *scn, int p){
    return (p/(scn->xsize*scn->ysize));
  }
  inline int GetVoxelAddress(sScene64 *scn, Voxel v){
    return (v.c.x + v.c.y*scn->xsize + 
	    v.c.z*scn->xsize*scn->ysize);
  }
  inline int GetVoxelAddress(sScene64 *scn, int x, int y, int z){
    return (x + y*scn->xsize + z*scn->xsize*scn->ysize);
  }


  //The dimensions of the window (i.e., xsize,ysize,zsize) 
  //should be given in millimeters.
  float ComputeWindowSum(sScene64 *Iscn,
			 float xsize, 
			 float ysize, 
			 float zsize,
			 Voxel u);

  //The dimensions of the window (i.e., xsize,ysize,zsize) 
  //should be given in millimeters.
  float ComputeWindowDensity(sScene64 *Iscn,
			     float xsize, 
			     float ysize, 
			     float zsize,
			     Voxel u);
  
  //The dimensions of the window (i.e., xsize,ysize,zsize) 
  //should be given in millimeters.
  Voxel  FindWindowOfMaximumSum(sScene64 *Iscn, 
				float xsize, 
				float ysize, 
				float zsize);
  
  
  } //end Scene64 namespace

  typedef Scene64::sScene64 sScene64;

} //end gft namespace


#endif
