
#ifndef _GFT_SCENE8_H_
#define _GFT_SCENE8_H_

#include "gft_common.h"

namespace gft{
  namespace Scene8{

  /**
   * It supports both linear and three-dimensional access 
   * (i.e., scn->data[p] or scn->array[z][y][x] for a voxel
   * (x,y,z) at address p=x+y*xsize+z*xsize*ysize).
   */
  struct sScene8 {
    uchar *data;
    uchar ***array;
    int xsize,ysize,zsize;
    float dx,dy,dz;
    uchar maxval;
    int n;
  };

  /**
   * \brief A constructor.
   */
  sScene8 *Create(int xsize,int ysize,int zsize);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  sScene8 *Create(sScene8 *scn);
  /**
   * \brief A destructor.
   */
  void    Destroy(sScene8 **scn);
  void    Copy(sScene8 *dest, sScene8 *src);
  void    Copy(sScene8 *dest, sScene8 *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  sScene8 *Clone(sScene8 *scn);
  sScene8 *SubScene(sScene8 *scn, Voxel l, Voxel h);
  sScene8 *SubScene(sScene8 *scn,
		    int xl, int yl, int zl,
		    int xh, int yh, int zh);
  void    Fill(sScene8 *scn, uchar value);

  void    Write(sScene8 *scn, char *filename);

  inline uchar GetValue(sScene8 *scn, Voxel v);
  inline uchar GetValue(sScene8 *scn, int p);
  inline uchar GetValue(sScene8 *scn, int x, int y, int z);
  uchar  GetValue_nn(sScene8 *scn, float x, float y, float z);

  inline int GetAddressX(sScene8 *scn, int p);
  inline int GetAddressY(sScene8 *scn, int p);
  inline int GetAddressZ(sScene8 *scn, int p);
  inline int GetVoxelAddress(sScene8 *scn, Voxel v);
  inline int GetVoxelAddress(sScene8 *scn, int x, int y, int z);

  inline bool IsValidVoxel(sScene8 *scn, int x, int y, int z);
  inline bool IsValidVoxel(sScene8 *scn, Voxel v);
  
  uchar   GetMaximumValue(sScene8 *scn);
  uchar   GetMinimumValue(sScene8 *scn);

  sScene8 *MBB(sScene8 *scn);
  void     MBB(sScene8 *scn, Voxel *l, Voxel *h);

  Voxel   Centroid(sScene8 *scn);

  sScene8 *AddFrame(sScene8 *scn,  int sz, uchar value);
  sScene8 *RemFrame(sScene8 *fscn, int sz);


  //---------inline definitions------------------
  inline bool  IsValidVoxel(sScene8 *scn, int x, int y, int z){
    if((x >= 0)&&(x < scn->xsize)&&
       (y >= 0)&&(y < scn->ysize)&&
       (z >= 0)&&(z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline bool IsValidVoxel(sScene8 *scn, Voxel v){
    if((v.c.x >= 0)&&(v.c.x < scn->xsize)&&
       (v.c.y >= 0)&&(v.c.y < scn->ysize)&&
       (v.c.z >= 0)&&(v.c.z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline uchar GetValue(sScene8 *scn, Voxel v){
    return (scn->array[v.c.z][v.c.y][v.c.x]);
  }
  inline uchar GetValue(sScene8 *scn, int p){
    return (scn->data[p]);
  }
  inline uchar GetValue(sScene8 *scn, int x, int y, int z){
    return (scn->array[z][y][x]);
  }

  inline int GetAddressX(sScene8 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))%scn->xsize);
  }
  inline int GetAddressY(sScene8 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))/scn->xsize);
  }
  inline int GetAddressZ(sScene8 *scn, int p){
    return (p/(scn->xsize*scn->ysize));
  }
  inline int GetVoxelAddress(sScene8 *scn, Voxel v){
    return (v.c.x + v.c.y*scn->xsize + 
	    v.c.z*scn->xsize*scn->ysize);
  }
  inline int GetVoxelAddress(sScene8 *scn, int x, int y, int z){
    return (x + y*scn->xsize + z*scn->xsize*scn->ysize);
  }

  
  } //end Scene8 namespace

  typedef Scene8::sScene8 sScene8;

} //end gft namespace


#include "gft_scene32.h"
#include "gft_scene16.h"

namespace gft{
  namespace Scene8{

    sScene32 *ConvertTo32(sScene8 *scn);
    sScene16 *ConvertTo16(sScene8 *scn);

  } //end Scene8 namespace
} //end gft namespace


#endif
