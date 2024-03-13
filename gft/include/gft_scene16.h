
#ifndef _GFT_SCENE16_H_
#define _GFT_SCENE16_H_

#include "gft_common.h"

namespace gft{
  namespace Scene16{

  /**
   * It supports both linear and three-dimensional access 
   * (i.e., scn->data[p] or scn->array[z][y][x] for a voxel
   * (x,y,z) at address p=x+y*xsize+z*xsize*ysize).
   */
  struct sScene16 {
    ushort *data;
    ushort ***array;
    int xsize,ysize,zsize;
    float dx,dy,dz;
    ushort maxval;
    int n;
  };


  /**
   * \brief A constructor.
   */
  sScene16 *Create(int xsize,int ysize,int zsize);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  sScene16 *Create(sScene16 *scn);
  /**
   * \brief A destructor.
   */
  void     Destroy(sScene16 **scn);
  void     Copy(sScene16 *dest, sScene16 *src);
  void     Copy(sScene16 *dest, sScene16 *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  sScene16 *Clone(sScene16 *scn);
  sScene16 *SubScene(sScene16 *scn, Voxel l, Voxel h);
  sScene16 *SubScene(sScene16 *scn,
		     int xl, int yl, int zl,
		     int xh, int yh, int zh);
  void     Fill(sScene16 *scn, ushort value);

  void     Write(sScene16 *scn, char *filename);

  inline ushort GetValue(sScene16 *scn, Voxel v);
  inline ushort GetValue(sScene16 *scn, int p);
  inline ushort GetValue(sScene16 *scn, int x, int y, int z);
  ushort GetValue_nn(sScene16 *scn, float x, float y, float z);

  inline int GetAddressX(sScene16 *scn, int p);
  inline int GetAddressY(sScene16 *scn, int p);
  inline int GetAddressZ(sScene16 *scn, int p);
  inline int GetVoxelAddress(sScene16 *scn, Voxel v);
  inline int GetVoxelAddress(sScene16 *scn, int x, int y, int z);

  inline bool IsValidVoxel(sScene16 *scn, int x, int y, int z);
  inline bool IsValidVoxel(sScene16 *scn, Voxel v);
  
  ushort GetMaximumValue(sScene16 *scn);
  ushort GetMinimumValue(sScene16 *scn);

  sScene16 *MBB(sScene16 *scn);
  void      MBB(sScene16 *scn, Voxel *l, Voxel *h);

  sScene16 *AddFrame(sScene16 *scn,  int sz, ushort value);
  sScene16 *RemFrame(sScene16 *fscn, int sz);


  //---------inline definitions------------------
  inline bool  IsValidVoxel(sScene16 *scn, int x, int y, int z){
    if((x >= 0)&&(x < scn->xsize)&&
       (y >= 0)&&(y < scn->ysize)&&
       (z >= 0)&&(z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline bool IsValidVoxel(sScene16 *scn, Voxel v){
    if((v.c.x >= 0)&&(v.c.x < scn->xsize)&&
       (v.c.y >= 0)&&(v.c.y < scn->ysize)&&
       (v.c.z >= 0)&&(v.c.z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline ushort  GetValue(sScene16 *scn, Voxel v){
    return (scn->array[v.c.z][v.c.y][v.c.x]);
  }
  inline ushort  GetValue(sScene16 *scn, int p){
    return (scn->data[p]);
  }
  inline ushort  GetValue(sScene16 *scn, int x, int y, int z){
    return (scn->array[z][y][x]);
  }

  inline int GetAddressX(sScene16 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))%scn->xsize);
  }
  inline int GetAddressY(sScene16 *scn, int p){
    return ((p%(scn->xsize*scn->ysize))/scn->xsize);
  }
  inline int GetAddressZ(sScene16 *scn, int p){
    return (p/(scn->xsize*scn->ysize));
  }
  inline int GetVoxelAddress(sScene16 *scn, Voxel v){
    return (v.c.x + v.c.y*scn->xsize + 
	    v.c.z*scn->xsize*scn->ysize);
  }
  inline int GetVoxelAddress(sScene16 *scn, int x, int y, int z){
    return (x + y*scn->xsize + z*scn->xsize*scn->ysize);
  }

  
  } //end Scene16 namespace

  typedef Scene16::sScene16 sScene16;  

} //end gft namespace


#include "gft_scene32.h"
#include "gft_scene8.h"

namespace gft{
  namespace Scene16{

    sScene32 *ConvertTo32(sScene16 *scn);
    sScene8  *ConvertTo8(sScene16 *scn);

  } //end Scene16 namespace
} //end gft namespace


#endif
