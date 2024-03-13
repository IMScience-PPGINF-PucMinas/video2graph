
#ifndef _GFT_SCENE32F_H_
#define _GFT_SCENE32F_H_

#include "gft_common.h"

namespace gft{
  namespace Scene32f{

  /**
   * It supports both linear and three-dimensional access 
   * (i.e., scn->data[p] or scn->array[z][y][x] for a voxel
   * (x,y,z) at address p=x+y*xsize+z*xsize*ysize).
   */
  struct sScene32f {
    float *data;
    float ***array;
    int xsize,ysize,zsize;
    float dx,dy,dz;
    float maxval;
    int n;
  };

  /**
   * \brief A constructor.
   */
  sScene32f *Create(int xsize,int ysize,int zsize);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  sScene32f *Create(sScene32f *scn);
  /**
   * \brief A destructor.
   */
  void     Destroy(sScene32f **scn);
  void     Copy(sScene32f *dest, sScene32f *src);
  void     Copy(sScene32f *dest, sScene32f *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  sScene32f *Clone(sScene32f *scn);
  sScene32f *SubScene(sScene32f *scn, Voxel l, Voxel h);
  sScene32f *SubScene(sScene32f *scn,
		      int xl, int yl, int zl,
		      int xh, int yh, int zh);
  void     Fill(sScene32f *scn, float value);

  //Scene32f *Read(char *filename);
  void     Write(sScene32f *scn, char *filename);

  inline float   GetValue(sScene32f *scn, Voxel v);
  inline float   GetValue(sScene32f *scn, int p);
  inline float   GetValue(sScene32f *scn, int x, int y, int z);
  float GetValue_trilinear(sScene32f *scn, float x, float y, float z);
  float   GetValue_nn(sScene32f *scn, float x, float y, float z);

  inline int GetAddressX(sScene32f *scn, int p);
  inline int GetAddressY(sScene32f *scn, int p);
  inline int GetAddressZ(sScene32f *scn, int p);
  inline int GetVoxelAddress(sScene32f *scn, Voxel v);
  inline int GetVoxelAddress(sScene32f *scn, int x, int y, int z);

  inline bool  IsValidVoxel(sScene32f *scn, int x, int y, int z);
  inline bool  IsValidVoxel(sScene32f *scn, Voxel v);
  
  float      GetMaximumValue(sScene32f *scn);
  float      GetMinimumValue(sScene32f *scn);

  sScene32f *MBB(sScene32f *scn);
  void       MBB(sScene32f *scn, Voxel *l, Voxel *h);

  sScene32f *AddFrame(sScene32f *scn,  int sz, float value);
  sScene32f *RemFrame(sScene32f *fscn, int sz);


  //---------inline definitions------------------
  inline bool  IsValidVoxel(sScene32f *scn, int x, int y, int z){
    if((x >= 0)&&(x < scn->xsize)&&
       (y >= 0)&&(y < scn->ysize)&&
       (z >= 0)&&(z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline bool  IsValidVoxel(sScene32f *scn, Voxel v){
    if((v.c.x >= 0)&&(v.c.x < scn->xsize)&&
       (v.c.y >= 0)&&(v.c.y < scn->ysize)&&
       (v.c.z >= 0)&&(v.c.z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline float   GetValue(sScene32f *scn, Voxel v){
    return (scn->array[v.c.z][v.c.y][v.c.x]);
  }
  inline float   GetValue(sScene32f *scn, int p){
    return (scn->data[p]);
  }
  inline float   GetValue(sScene32f *scn, int x, int y, int z){
    return (scn->array[z][y][x]);
  }

  inline int GetAddressX(sScene32f *scn, int p){
    return ((p%(scn->xsize*scn->ysize))%scn->xsize);
  }
  inline int GetAddressY(sScene32f *scn, int p){
    return ((p%(scn->xsize*scn->ysize))/scn->xsize);
  }
  inline int GetAddressZ(sScene32f *scn, int p){
    return (p/(scn->xsize*scn->ysize));
  }
  inline int GetVoxelAddress(sScene32f *scn, Voxel v){
    return (v.c.x + v.c.y*scn->xsize + 
	    v.c.z*scn->xsize*scn->ysize);
  }
  inline int GetVoxelAddress(sScene32f *scn, int x, int y, int z){
    return (x + y*scn->xsize + z*scn->xsize*scn->ysize);
  }

  
  } //end Scene32f namespace

  typedef Scene32f::sScene32f sScene32f;

} //end gft namespace


#include "gft_scene16.h"
#include "gft_scene8.h"

namespace gft{
  namespace Scene32f{

    sScene16 *ConvertTo16(sScene32f *scn);
    sScene8  *ConvertTo8(sScene32f *scn);

  } //end Scene32f namespace
} //end gft namespace


#endif
