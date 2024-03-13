
#ifndef _GFT_SCENE64F_H_
#define _GFT_SCENE64F_H_

#include "gft_common.h"

namespace gft{
  namespace Scene64f{

  /**
   * It supports both linear and three-dimensional access 
   * (i.e., scn->data[p] or scn->array[z][y][x] for a voxel
   * (x,y,z) at address p=x+y*xsize+z*xsize*ysize).
   */
  struct sScene64f {
    double *data;
    double ***array;
    int xsize,ysize,zsize;
    float dx,dy,dz;
    double maxval;
    int n;
  };

  /**
   * \brief A constructor.
   */
  sScene64f *Create(int xsize,int ysize,int zsize);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  sScene64f *Create(sScene64f *scn);
  /**
   * \brief A destructor.
   */
  void     Destroy(sScene64f **scn);
  void     Copy(sScene64f *dest, sScene64f *src);
  void     Copy(sScene64f *dest, sScene64f *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  sScene64f *Clone(sScene64f *scn);
  sScene64f *SubScene(sScene64f *scn, Voxel l, Voxel h);
  sScene64f *SubScene(sScene64f *scn,
		      int xl, int yl, int zl,
		      int xh, int yh, int zh);
  void     Fill(sScene64f *scn, double value);

  //Scene64f *Read(char *filename);
  void     Write(sScene64f *scn, char *filename);

  inline double   GetValue(sScene64f *scn, Voxel v);
  inline double   GetValue(sScene64f *scn, int p);
  inline double   GetValue(sScene64f *scn, int x, int y, int z);
  double GetValue_trilinear(sScene64f *scn, float x, float y, float z);
  double   GetValue_nn(sScene64f *scn, float x, float y, float z);

  inline int GetAddressX(sScene64f *scn, int p);
  inline int GetAddressY(sScene64f *scn, int p);
  inline int GetAddressZ(sScene64f *scn, int p);
  inline int GetVoxelAddress(sScene64f *scn, Voxel v);
  inline int GetVoxelAddress(sScene64f *scn, int x, int y, int z);

  inline bool  IsValidVoxel(sScene64f *scn, int x, int y, int z);
  inline bool  IsValidVoxel(sScene64f *scn, Voxel v);
  
  double GetMaximumValue(sScene64f *scn);
  double GetMinimumValue(sScene64f *scn);

  sScene64f *MBB(sScene64f *scn);
  void       MBB(sScene64f *scn, Voxel *l, Voxel *h);

  sScene64f *AddFrame(sScene64f *scn,  int sz, double value);
  sScene64f *RemFrame(sScene64f *fscn, int sz);


  //---------inline definitions------------------
  inline bool  IsValidVoxel(sScene64f *scn, int x, int y, int z){
    if((x >= 0)&&(x < scn->xsize)&&
       (y >= 0)&&(y < scn->ysize)&&
       (z >= 0)&&(z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline bool  IsValidVoxel(sScene64f *scn, Voxel v){
    if((v.c.x >= 0)&&(v.c.x < scn->xsize)&&
       (v.c.y >= 0)&&(v.c.y < scn->ysize)&&
       (v.c.z >= 0)&&(v.c.z < scn->zsize))
      return(true);
    else
      return(false);
  }

  inline double   GetValue(sScene64f *scn, Voxel v){
    return (scn->array[v.c.z][v.c.y][v.c.x]);
  }
  inline double   GetValue(sScene64f *scn, int p){
    return (scn->data[p]);
  }
  inline double   GetValue(sScene64f *scn, int x, int y, int z){
    return (scn->array[z][y][x]);
  }

  inline int GetAddressX(sScene64f *scn, int p){
    return ((p%(scn->xsize*scn->ysize))%scn->xsize);
  }
  inline int GetAddressY(sScene64f *scn, int p){
    return ((p%(scn->xsize*scn->ysize))/scn->xsize);
  }
  inline int GetAddressZ(sScene64f *scn, int p){
    return (p/(scn->xsize*scn->ysize));
  }
  inline int GetVoxelAddress(sScene64f *scn, Voxel v){
    return (v.c.x + v.c.y*scn->xsize + 
	    v.c.z*scn->xsize*scn->ysize);
  }
  inline int GetVoxelAddress(sScene64f *scn, int x, int y, int z){
    return (x + y*scn->xsize + z*scn->xsize*scn->ysize);
  }

  
  } //end Scene64f namespace

  typedef Scene64f::sScene64f sScene64f;

} //end gft namespace


#include "gft_scene16.h"
#include "gft_scene8.h"

namespace gft{
  namespace Scene64f{

    sScene16 *ConvertTo16(sScene64f *scn);
    sScene8  *ConvertTo8(sScene64f *scn);

  } //end Scene64f namespace
} //end gft namespace


#endif
