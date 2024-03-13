
#ifndef _GFT_SCENE_H_
#define _GFT_SCENE_H_

#include "gft_common.h"
#include "gft_scene8.h"
#include "gft_scene16.h"
#include "gft_scene32.h"
#include "gft_scene32f.h"
#include "gft_scene64f.h"

typedef enum {integer, reals} gft_SceneType;

namespace gft{
  namespace Scene{

  struct sScene {
    uchar nbits;
    gft_SceneType type;
    union{
      sScene8   *scn8;
      sScene16  *scn16;
      sScene32  *scn32;
      sScene32f *scn32f;
      sScene64f *scn64f;
    } ptr;
  };


  /**
   * \brief A constructor.
   */
  sScene *Create(int xsize,int ysize,int zsize, int nbits, gft_SceneType type);
  /**
   * \brief A constructor for integer data.
   */
  sScene *Create(int xsize,int ysize,int zsize, int nbits);
  /**
   * \brief A constructor taking a reference scene as template.
   */
  sScene *Create(sScene *scn);
  /**
   * \brief A destructor.
   */
  void   Destroy(sScene **scn);
  void   Copy(sScene *dest, sScene *src);
  void   Copy(sScene *dest, sScene *src, Voxel v);
  /**
   * \brief A copy constructor.
   */
  sScene *Clone(sScene *scn);
  sScene *SubScene(sScene *scn, Voxel l, Voxel h);
  sScene *SubScene(sScene *scn,
		   int xl, int yl, int zl,
		   int xh, int yh, int zh);
  void   Fill(sScene *scn, int value);

  sScene *Read(char *filename);
  void    Write(sScene *scn, char *filename);

  void   SetValue(sScene *scn, int p, int value);

  int    GetValue(sScene *scn, Voxel v);
  int    GetValue(sScene *scn, int p);
  int    GetValue(sScene *scn, int x, int y, int z);
  int    GetValue_nn(sScene *scn, float x, float y, float z);

  double GetValue64f(sScene *scn, Voxel v);
  double GetValue64f(sScene *scn, int p);
  double GetValue64f(sScene *scn, int x, int y, int z);
  double GetValue64f_nn(sScene *scn, float x, float y, float z);

  int    GetNumberOfVoxels(sScene *scn);

  int    GetAddressX(sScene *scn, int p);
  int    GetAddressY(sScene *scn, int p);
  int    GetAddressZ(sScene *scn, int p);
  int    GetVoxelAddress(sScene *scn, Voxel v);
  int    GetVoxelAddress(sScene *scn, int x, int y, int z);

  bool   IsValidVoxel(sScene *scn, int x, int y, int z);
  bool   IsValidVoxel(sScene *scn, Voxel v);
  
  int    GetMaximumValue(sScene *scn);
  int    GetMinimumValue(sScene *scn);

  sScene *MBB(sScene *scn);
  void    MBB(sScene *scn, Voxel *l, Voxel *h);

  sScene *AddFrame(sScene *scn,  int sz, int value);
  sScene *RemFrame(sScene *fscn, int sz);

  } //end Scene namespace

  typedef Scene::sScene sScene;
  
} //end gft namespace

#endif

