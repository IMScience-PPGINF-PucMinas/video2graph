#ifndef _GFT_ADJREL_H_
#define _GFT_ADJREL_H_

#include "gft_common.h"

namespace gft
{
  namespace AdjRel
  {

    struct sAdjRel
    {
      int *dx;
      int *dy;
      int *dz;
      int n;
    };

    struct sAdjPxl
    {
      int *dp;
      int n;
    };

    sAdjRel *Create(int n);
    void Destroy(sAdjRel **A);
    sAdjRel *Clone(sAdjRel *A);

    sAdjRel *Neighborhood_4(); /* 4-neighborhood */
    sAdjRel *Neighborhood_6(); /* 6-neighborhood */
    sAdjRel *Neighborhood_8(); /* 8-neighborhood */
    sAdjRel *Neighborhood_8_counterclockwise();
    sAdjRel *Neighborhood_8_clockwise();

    sAdjRel *Circular(float r);
    sAdjRel *Box(int ncols, int nrows);

    //-----------------------------------
    int GetFrameSize(sAdjRel *A);

    void Mult(sAdjRel *A, int val);

    sAdjRel *RightSide(sAdjRel *A);
    sAdjRel *LeftSide(sAdjRel *A);

    sAdjRel *RightSide8(sAdjRel *A);
    sAdjRel *LeftSide8(sAdjRel *A);

  } // namespace AdjRel

  typedef AdjRel::sAdjPxl sAdjPxl;
  typedef AdjRel::sAdjRel sAdjRel;

} // namespace gft

#include "gft_image32.h"

namespace gft
{

  namespace Image32
  {

    sImage32 *Render(sAdjRel *A);

    void DrawAdjRel(sImage32 *img,
                    sAdjRel *A,
                    int p, int val);

    sImage32 *GetBoundaries(sImage32 *img, sAdjRel *A);

  } // namespace Image32

} // namespace gft

#include "gft_cimage.h"

namespace gft
{
  namespace CImage
  {

    void DrawAdjRel(sCImage *cimg,
                    sAdjRel *A,
                    int p, int color);

  } // namespace CImage
} // namespace gft

namespace gft
{
  namespace AdjRel
  {

    sAdjPxl *AdjPixels(sAdjRel *A, int ncols);
    sAdjPxl *AdjPixels(sAdjRel *A, sImage32 *img);

    sAdjPxl *AdjPixels(sAdjRel *A, sCImage *cimg);

    void DestroyAdjPxl(sAdjPxl **N);

  } // namespace AdjRel
} // namespace gft

#endif
