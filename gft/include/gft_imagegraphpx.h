/*
IFT:
custo favoravel   0 
custo normal      > 0, < INT_MAX
aresta impossivel INT_MAX

GC:
custo favoravel   INT_MAX = CAP_MAX
custo normal      > 0, < INT_MAX
aresta impossivel 0
*/

#ifndef _GFT_IMAGEGRAPHPX_H_
#define _GFT_IMAGEGRAPHPX_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_cimage.h"
#include "gft_adjrel.h"
#include "gft_filtering.h"


namespace gft{
  namespace ImageGraphPx{

#define DISSIMILARITY 0
#define CAPACITY      1
    
    // Image Sparse Graph:
    // n-links connect pairs of neighboring pixels.
    struct sImageGraphPx {
      int type;
      int Wmax;
      int ncols, nrows;
      sAdjRel *A;
      /*private: References to external images*/
      sImage32 *W;
      sImage32 *img;
      sCImage *cimg;
    };
    
    sImageGraphPx *ByAbsDiff(sImage32 *img, float r);
    sImageGraphPx *ByAbsDiff(sCImage *cimg,  float r);
    sImageGraphPx *ByWeightImage(sImage32 *W, float r);
    
    //------ ImageGraphPx Functions ------------
    
    void           Destroy(sImageGraphPx **g);
    sImageGraphPx  *Create(int ncols, int nrows, 
			   sAdjRel *A);
    sImageGraphPx  *Clone(sImageGraphPx *g);

    /*
    void          ChangeType(sImageGraphPx *sg, 
                             int type);
    
    void          Pow(sImageGraphPx *sg, int power, int max);
    
    void Orient2Digraph(sImageGraphPx *sg, 
			sImage32 *img,
			float per);
    
    void Orient2DigraphInner(sImageGraphPx *sg, 
			     sImage32 *P_sum);
    void Orient2DigraphOuter(sImageGraphPx *sg, 
			     sImage32 *P_sum);
    */
    
  } //end ImageGraphPx namespace

  typedef ImageGraphPx::sImageGraphPx sImageGraphPx;

  
} //end gft namespace


#endif


