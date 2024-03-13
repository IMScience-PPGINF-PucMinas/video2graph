
#ifndef _GFT_SUPERPIXELS_H_
#define _GFT_SUPERPIXELS_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_cimage.h"
#include "gft_cimage32f.h"
#include "gft_color.h"
#include "gft_adjrel.h"
#include "gft_filtering.h"
#include "gft_queue.h"
#include "gft_heap.h"
#include "gft_scene32.h"
#include "gft_adjrel3.h"


namespace gft{
  namespace Superpixels{

    struct sFeature{
      int j;
      int i;
      float l;
      float a;
      float b;
    };

    struct sSPixelSeed{
      sFeature actual;
      sFeature last_computed;
      bool recompute;
      bool updateseeds;
      int n;
    };

    //--------------------------

    float ColorSquaredDistance(sCImage32f *cimg_lab, 
			       const sSPixelSeed& seed, int q);
    
    float LastColorSquaredDistance(const sSPixelSeed& seed);
    float LastSpatialSquaredDistance(const sSPixelSeed& seed);
    
    void RunSPixelsByIFT(sCImage32f *cimg_lab,
			 sImage32 *label,
			 sHeap *Q,
			 sAdjRel *A,
			 sAdjPxl *P,
			 sImage32 *pred,
			 float *cost,
			 float alpha,
			 //float beta,
			 sSPixelSeed *seeds, int nseeds);

    void RunSPixelsByDIFT(sCImage32f *cimg_lab,
			  sImage32 *label,
			  sHeap *Q,
			  sAdjRel *A,
			  sAdjPxl *P,
			  sImage32 *pred,
			  float *cost,
			  float alpha,
			  //float beta,
			  sSPixelSeed *seeds, int nseeds);
    
    sImage32 *IFT_SLIC(sCImage *cimg,
		       int k, float alpha, //float beta,
		       float err_dc, float err_ds, int itMax);
    sImage32 *IFT_SLIC(sImage32 *img,
		       int k, float alpha, //float beta,
		       float err_dc, float err_ds, int itMax);
    //--------------------------
    float WeightedDistanceMeasure(sCImage32f *cimg_lab, 
				  sSPixelSeed seed, int i, int j,
				  float m, float S);

    sImage32 *mySLIC(sCImage *cimg, 
		     int k, float m);
    sImage32 *mySLIC(sImage32 *img, 
		     int k, float m);
    //--------------------------
    sSPixelSeed *GetSPixelsSeeds(sCImage32f *cimg_lab,
				 sImage32 *grad,
				 int k,
				 int *nseeds);
    sSPixelSeed *GetSPixelsSeedsHexagon(sCImage32f *cimg_lab,
					sImage32 *grad,
					int k,
					int *nseeds);
    
    void Move2LowerGradient(sImage32 *grad, 
			    int *j, int *i);

    void UpdateSPixelsSeeds(sCImage32f *cimg_lab,
			    sImage32 *label,
			    sImage32 *pred,
			    sAdjPxl *P,
			    sSPixelSeed *seeds, int nseeds,
			    int inside);

    int GetNumberOfSuperPixels(sImage32 *label);


    void Postprocessing_CC(sImage32 *label,
			   sSPixelSeed *seeds, int nseeds);
    void Postprocessing(sImage32 *label,
			sSPixelSeed *seeds, int nseeds);


    //--------------------3D--------------------------------
    struct sFeature3D{
      int j;
      int i;
      int k;
      float l;
    };

    struct sSVoxelSeed{
      sFeature3D actual;
      sFeature3D last_computed;
      bool recompute;
      bool updateseeds;
      int n;
    };


    float FeatureSquaredDistance(sScene32 *scn, 
				 const sSVoxelSeed& seed, int q);
    
    float LastFeatureSquaredDistance(const sSVoxelSeed& seed);
    float LastSpatialSquaredDistance(const sSVoxelSeed& seed);
    
    void RunSVoxelsByIFT(sScene32 *scn,
			 sScene32 *label,
			 sHeap *Q,
			 sAdjRel3 *A,
			 sAdjVxl *V,
			 sScene32 *pred,
			 float *cost,
			 float alpha,
			 //float beta,
			 sSVoxelSeed *seeds, int nseeds);

    void RunSVoxelsByDIFT(sScene32 *scn,
			  sScene32 *label,
			  sHeap *Q,
			  sAdjRel3 *A,
			  sAdjVxl *V,
			  sScene32 *pred,
			  float *cost,
			  float alpha,
			  //float beta,
			  sSVoxelSeed *seeds, int nseeds);
    
    sScene32 *IFT_SLIC(sScene32 *scn,
		       int k, float alpha, //float beta,
		       float err_df, float err_ds, int itMax);

    sSVoxelSeed *GetSVoxelsSeeds(sScene32 *scn,
				 int k,
				 int *nseeds);
    
    void UpdateSVoxelsSeeds(sScene32 *scn,
			    sScene32 *label,
			    sScene32 *pred,
			    sAdjVxl *V,
			    sSVoxelSeed *seeds, int nseeds,
			    int inside);
    
  } //end Superpixels namespace
} //end gft namespace


#endif

