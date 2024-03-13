
#ifndef _GFT_RADIOMETRIC3_H_
#define _GFT_RADIOMETRIC3_H_

#include "gft_common.h"
#include "gft_scene.h"
#include "gft_curve.h"

namespace gft{


  namespace Scene32{

    sCurve *NormHistogram(sScene32 *scn);
    sCurve *Histogram(sScene32 *scn);
    sCurve *Histogram(sScene32 *scn, int binwidth);

    void LinearStretchinplace(sScene32 *scn, 
			      int omin,int omax,
			      int nmin,int nmax);
    sScene32 *LinearStretch(sScene32 *scn, 
			    int omin,int omax,
			    int nmin,int nmax);
  } //end Scene32 namespace


  namespace Scene16{
    sCurve *Histogram(sScene16 *scn, int binwidth);

    void LinearStretchinplace(sScene16 *scn, 
			      int omin,int omax,
			      int nmin,int nmax);
    sScene16 *LinearStretch(sScene16 *scn, 
			    int omin,int omax,
			    int nmin,int nmax);
  } //end Scene16 namespace


  namespace Scene8{
    sCurve *Histogram(sScene8 *scn, int binwidth);

    void LinearStretchinplace(sScene8 *scn, 
			      int omin,int omax,
			      int nmin,int nmax);
    sScene8 *LinearStretch(sScene8 *scn, 
			   int omin,int omax,
			   int nmin,int nmax);
  } //end Scene8 namespace


  namespace Scene{
    sCurve *Histogram(sScene *scn, int binwidth);

    void LinearStretchinplace(sScene *scn, 
			      int omin,int omax,
			      int nmin,int nmax);
    sScene *LinearStretch(sScene *scn, 
			  int omin,int omax,
			  int nmin,int nmax);
  } //end Scene namespace


} //end gft namespace

#endif

