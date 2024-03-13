
#ifndef _GFT_CURVE_H_
#define _GFT_CURVE_H_

#include "gft_common.h"

namespace gft{
  namespace Curve{

  struct sCurve {
    float *X;
    float *Y;
    int n;
  };

  sCurve *Create(int n);
  void    Destroy(sCurve **curve);
  sCurve *Clone(sCurve *curve);

  /**
   * Compatible with Gnuplot.
   */
  sCurve *Read(char *filename);

  /**
   * Compatible with Gnuplot.
   */
  void   Write(sCurve *curve, char *filename);

  sCurve *Normalize(sCurve *curve);
  void   Normalizeinplace(sCurve *curve);

  int    LowerPercentage(sCurve *curve,
			 float perc);
  int    HigherPercentage(sCurve *curve,
			  float perc);
  int    Median(sCurve *curve);
  int    Otsu(sCurve *curve);

  void   Sort(sCurve *curve);

  void   InvertXY(sCurve *curve);
  
  } //end Curve namespace

  typedef Curve::sCurve sCurve;

} //end gft namespace


#endif
