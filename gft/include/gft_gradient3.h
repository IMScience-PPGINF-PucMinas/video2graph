
#ifndef _GFT_GRADIENT3_H_
#define _GFT_GRADIENT3_H_

#include "gft_common.h"
#include "gft_scene32.h"
#include "gft_adjrel3.h"
#include "gft_scnmath.h"

namespace gft{
  namespace Gradient3{
    
    struct sGradient3 {
      sScene32 *Gx;
      sScene32 *Gy;
      sScene32 *Gz;
      
      //Available upon request:
      //--> Must call "ComputeMagnitude".
      sScene32 *mag;
    };


    sGradient3 *Create(int xsize,int ysize,int zsize);
    void        Destroy(sGradient3 **grad);
    sGradient3 *RemFrame(sGradient3 *fgrad, int sz);
    sGradient3 *LinearInterpCentr(sGradient3 *grad,
				  float dx, float dy, float dz);
    sGradient3 *ChangeOrientationToLPS(sGradient3 *grad,
				       char *ori);
    
    sGradient3 *Read(char *filename);
    void       Write(sGradient3 *grad, char *filename);
    
    sGradient3 *Spherical(sScene32 *scn, float r);
    
    void    ComputeMagnitude(sGradient3 *grad);
    int     MaximumMag(sGradient3 *grad);
    void    Normalize(sGradient3 *grad,
		      int omin,int omax,
		      int nmin,int nmax);
    
    //void    PowerEnhancement(sGradient3 *grad);

  } //end Gradient3 namespace

  typedef Gradient3::sGradient3 sGradient3;

} //end gft namespace


#endif

