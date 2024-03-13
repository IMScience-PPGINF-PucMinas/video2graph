
#ifndef _GFT_HIGHLIGHT_H_
#define _GFT_HIGHLIGHT_H_

#include "gft_common.h"
#include "gft_image32.h"
#include "gft_cimage.h"
#include "gft_adjrel.h"
#include "gft_color.h"

namespace gft{
  namespace Highlight{
    
    sImage32 *Wide(sImage32 *img, 
		   sImage32 *label, 
		   float radius, int value, bool fill);
    sCImage *CWide(sCImage *cimg, 
		   sImage32 *label, 
		   float radius, int color, bool fill);
    
    sCImage *CWideLabels(sCImage *cimg, 
			 sImage32 *label,
			 float radius, int *colormap, float fill,
			 bool thickborder, bool singlecolor);

    bool    HStripedTexture(int x, int y, int w, int h);
    bool    VStripedTexture(int x, int y, int w, int h);
    bool    BackslashTexture(int x, int y, int w, int h);
    bool    SlashTexture(int x, int y, int w, int h);
    bool    GridTexture(int x, int y, int w, int h);
    bool    RGridTexture(int x, int y, int w, int h);
    
    sImage32 *Texture(sImage32 *img,
		      sImage32 *label,
		      float radius, int value, bool fill,
		      bool (*texture)(int,int,int,int), 
		      int w, int h);
    sCImage *CTexture(sCImage *cimg,
		      sImage32 *label,
		      float radius, int color, bool fill,
		      bool (*texture)(int,int,int,int), int w, int h);

  } //end Highlight namespace
} //end gft namespace

#endif

