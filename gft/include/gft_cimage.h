#ifndef _GFT_CIMAGE_H_
#define _GFT_CIMAGE_H_

#include "gft_image32.h"
#include "gft_color.h"

namespace gft{
  namespace CImage{

    struct sCImage {
      sImage32 *C[3];
    };

    sCImage *Create(int ncols, int nrows);
    sCImage *Create(sCImage *cimg);
    sCImage *Create(sImage32 *img);
    void    Destroy(sCImage **cimg);
    sCImage *Clone(sCImage *cimg);
    sCImage *Clone(sCImage *cimg, Pixel l, Pixel h);
    sCImage *Clone(sImage32 *img);
    sCImage *Clone(sImage32 *img, int Imin, int Imax);
    void    Copy(sCImage *cimg, sCImage *sub, Pixel l);

    sCImage *RandomColorize(sImage32 *img);
    
    sCImage *Read(char *filename);
    void    Write(sCImage *cimg, char *filename);
    void    Set(sCImage *cimg, int r, int g, int b);
    void    Set(sCImage *cimg, int color);
    
    sCImage *ColorizeLabel(sImage32 *label);

    sCImage *RGB2Lab(sCImage *cimg);
 
    sImage32 *Lightness(sCImage *cimg);

    /*The luminosity method works best overall and is the 
      default method used if you ask GIMP to change an image 
      from RGB to grayscale */
    sImage32 *Luminosity(sCImage *cimg);

    void MBB(sCImage *cimg, int bkgcolor, Pixel *l, Pixel *h);

    sCImage *AddFrame(sCImage *cimg, int sz, int r, int g, int b);
    sCImage *RemFrame(sCImage *cimg, int sz);

    void    DrawRectangle(sCImage *cimg, 
			  int x1, int y1, 
			  int x2, int y2, int color);
    void    DrawLineDDA(sCImage *cimg,
			int x1, int y1, 
			int xn, int yn, int color);

    void    DrawCircle(sCImage *cimg,
		       int x1, int y1,
		       float r,
		       int color);

    sCImage *Scale(sCImage *cimg, float Sx, float Sy,
		   gft::InterpolationType interpolation);
    
  } //end CImage namespace

  typedef CImage::sCImage sCImage;

} //end gft namespace

#endif

