#ifndef _GFT_IMAGE32_H_
#define _GFT_IMAGE32_H_

#include "gft_common.h"

namespace gft{

  typedef enum {none, linear, cubic} InterpolationType;

  namespace Image32{

    /**
     * It supports both linear and two-dimensional access 
     * (i.e., img->data[p] or img->array[y][x] for a pixel
     * (x,y) at address p=x+y*xsize).
     */
    struct sImage32 {
      int *data;
      int **array;
      int nrows; /* numero de linhas (altura) */
      int ncols; /* numero de colunas (largura) */
      int n;     /* numero de pixels */
      float dx;
      float dy;
    };

    /**
     * \brief A constructor.
     */
    sImage32 *Create(int ncols,int nrows);
    sImage32 *Create(sImage32 *img);
    
    /**
     * \brief A destructor.
     */
    void    Destroy(sImage32 **img);

    /**
     * \brief A copy constructor.
     */
    sImage32 *Clone(sImage32 *img);
    sImage32 *Clone(sImage32 *img, Pixel l, Pixel h);

    void     Copy(sImage32 *img,  sImage32 *sub, Pixel l);
    void     Copy(sImage32 *dest, sImage32 *src);

    sImage32 *Add( sImage32 *img1, sImage32 *img2);
    sImage32 *Add( sImage32 *img,  int value);
    sImage32 *Mult(sImage32 *img1, sImage32 *img2);
    
    sImage32 *Read(char *filename);
    void      Write(sImage32 *img, char *filename);

    sImage32 *ConvertToNbits(sImage32 *img, int N);
    sImage32 *ConvertToNbits(sImage32 *img, int N, bool Imin);
    sImage32 *Complement(sImage32 *img);
    
    int     GetMinVal(sImage32 *img);
    int     GetMaxVal(sImage32 *img);
    int     GetMaxVal(sImage32 *img, int ignoredvalue);

    int     GetFreqVal(sImage32 *img, int val);
    
    void    Set(sImage32 *img, int value);

    bool    IsValidPixel(sImage32 *img, int x, int y);

    sImage32 *Threshold(sImage32 *img, int L, int H);
    
    void    DrawRectangle(sImage32 *img, 
			  int x1, int y1, 
			  int x2, int y2, int val);
    void    DrawLineDDA(sImage32 *img, 
			int x1, int y1, 
			int xn, int yn, int val);

    void    DrawCircle(sImage32 *img,
		       int x1, int y1,
		       float r,
		       int val);

    //------------------------------------

    sImage32 *AddFrame(sImage32 *img, int sz, int value);
    sImage32 *RemFrame(sImage32 *fimg, int sz);

    //------------------------------------

    sImage32 *Scale(sImage32 *img, float Sx, float Sy,
		    InterpolationType interpolation);

    sImage32 *MBB(sImage32 *img);
    void      MBB(sImage32 *img, Pixel *l, Pixel *h);
 
   
  } //end Image32 namespace

  typedef Image32::sImage32 sImage32;

} //end gft namespace



#endif

