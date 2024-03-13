#ifndef _GFT_IMAGE64_H_
#define _GFT_IMAGE64_H_

#include "gft_common.h"
#include "gft_image32.h"

namespace gft{
  namespace Image64{

    struct sImage64 {
      long long *data;
      long long **array; /* long long signed integer type.
			    At least 64 bits in size.
			    Specified since the C99 version
			    of the standard. */
      int nrows; /* numero de linhas (altura) */
      int ncols; /* numero de colunas (largura) */
      int n;     /* numero de pixels */
    };


    sImage64 *Create(int ncols, int nrows);
    void      Destroy(sImage64 **img);
    
    sImage64 *ConvertToImage64(sImage32 *img);
    
    sImage64 *Read(char *filename);
    void      Write(sImage64 *img, char *filename);
    sImage64 *ComputeIntegralImage(sImage64 *img);
    
    /* recebe uma imagem integral e calcula
       a soma dentro do retangulo com x em [x1, x2] e y em [y1, y2] */
    long long ComputeSum(sImage64 *iimg,
			 int x1, int y1,
			 int x2, int y2);
    
    /* Calcula a variancia do retangulo
       com x no intervalo fechado [x1, x2], e y em [y1, y2] */
    double ComputeVariance(sImage64 *img,
			   int x1, int y1,
			   int x2, int y2);
    
    
    /* Calcula a variancia do retangulo
       com x no intervalo fechado [x1, x2], e y em [y1, y2] */
    double ComputeVariance(sImage64 *iimg,
			   sImage64 *iimg2,
			   int x1, int y1,
			   int x2, int y2);
    
    sImage64 *Squared(sImage64 *img);
    
    
    /* Tamanho do lado do quadrado.*/
    int GetPosWindowMaxSum(sImage32 *img, int tam_x, int tam_y);

  } /*end Image64 namespace*/

  typedef Image64::sImage64 sImage64;

} /*end gft namespace*/

#endif

