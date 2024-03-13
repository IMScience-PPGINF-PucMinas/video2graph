
#ifndef _GFT_MATRIX_H_
#define _GFT_MATRIX_H_

#include "gft_common.h"
#include "gft_image32.h"

namespace gft{
  namespace Matrix{

    /**
     * It supports both linear and two-dimensional access 
     * (i.e., M->data[p] or M->array[i][j] for an entry
     * (i,j) at address p=j+i*ncols).
     */
    struct sMatrix {
      float *data;
      float **array;
      int ncols,nrows;
    };


    sMatrix *Create(int ncols,int nrows);
    void     Destroy(sMatrix **mat);
    sMatrix *Clone(sMatrix *mat);
    void     Copy(sMatrix *dest, 
		  sMatrix *src);
    
    sMatrix *Invert(sMatrix *A);
    sMatrix *Transpose(sMatrix *A);
    
    sMatrix *Mult(sMatrix *A, 
		  sMatrix *B);
    sMatrix *MultByScalar(sMatrix *A, float k);

    sMatrix *Sub(sMatrix *A, 
		 sMatrix *B);
    sMatrix *Add(sMatrix *A, 
		 sMatrix *B);

    float   GetTrace(sMatrix *M);
    
    void    Print(sMatrix *M);
    void    PrintDimension(sMatrix *M);
    
    float   ComputeDistanceL2(sMatrix *Y, 
			      sMatrix *X);
    
    void    Fill(sMatrix *M, float value);
    void    ChangeValue(sMatrix *M, 
			float old_value,
			float new_value);
    
    bool    IsValidEntry(sMatrix *M,
			 int i, int j);
    
    sImage32 *Convert2Image(sMatrix *M);
    
    sMatrix *Read(char *filename);
    void     Write(sMatrix *M,
		   char *filename);
    
    float   GetMinimumValue(sMatrix *M);
    float   GetMaximumValue(sMatrix *M);

    sImage32 *Threshold(sMatrix *M,
			float lower, 
			float higher);
    
    /**
     * @param axis an option (0->x / 1->y / 2->z).
     */
    sMatrix* RotationMatrix3(int axis, 
			     float th); 
    
    sMatrix* TranslationMatrix3(float dx, float dy, float dz);
    
    sMatrix* TransformVoxel(sMatrix *m, Voxel v);

  } //end Matrix namespace

  typedef Matrix::sMatrix sMatrix;

} //end gft namespace


#endif


