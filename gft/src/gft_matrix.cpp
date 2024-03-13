
#include "gft_matrix.h"

namespace gft{
  namespace Matrix{

    sMatrix *Create(int ncols,int nrows){
      sMatrix *mat=NULL;
      float *aux;
      int i;
      
      mat = (sMatrix *) calloc(1,sizeof(sMatrix));
      if(mat == NULL)
	gft::Error((char *)MSG1,(char *)"Matrix::Create");
      
      aux = (float *)calloc(nrows*ncols, sizeof(float));
      mat->array = (float **) calloc(nrows, sizeof(float *));
      if(mat->array == NULL || aux == NULL)
	gft::Error((char *)MSG1,(char *)"Matrix::Create");

      mat->data = aux;
      mat->array[0] = aux;
      for(i=1; i<nrows; i++) 
	mat->array[i] = mat->array[i-1] + ncols;
    
      mat->ncols = ncols;
      mat->nrows = nrows;
      return(mat);
    }


    void   Destroy(sMatrix **mat){
      sMatrix *aux;
      aux = *mat;
      if(aux != NULL){
	if(aux->array != NULL){
	  if(*(aux->array) != NULL)
	    free(*(aux->array));
	  free(aux->array); 
	}
	free(aux);
	*mat = NULL;
      }
    }


    sMatrix *Clone(sMatrix *mat){
      sMatrix *matc;
      matc = Create(mat->ncols, mat->nrows);
      memcpy(matc->array[0], mat->array[0],
	     mat->ncols*mat->nrows*sizeof(float));
      return(matc);
    }
    

    void     Copy(sMatrix *dest, 
		  sMatrix *src){
      if(dest->ncols!=src->ncols ||
	 dest->nrows!=src->nrows)
	gft::Error((char *)"Incompatible matrices",
		   (char *)"Matrix::Copy");
      
      memcpy(dest->array[0], src->array[0], 
	     src->ncols*src->nrows*sizeof(float));
    }
    

    void   Print(sMatrix *M){
      int x, y;
  
      printf("\n");
      for(y=0; y<M->nrows; y++){
	for(x=0; x<M->ncols; x++)
	  printf("%3.5lf ", M->array[y][x]);
	printf("\n");
      }
      printf("\n");
    }


    void   PrintDimension(sMatrix *M){
      printf("\n%d x %d\n", M->nrows, M->ncols);
    }
    

    sMatrix *Invert(sMatrix *A){
      sMatrix *B=NULL;
      sMatrix *I=NULL;
      int i,j,k;
      float m;
  
      if(A->ncols!=A->nrows)
	gft::Error((char *)"Matrix dimension error",
		   (char *)"Matrix::Invert");
  
      I = Create(A->ncols, A->nrows);
      B = Clone(A);
      
      for(i=0; i<A->nrows; i++)
	I->array[i][i] = 1.0;

      for(k=0; k<A->nrows; k++){
	m = B->array[k][k];
	//if(m < 0.0000000000001)
	// gft::Error((char *)"Singular matrix",
	//            (char *)"Matrix::Invert"); THIS IS WRONG
	
	B->array[k][k] = 1.0;
	for(j=0; j<A->ncols; j++){
	  if(j!=k)
	    B->array[k][j] = B->array[k][j]/m;
	  I->array[k][j] = I->array[k][j]/m;
	}
	
	for(i=0; i<A->nrows; i++){
	  if(i!=k){
	    m = B->array[i][k]/B->array[k][k];
	    
	    B->array[i][k] = 0.0;
	    for(j=0; j<A->ncols; j++){
	      if(j!=k)
		B->array[i][j] = B->array[i][j] - m*B->array[k][j];
	      I->array[i][j] = I->array[i][j] - m*I->array[k][j];
	    }
	  }
	}
      }
      Destroy(&B);
      return I;
    }


    float        GetTrace(sMatrix *M){
      float sum;
      int i;
      
      if(M->ncols!=M->nrows)
	gft::Error((char *)"Matrix dimension error",
		   (char *)"Matrix::GetTrace");
      sum = 0.0;
      for(i=0; i<M->nrows; i++)
	sum += M->array[i][i];
      return sum;
    }
    
    
    sMatrix *Mult(sMatrix *A, 
		  sMatrix *B){
      sMatrix *M = NULL;
      int i,j,k;
      
      if(A->ncols!=B->nrows)
	gft::Error((char *)"Matrix dimension error",
		   (char *)"Matrix::Mult");
  
      M = Create(B->ncols, A->nrows);
      for(i=0; i<M->nrows; i++){
	for(j=0; j<M->ncols; j++){
	  M->array[i][j] = 0.0;
	  for (k=0; k<A->ncols; k++)
	    M->array[i][j] += A->array[i][k]*B->array[k][j];
	}
      }
      return(M);
    }
    
    
    sMatrix *MultByScalar(sMatrix *A, 
			  float k){
      sMatrix *M = NULL;
      int i,j;
      
      M = Create(A->ncols, A->nrows);
      for(i=0; i<M->nrows; i++){
	for(j=0; j<M->ncols; j++){
	  M->array[i][j] = k*A->array[i][j];
	}
      }
      return(M);
    }


    sMatrix *Sub(sMatrix *A, 
		 sMatrix *B){
      sMatrix *M = NULL;
      int i,j;
      
      if((A->ncols!=B->ncols)||(A->nrows!=B->nrows))
	gft::Error((char *)"Matrix dimension error",
		   (char *)"Matrix::Sub");
      M = Create(A->ncols, A->nrows);
      for(i=0; i<M->nrows; i++){
	for(j=0; j<M->ncols; j++){
	  M->array[i][j] = A->array[i][j] - B->array[i][j];
	}
      }
      return(M);
    }


    sMatrix *Add(sMatrix *A, 
		 sMatrix *B){
      sMatrix *M = NULL;
      int i,j;
      
      if((A->ncols!=B->ncols)||(A->nrows!=B->nrows))
	gft::Error((char *)"Matrix dimension error",
		   (char *)"Matrix::Add");
      
      M = Create(A->ncols, A->nrows);
      for(i=0; i<M->nrows; i++){
	for(j=0; j<M->ncols; j++){
	  M->array[i][j] = A->array[i][j] + B->array[i][j];
	}
      }
      return(M);
    }


    sMatrix *Transpose(sMatrix *A){
      sMatrix *M = NULL;
      int i,j;
      
      M = Create(A->nrows, A->ncols);
      for(i=0; i<M->nrows; i++){
	for(j=0; j<M->ncols; j++){
	  M->array[i][j] = A->array[j][i];
	}
      }
      return(M);
    }


    float ComputeDistanceL2(sMatrix *Y, 
			    sMatrix *X){
      sMatrix *A,*B,*R;
      float d;
  
      A = Sub(X, Y);
      B = Transpose(A);
      R = Mult(A, B);
      d = GetTrace(R);
      d = sqrtf(d);
      Destroy(&A);
      Destroy(&B);
      Destroy(&R);
      return (d);
    }

    
    sMatrix *Read(char *filename){
      sMatrix *M;
      char msg[512];
      int  ncols,nrows,size,n,p;
      double *daux=NULL;
      FILE *fp;

      fp = fopen(filename,"rb");
      if(fp == NULL){
	sprintf(msg,"Cannot open %s",filename);
	gft::Error((char *)msg,(char *)"Matrix::Read");
      }
      fread(&ncols, sizeof(int), 1, fp);
      fread(&nrows, sizeof(int), 1, fp);
      fread(&size,  sizeof(int), 1, fp);

      M = Create(ncols, nrows);
      n = ncols*nrows;
      if(size==sizeof(float))
	fread(M->array[0], sizeof(float), n, fp);
      else if(size==sizeof(double)){
	daux = AllocDoubleArray(n);
	fread(daux, sizeof(double), n, fp);
	for(p=0; p<n; p++)
	  M->array[0][p] = (float)daux[p];
	free(daux);
      }
      else
	gft::Error((char *)"Bad or corrupted file",
		   (char *)"Matrix::Read");
      fclose(fp);
      return M;
    }


    void    Write(sMatrix *M,
		  char *filename){
      char msg[512];
      int n,ncols,nrows,size;
      FILE *fp;
      
      fp = fopen(filename,"wb");
      if(fp == NULL){
	sprintf(msg,"Cannot open %s",filename);
	gft::Error((char *)msg,
		   (char *)"Matrix::Write");
      }
      size  = sizeof(float);
      ncols = M->ncols;
      nrows = M->nrows;
      n = ncols*nrows;
      fwrite(&ncols, sizeof(int),  1, fp);
      fwrite(&nrows, sizeof(int),  1, fp);
      fwrite(&size,  sizeof(int),  1, fp);
      fwrite(M->array[0], sizeof(float), n ,fp);
      fclose(fp);
    }


    bool   IsValidEntry(sMatrix *M,
			int i, int j){
      if((j >= 0)&&(j < M->ncols)&&
	 (i >= 0)&&(i < M->nrows))
	return(true);
      else
	return(false);
    }


    float   GetMaximumValue(sMatrix *M){
      float max;
      int p,n;
      
      n = M->ncols*M->nrows;
      max = M->array[0][0];
      for(p = 1; p < n; p++)
	if(M->array[0][p] > max)
	  max = M->array[0][p];
      return(max);
    }
    
    
    float    GetMinimumValue(sMatrix *M){
      float min;
      int p,n;
      n = M->ncols*M->nrows;
      min = M->array[0][0];
      for(p=1; p<n; p++)
	if(M->array[0][p] < min)
	  min = M->array[0][p];
      return(min);
    }


    void    Fill(sMatrix *M, float value){
      int p,n;
      n = M->ncols*M->nrows;
      for(p=0; p<n; p++)
	M->array[0][p] = value;
    }
    

    void    ChangeValue(sMatrix *M, 
			float old_value,
			float new_value){
      int p,n;
      
      n = M->ncols*M->nrows;
      for(p = 0; p < n; p++)
	if(M->array[0][p] == old_value)
	  M->array[0][p] = new_value;
    }


    sImage32 *Convert2Image(sMatrix *M){
      sImage32 *img;
      float max,min;
      int p,n;
      
      n = M->ncols*M->nrows;
      img = Image32::Create(M->ncols, M->nrows);
      max = GetMaximumValue(M);
      min = GetMinimumValue(M);
      for(p=0; p<n; p++)
	img->data[p] = ROUND(255.0*(M->array[0][p]-min)/(max-min));
      
      return img;
    }


    sImage32 *Threshold(sMatrix *M,
			float lower, 
			float higher){
      sImage32 *bin=NULL;
      int p,n;
      bin = Image32::Create(M->ncols,M->nrows);
      n = M->ncols*M->nrows;
      for (p=0; p < n; p++)
	if ((M->array[0][p] >= lower)&&(M->array[0][p] <= higher))
	  bin->data[p]=1;
      return(bin);
    }


    // options: 0 (x) / 1 (y) / 2 (z)
    sMatrix* RotationMatrix3(int axis, 
			     float th){
      sMatrix *m;
      m = Create(4,4);
      if (axis==0) {
	m->array[0][0] = 1.0;  m->array[0][1] = 0.0;     m->array[0][2] = 0.0;      m->array[0][3] = 0.0;
	m->array[1][0] = 0.0;  m->array[1][1] = cos(th); m->array[1][2] = -sin(th); m->array[1][3] = 0.0;
	m->array[2][0] = 0.0;  m->array[2][1] = sin(th); m->array[2][2] = cos(th);  m->array[2][3] = 0.0;
	m->array[3][0] = 0.0;  m->array[3][1] = 0.0;     m->array[3][2] = 0.0;      m->array[3][3] = 1.0;
      }
      if (axis==1) {
	m->array[0][0] = cos(th);  m->array[0][1] = 0.0;  m->array[0][2] = sin(th);  m->array[0][3] = 0.0;
	m->array[1][0] = 0.0;      m->array[1][1] = 1;    m->array[1][2] = 0.0;      m->array[1][3] = 0.0;
	m->array[2][0] = -sin(th); m->array[2][1] = 0;    m->array[2][2] = cos(th);  m->array[2][3] = 0.0;
	m->array[3][0] = 0.0;      m->array[3][1] = 0.0;  m->array[3][2] = 0.0;      m->array[3][3] = 1.0;
	
      }
      if (axis==2) {
	m->array[0][0] = cos(th); m->array[0][1] = -sin(th); m->array[0][2] = 0.0;  m->array[0][3] = 0.0;
	m->array[1][0] = sin(th); m->array[1][1] = cos(th);  m->array[1][2] = 0.0;  m->array[1][3] = 0.0;
	m->array[2][0] = 0.0;     m->array[2][1] = 0.0;      m->array[2][2] = 1.0;  m->array[2][3] = 0.0;
	m->array[3][0] = 0.0;     m->array[3][1] = 0.0;      m->array[3][2] = 0.0;  m->array[3][3] = 1.0;
      }
      return m;
    }



    sMatrix* TranslationMatrix3(float dx, float dy, float dz){
      sMatrix *m;
      m = Create(4,4);
      m->array[0][0] = 1.0;  m->array[0][1] = 0.0;  m->array[0][2] = 0.0;  m->array[0][3] = dx;
      m->array[1][0] = 0.0;  m->array[1][1] = 1.0;  m->array[1][2] = 0.0;  m->array[1][3] = dy;
      m->array[2][0] = 0.0;  m->array[2][1] = 0.0;  m->array[2][2] = 1.0;  m->array[2][3] = dz;
      m->array[3][0] = 0.0;  m->array[3][1] = 0.0;  m->array[3][2] = 0.0;  m->array[3][3] = 1.0;
      return m;
    }


    sMatrix* TransformVoxel(sMatrix *m, Voxel v){
      sMatrix *vm,*res;
      vm = Create(1,4);
      vm->array[0][0]=v.c.x;
      vm->array[1][0]=v.c.y;
      vm->array[2][0]=v.c.z;
      vm->array[3][0]=1.0;
      res=Mult(m,vm);
      Destroy(&vm);
      return res;
    }


  } //end Matrix namespace
} //end gft namespace

