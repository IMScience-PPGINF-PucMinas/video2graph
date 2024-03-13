
#include "gft_filtering.h"

namespace gft{

  namespace Kernel{

    sKernel *Make(char *coefs){
      sKernel *K;
      sAdjRel *A;
      int xsize,ysize,i;
      float val;

      sscanf(coefs,"%d",&xsize);
      coefs=strchr(coefs,',')+1;
      sscanf(coefs,"%d",&ysize);
      coefs=strchr(coefs,',')+1;
      
      A = AdjRel::Box(xsize, ysize);
      K = Create(A);
      for(i = 0; i < A->n; i++){
	sscanf(coefs, "%f", &K->val[i]);
	coefs = strchr(coefs,',')+1;
      }

      /* Put the middle value (corresponding to the origin) at the first
	 place to match the vector of coefficients with the vector of
	 displacements (adjacency relation) */
      
      for(i=A->n/2; i > 0; i--){
	val = K->val[i];
	K->val[i] = K->val[i-1];
	K->val[i-1] = val;
      }
      
      AdjRel::Destroy(&A);
      return(K);
    }

    
    sKernel *Create(sAdjRel *A){
      sKernel *K=NULL;
      int i;
      
      K = (sKernel *) calloc(1,sizeof(sKernel));
      if(K == NULL){
	gft::Error(MSG1, "Kernel::Create");
      }
      K->val = gft::AllocFloatArray(A->n);
      K->adj = AdjRel::Clone(A);
      return(K);
    }


    sKernel *Clone(sKernel *K){
      sKernel *C;
      
      C = (sKernel *) calloc(1,sizeof(sKernel));
      if(C == NULL)
	gft::Error(MSG1, "Kernel::Clone");
      
      C->val = gft::AllocFloatArray(K->adj->n);
      memcpy(C->val, K->val,
	     sizeof(float)*K->adj->n);
      C->adj = AdjRel::Clone(K->adj);
      return C;
    }


    sKernel *Normalize(sKernel *K){
      sKernel *C = Clone(K);
      float wt=0.0;
      int i;
      for(i = 0; i < K->adj->n; i++)
	wt += K->val[i];
      for(i = 0; i < C->adj->n; i++)
	C->val[i] /= wt;
      
      return C;
    }


    void    Destroy(sKernel **K){
      sKernel *aux;
      aux = *K;
      if(aux != NULL){
	if (aux->val != NULL) gft::FreeFloatArray(&(aux->val));
	AdjRel::Destroy(&(aux->adj));
	free(aux);
	*K = NULL;
      }
    }


    sKernel *Gaussian(sAdjRel *A, float stddev){
      double k,k1,d2,sigma,sigma2;
      sKernel *K;
      int i;
      d2 = 0.0;
      for (i=0;i<A->n;i++) {
	d2 = MAX((double)(A->dx[i]*A->dx[i]+A->dy[i]*A->dy[i]),d2);
      }
      sigma2 = d2 / (stddev*stddev);
      sigma = sqrt(sigma2);
      k = 2.0*sigma2;
      k1 = 1.0/(sqrt(2*PI)*sigma);
      K = Create(A);
      for(i = 0; i < A->n; i++){
	d2 = (double)(A->dx[i] * A->dx[i] + A->dy[i] * A->dy[i]);
	K->val[i] = (float)(k1*exp(-d2/k)); // gaussian
      }
      return K;
    }


  } //end Kernel namespace


  namespace Image32{


    sImage32 *GaussianBlur(sImage32 *img,
			   float stddev){
      sImage32 *blur;
      sAdjRel *A;
      sKernel *K,*N;
      A = AdjRel::Circular(2.0);
      K = Kernel::Gaussian(A, stddev);
      N = Kernel::Normalize(K);
      blur = LinearFilter(img, N);
      AdjRel::Destroy(&A);
      Kernel::Destroy(&K);
      Kernel::Destroy(&N);
      return blur;
    }


    sImage32 *GaussianBlur(sImage32 *img){
      sImage32 *blur;
      sKernel *K,*N;
      K = Kernel::Make("5,5, 1, 4, 7, 4, 1, 4,16,26,16, 4, 7,26,41,26, 7, 4,16,26,16, 4, 1, 4, 7, 4, 1");
      N = Kernel::Normalize(K);
      blur = LinearFilter(img, N);
      Kernel::Destroy(&K);
      Kernel::Destroy(&N);
      return blur;
    }
    
    /*
    void ModeFilterLabel(sImage32 *label, float r){
    }
    */

    sImage32 *SobelFilter(sImage32 *img){
      sImage32 *gradx=NULL,*grady=NULL,*grad=NULL;
      sKernel *Kx,*Ky;

      Ky = Kernel::Make("3,3,-1.0,-2.0,-1.0,0.0,0.0,0.0,1.0,2.0,1.0");
      Kx = Kernel::Make("3,3,-1.0,0.0,1.0,-2.0,0.0,2.0,-1.0,0.0,1.0");
      gradx = LinearFilter(img, Kx);
      grady = LinearFilter(img, Ky);
      grad  = ImageMagnitude(gradx, grady);
      Destroy(&gradx);
      Destroy(&grady);
      Kernel::Destroy(&Kx);
      Kernel::Destroy(&Ky);
      return(grad);
    }
    

    sImage32 *LinearFilter(sImage32 *img, sKernel *K){
      sImage32 *cimg;
      int u_x,u_y;
      float conv;
      unsigned int p, i, x, y;
      
      cimg = Create(img->ncols, img->nrows);
      
      for(u_y = 0; u_y < img->nrows; ++u_y){
	for(u_x = 0; u_x < img->ncols; ++u_x){
	  conv = 0.0;
	  for(i = 0; i < K->adj->n; ++i){
	    x = u_x + K->adj->dx[i];
	    y = u_y + K->adj->dy[i];
	    if ((x >= 0) && (x < img->ncols) && (y >= 0) && (y < img->nrows))
	      conv += (float)img->array[y][x] * K->val[i];
	  }
	  cimg->array[u_y][u_x] = (int)conv;
	}
      }
      return(cimg);
    }


    sImage32  *ImageMagnitude(sImage32 *imgx, sImage32 *imgy){
      sImage32 *mag = Create(imgx->ncols, imgx->nrows);
      int p, n = imgx->ncols*imgx->nrows;
      
      for(p = 0; p < n; p++)
        mag->data[p] = ROUND(sqrt(imgx->data[p]*imgx->data[p] + imgy->data[p]*imgy->data[p]));
 
      return(mag);
    }


    /**
     * Comparison function to be used in QuickSort
     */
    static int comp(const void *a, const void *b) {
      return *(int*)a - *(int*)b;
    }

    sImage32 *MedianFilter(sImage32 *img, sAdjRel *A){
      int *val,*oval,n,i,p,q;
      Pixel u,v;
      sImage32 *med;
      val = gft::AllocIntArray(A->n);
      med = Create(img);
      for (u.y=0; u.y < img->nrows; u.y++){
	for (u.x=0; u.x < img->ncols; u.x++) {
	  p = u.x + u.y*img->ncols;
	  n = 0;
	  for (i=0; i < A->n; i++) {
	    v.x = u.x + A->dx[i];
	    v.y = u.y + A->dy[i];
	    if (IsValidPixel(img,v.x,v.y)){
	      q = v.x + v.y*img->ncols;
	      val[n] = img->data[q];
	      n++;
	    }
	  }
	  qsort(val, n, sizeof(int), comp);
	  med->data[p] = val[n/2];
	}
      }
      gft::FreeIntArray(&val);
      return(med);
    }


    void ModeFilterLabel(sImage32 *label, float r){
      sAdjRel *A = AdjRel::Circular(r);
      sAdjPxl *N;
      sImage32 *sub = NULL, *mode = NULL, *tmp = NULL;
      int f,Lmax,dp,p,q,i,l,lmax;
      int *frequency;
      struct{
	Pixel p1;
	Pixel p2;
      } box;
      
      MBB(label, &box.p1, &box.p2);
      sub = Clone(label, box.p1, box.p2);

      f = AdjRel::GetFrameSize(A);
      
      tmp = AddFrame(sub, f, 0);
      Destroy(&sub);
      sub = tmp;
      mode = Clone(sub);

      Lmax = GetMaxVal(sub);
      frequency = (int *)calloc(Lmax+1,sizeof(int));

      N = AdjRel::AdjPixels(A, sub);

      dp  = f*1;
      dp += f*sub->ncols;

      for(p=dp; p<sub->n-dp; p++){
	
	memset(frequency, 0, (Lmax+1)*sizeof(int));
	
	for(i=0; i<N->n; i++){
	  q = p + N->dp[i];
	  frequency[sub->data[q]]++;
	}
	lmax = sub->data[p];
	for(l=0; l<=Lmax; l++){
	  if(frequency[l]>frequency[lmax])
	    lmax = l;
	}
	mode->data[p] = lmax;
      }
      Destroy(&sub);

      sub = RemFrame(mode, f);
      Copy(label, sub, box.p1);

      Destroy(&sub);
      Destroy(&mode);
      AdjRel::Destroy(&A);
      AdjRel::DestroyAdjPxl(&N);
      free(frequency);
    }
    

  } //end Image32 namespace





  namespace Image32f{  

    sImage32f *SobelFilter(sImage32f *img){
      sImage32f *gradx=NULL,*grady=NULL,*grad=NULL;
      sKernel *Kx,*Ky;
      
      Ky = Kernel::Make("3,3,-1.0,-2.0,-1.0,0.0,0.0,0.0,1.0,2.0,1.0");
      Kx = Kernel::Make("3,3,-1.0,0.0,1.0,-2.0,0.0,2.0,-1.0,0.0,1.0");
      gradx = LinearFilter(img, Kx);
      grady = LinearFilter(img, Ky);
      grad  = ImageMagnitude(gradx, grady);
      Destroy(&gradx);
      Destroy(&grady);
      Kernel::Destroy(&Kx);
      Kernel::Destroy(&Ky);
      return(grad);
    }

    sImage32f *LinearFilter(sImage32f *img, sKernel *K){
      sImage32f *cimg;
      int u_x,u_y;
      float conv;
      unsigned int p, i, x, y;
      
      cimg = Create(img->ncols, img->nrows);
      
      for(u_y = 0; u_y < img->nrows; ++u_y){
	for(u_x = 0; u_x < img->ncols; ++u_x){
	  conv = 0.0;
	  for(i = 0; i < K->adj->n; ++i){
	    x = u_x + K->adj->dx[i];
	    y = u_y + K->adj->dy[i];
	    if ((x >= 0) && (x < img->ncols) && (y >= 0) && (y < img->nrows))
	      conv += img->array[y][x] * K->val[i];
	  }
	  cimg->array[u_y][u_x] = conv;
	}
      }
      return(cimg);
    }


    sImage32f  *ImageMagnitude(sImage32f *imgx, sImage32f *imgy){
      sImage32f *mag = Create(imgx->ncols, imgx->nrows);
      int p, n = imgx->ncols*imgx->nrows;
      
      for(p = 0; p < n; p++)
        mag->data[p] = sqrt(imgx->data[p]*imgx->data[p] +
			    imgy->data[p]*imgy->data[p]);
      
      return(mag);
    }
    
    
  } //end Image32f namespace

  
  namespace CImage{

    sImage32 *SobelFilter(sCImage *cimg){
      sImage32 *g0,*g1,*g2,*grad;
      int p;
      g0 = gft::Image32::SobelFilter(cimg->C[0]);
      g1 = gft::Image32::SobelFilter(cimg->C[1]);
      g2 = gft::Image32::SobelFilter(cimg->C[2]);
      grad = gft::Image32::Create(cimg->C[0]);
      for(p = 0; p < cimg->C[0]->n; p++){
	/*
	grad->data[p] = MAX(g0->data[p],
			    MAX(g1->data[p], g2->data[p]));
	*/
	grad->data[p] = ROUND(sqrtf(g0->data[p]*g0->data[p] +
				    g1->data[p]*g1->data[p] +
				    g2->data[p]*g2->data[p]));
      }
      gft::Image32::Destroy(&g0);
      gft::Image32::Destroy(&g1);
      gft::Image32::Destroy(&g2);
      return grad;
    }
    
  } //end CImage namespace



  namespace CImage32f{

    sImage32f *SobelFilter(sCImage32f *cimg){
      sImage32f *g0,*g1,*g2,*grad;
      int p;
      g0 = gft::Image32f::SobelFilter(cimg->C[0]);
      g1 = gft::Image32f::SobelFilter(cimg->C[1]);
      g2 = gft::Image32f::SobelFilter(cimg->C[2]);
      grad = gft::Image32f::Create(cimg->C[0]);
      for(p = 0; p < cimg->C[0]->n; p++){
	/*
	grad->data[p] = MAX(g0->data[p],
			    MAX(g1->data[p], g2->data[p]));
	*/
	grad->data[p] = sqrtf(g0->data[p]*g0->data[p] +
			      g1->data[p]*g1->data[p] +
			      g2->data[p]*g2->data[p]);
      }
      gft::Image32f::Destroy(&g0);
      gft::Image32f::Destroy(&g1);
      gft::Image32f::Destroy(&g2);
      return grad;
    }
    
  } //end CImage32f namespace
  
  
} //end gft namespace

