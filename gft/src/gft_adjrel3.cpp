
#include "gft_adjrel3.h"

namespace gft{
  namespace AdjRel3{

    sAdjRel3 *Create(int n){
      sAdjRel3 *A=NULL;

      A = (sAdjRel3 *) calloc(1,sizeof(sAdjRel3));
      if(A != NULL){
	A->d = (sDisplacement3 *)_mm_malloc(sizeof(sDisplacement3)*n, 16);
	if(A->d==NULL)
	  Error((char *)MSG1,(char *)"AdjRel3::Create");
	memset((void *)A->d, 0, sizeof(sDisplacement3)*n);
	A->n = n;
      }
      else
	Error((char *)MSG1,(char *)"AdjRel3::Create");
      return(A);
    }


    void     Destroy(sAdjRel3 **A){
      sAdjRel3 *aux;
      aux = *A;
      if(aux != NULL){
	if(aux->d != NULL) _mm_free(aux->d);
	free(aux);
	*A = NULL;
      }
    }


    sAdjRel3 *Clone(sAdjRel3 *A){
      sAdjRel3 *C;
      int i;

      C = Create(A->n);
      for(i=0; i < A->n; i++){
	C->d[i] = A->d[i];
      }
      return C;
    }

    
    sAdjRel3 *Spheric(float r){
      sAdjRel3 *A=NULL;
      int i,n,r0,r2,dx,dy,dz,i0=0;
      n=0;
      r0 = (int)r;
      r2  = (int)(r*r);
      for(dz=-r0;dz<=r0;dz++)
	for(dy=-r0;dy<=r0;dy++)
	  for(dx=-r0;dx<=r0;dx++)
	    if(((dx*dx)+(dy*dy)+(dz*dz)) <= r2)
	      n++;
	
      A = Create(n);
      i=0;
      for(dz=-r0;dz<=r0;dz++)
	for(dy=-r0;dy<=r0;dy++)
	  for(dx=-r0;dx<=r0;dx++)
	    if(((dx*dx)+(dy*dy)+(dz*dz)) <= r2){
	      A->d[i].axis.x=dx;
	      A->d[i].axis.y=dy;
	      A->d[i].axis.z=dz;
	      if ((dx==0)&&(dy==0)&&(dz==0))
		i0 = i;
	      i++;	  
	    }
      /* shift to right and place central voxel at first */
        for (i=i0; i > 0; i--) {
	  dx = A->d[i].axis.x;
	  dy = A->d[i].axis.y;
	  dz = A->d[i].axis.z;
	  A->d[i].axis.x = A->d[i-1].axis.x;
	  A->d[i].axis.y = A->d[i-1].axis.y;
	  A->d[i].axis.z = A->d[i-1].axis.z;
	  A->d[i-1].axis.x = dx;
	  A->d[i-1].axis.y = dy;
	  A->d[i-1].axis.z = dz;
	}
	return(A);
    }


    sAdjRel3 *Ellipsoid(float rx, float ry, float rz){
      sAdjRel3 *A=NULL;
      int i,n,dx,dy,dz,i0=0;
      float rx2,ry2,rz2;
      n=0;
      rx2  = rx * rx;
      ry2  = ry * ry;
      rz2  = rz * rz;
      for(dz=-rz;dz<=rz;dz++)
	for(dy=-ry;dy<=ry;dy++)
	  for(dx=-rx;dx<=rx;dx++)
	    if(((dx*dx)/rx2+(dy*dy)/ry2+(dz*dz)/rz2) <= 1.001)
	      n++;
  
      A = Create(n);

      i=0;
      for(dz=-rz;dz<=rz;dz++)
	for(dy=-ry;dy<=ry;dy++)
	  for(dx=-rx;dx<=rx;dx++)
	    if(((dx*dx)/rx2+(dy*dy)/ry2+(dz*dz)/rz2) <= 1.001) {
	      A->d[i].axis.x=dx;
	      A->d[i].axis.y=dy;
	      A->d[i].axis.z=dz;
	      if ((dx==0)&&(dy==0)&&(dz==0))
		i0 = i;
	      i++;	  
	    }
      /* shift to right and place central voxel at first */
      for (i=i0; i > 0; i--) {
	dx = A->d[i].axis.x;
	dy = A->d[i].axis.y;
	dz = A->d[i].axis.z;
	A->d[i].axis.x = A->d[i-1].axis.x;
	A->d[i].axis.y = A->d[i-1].axis.y;
	A->d[i].axis.z = A->d[i-1].axis.z;
	A->d[i-1].axis.x = dx;
	A->d[i-1].axis.y = dy;
	A->d[i-1].axis.z = dz;
      }
      return(A);
    }


    sAdjRel3 *SphericalShell(float inner_radius,
			     float outer_radius){
      sAdjRel3 *A=NULL;
      int i,n,dx,dy,dz,r0i,r2i,r0o,r2o,d;

      if(outer_radius <= inner_radius)
	Error((char *)"outer_radius must be greater than inner_radius",
	      (char *)"AdjRel3::SphericalShell");
      n=0;
      r0i  = (int)inner_radius;
      r2i  = (int)(inner_radius*inner_radius);
      r0o  = (int)outer_radius;
      r2o  = (int)(outer_radius*outer_radius);

      for(dz=-r0o;dz<=r0o;dz++){
	for(dy=-r0o;dy<=r0o;dy++){
	  for(dx=-r0o;dx<=r0o;dx++){
	    d = (dx*dx)+(dy*dy)+(dz*dz);
	    if((d <= r2o)&&(d>=r2i))
	      n++;
	  }
	}
      }
      
      A = Create(n);
      i=0;

      for(dz=-r0o;dz<=r0o;dz++){
	for(dy=-r0o;dy<=r0o;dy++){
	  for(dx=-r0o;dx<=r0o;dx++){
	    d = (dx*dx)+(dy*dy)+(dz*dz);
	    if((d <= r2o)&&(d>=r2i)){
	      A->d[i].axis.x=dx;
	      A->d[i].axis.y=dy;
	      A->d[i].axis.z=dz;
	      i++;
	    }
	  }
	}
      }
      return(A);
    }


    sAdjRel3 *Box(int xsize, int ysize, int zsize){
      sAdjRel3 *A=NULL;
      int i,dx,dy,dz;

      if(xsize%2 == 0) xsize++;
      if(ysize%2 == 0) ysize++;
      if(zsize%2 == 0) zsize++;
      
      A = Create(xsize*ysize*zsize);
      i=1;
      for(dz=-zsize/2;dz<=zsize/2;dz++)
	for(dy=-ysize/2;dy<=ysize/2;dy++)
	  for(dx=-xsize/2;dx<=xsize/2;dx++)
	    if (!((dx == 0)&&(dy == 0)&&(dz == 0))){
	      A->d[i].axis.x=dx;
	      A->d[i].axis.y=dy;
	      A->d[i].axis.z=dz;
	      i++;
	    }
      /* place the central pixel at first */
      A->d[0].axis.x = 0;
      A->d[0].axis.y = 0;
      A->d[0].axis.z = 0;
      return(A);
    }


    void     Scale(sAdjRel3 *A,
		   float Sx, float Sy, float Sz){
      int i;
      for(i=0; i<A->n; i++){
	A->d[i].axis.x = ROUND((float)A->d[i].axis.x*Sx);
	A->d[i].axis.y = ROUND((float)A->d[i].axis.y*Sy);
	A->d[i].axis.z = ROUND((float)A->d[i].axis.z*Sz);
      }
    }


    void     ClipX(sAdjRel3 *A, int lower, int higher){
      sAdjRel3 *B=NULL;
      int i,j,n=0;
      sDisplacement3 *tmp;

      for(i=0; i<A->n; i++)
	if(A->d[i].axis.x>=lower && A->d[i].axis.x<=higher) 
	  n++;

      B = Create(n);
      j = 0;
      for(i=0; i<A->n; i++){
	if(A->d[i].axis.x>=lower && A->d[i].axis.x<=higher){
	  B->d[j].axis.x = A->d[i].axis.x;
	  B->d[j].axis.y = A->d[i].axis.y;
	  B->d[j].axis.z = A->d[i].axis.z;
	  j++;
	}
      }
      A->n = B->n;
      tmp  = A->d;
      A->d = B->d;
      B->d = tmp;
      Destroy(&B);
    }


    void     ClipY(sAdjRel3 *A, int lower, int higher){
      sAdjRel3 *B=NULL;
      int i,j,n=0;
      sDisplacement3 *tmp;

      for(i=0; i<A->n; i++)
	if(A->d[i].axis.y>=lower && A->d[i].axis.y<=higher) 
	  n++;

      B = Create(n);
      j = 0;
      for(i=0; i<A->n; i++){
	if(A->d[i].axis.y>=lower && A->d[i].axis.y<=higher){
	  B->d[j].axis.x = A->d[i].axis.x;
	  B->d[j].axis.y = A->d[i].axis.y;
	  B->d[j].axis.z = A->d[i].axis.z;
	  j++;
	}
      }
      A->n = B->n;
      tmp  = A->d;
      A->d = B->d;
      B->d = tmp;
      Destroy(&B);
    }


    void     ClipZ(sAdjRel3 *A, int lower, int higher){
      sAdjRel3 *B=NULL;
      int i,j,n=0;
      sDisplacement3 *tmp;

      for(i=0; i<A->n; i++)
	if(A->d[i].axis.z>=lower && A->d[i].axis.z<=higher) 
	  n++;

      B = Create(n);
      j = 0;
      for(i=0; i<A->n; i++){
	if(A->d[i].axis.z>=lower && A->d[i].axis.z<=higher){
	  B->d[j].axis.x = A->d[i].axis.x;
	  B->d[j].axis.y = A->d[i].axis.y;
	  B->d[j].axis.z = A->d[i].axis.z;
	  j++;
	}
      }
      A->n = B->n;
      tmp  = A->d;
      A->d = B->d;
      B->d = tmp;
      Destroy(&B);
    }

    
    int      GetFrameSize(sAdjRel3 *A){
      int sz=INT_MIN,i=0;

      for (i=0; i < A->n; i++){
	if (abs(A->d[i].axis.x) > sz) sz = abs(A->d[i].axis.x);
	if (abs(A->d[i].axis.y) > sz) sz = abs(A->d[i].axis.y);
	if (abs(A->d[i].axis.z) > sz) sz = abs(A->d[i].axis.z);
      }
      return(sz);
    }


    void     GetFrameSize(sAdjRel3 *A, 
			  int *sz_x, 
			  int *sz_y, 
			  int *sz_z){
      int i=0;
      *sz_x = INT_MIN;
      *sz_y = INT_MIN;
      *sz_z = INT_MIN;
      for(i=0; i<A->n; i++){
	if(abs(A->d[i].axis.x) > *sz_x) *sz_x = abs(A->d[i].axis.x);
	if(abs(A->d[i].axis.y) > *sz_y) *sz_y = abs(A->d[i].axis.y);
	if(abs(A->d[i].axis.z) > *sz_z) *sz_z = abs(A->d[i].axis.z);
      }
    }


    float   *GetDistanceArray(sAdjRel3 *A){
      float *mg=NULL;
      int i;
      mg = gft::AllocFloatArray(A->n);
      for(i=0; i<A->n; i++){
	mg[i] = sqrt(A->d[i].axis.x*A->d[i].axis.x +
		     A->d[i].axis.y*A->d[i].axis.y +
		     A->d[i].axis.z*A->d[i].axis.z);
      }
      return mg;
    }


    float   *GetDistanceArray_mm(sAdjRel3 *A, sScene32 *scn){
      float *mg=NULL;
      float dx,dy,dz;
      int i;
      mg = gft::AllocFloatArray(A->n);
      for(i=0; i<A->n; i++){
	dx = A->d[i].axis.x*scn->dx;
	dy = A->d[i].axis.y*scn->dy;
	dz = A->d[i].axis.z*scn->dz;
	mg[i] = sqrt(dx*dx + dy*dy + dz*dz);
      }
      return mg;
    }


    sAdjVxl  *AdjVoxels(sAdjRel3 *A, int xsize, int ysize){
      sAdjVxl *N;
      int i,xysize;

      xysize = xsize * ysize;

      N = (sAdjVxl *) calloc(1,sizeof(sAdjVxl));
      if(N != NULL){
	N->dp = gft::AllocIntArray(A->n);
	N->n  = A->n;
	for (i=0; i < N->n; i++)
	  N->dp[i] = (A->d[i].axis.x + 
		      A->d[i].axis.y * xsize + 
		      A->d[i].axis.z * xysize);
      }
      else
	Error((char *)MSG1,(char *)"AdjRel3::AdjVoxels");  
      return(N);
    }


    sAdjVxl  *AdjVoxels(sAdjRel3 *A, sScene32 *scn){
      return AdjVoxels(A, scn->xsize, scn->ysize);
    }

    sAdjVxl  *AdjVoxels(sAdjRel3 *A, sScene16 *scn){
      return AdjVoxels(A, scn->xsize, scn->ysize);
    }

    sAdjVxl  *AdjVoxels(sAdjRel3 *A, sScene8 *scn){
      return AdjVoxels(A, scn->xsize, scn->ysize);
    }

    sAdjVxl  *AdjVoxels(sAdjRel3 *A, sScene *scn){
      switch(scn->nbits){
      case 8:
	return AdjVoxels(A, scn->ptr.scn8);
      case 16:
	return AdjVoxels(A, scn->ptr.scn16);
      case 32:
	return AdjVoxels(A, scn->ptr.scn32);
      }
      return NULL;
    }


    void     DestroyAdjVxl(sAdjVxl **N){
      sAdjVxl *aux;
      aux = *N;
      if(aux != NULL){
	if(aux->dp != NULL) gft::FreeIntArray(&aux->dp);
	free(aux);
	*N = NULL;
      }
    }


    sAdjRel3 *Spherical_mm(sScene32 *scn, float r_mm){
      sAdjRel3 *A=NULL;
      float rx,ry,rz;
      rx = MAX(r_mm/scn->dx, 1.0);
      ry = MAX(r_mm/scn->dy, 1.0);
      rz = MAX(r_mm/scn->dz, 1.0);
      if(scn->dx==scn->dy && scn->dx==scn->dz)
	A = Spheric(rx);
      else
	A = Ellipsoid(rx, ry, rz);
      return A;
    }


    sAdjRel3 *SphericalGrid_mm(float dx, float dy, float dz,
			       float r_mm,
			       float spacement_mm){
      sAdjRel3 *A=NULL;
      int ri,i;

      ri = ROUND(r_mm/spacement_mm);
      A = Spheric((float)ri);

      for(i=1; i<A->n; i++){
	A->d[i].axis.x = ROUND((float)A->d[i].axis.x*
			       spacement_mm/dx);
	A->d[i].axis.y = ROUND((float)A->d[i].axis.y*
			       spacement_mm/dy);
	A->d[i].axis.z = ROUND((float)A->d[i].axis.z*
			       spacement_mm/dz);
      }
      return A;
    }


    sAdjRel3 *SphericalGrid_mm(sScene32 *scn, float r_mm,
			       float spacement_mm){
      return SphericalGrid_mm(scn->dx, scn->dy, scn->dz,
			      r_mm, spacement_mm);
    }


  } //end AdjRel3 namespace
} //end gft namespace


