
#include "gft_filtering3.h"

namespace gft{

    namespace Kernel3{

    sKernel3 *Create(sAdjRel3 *A){
      sKernel3 *K=NULL;
      int i;
      Voxel max, min;
      
      max.c.x = max.c.y = max.c.z = INT_MIN;
      min.c.x = min.c.y = min.c.z = INT_MAX;
      
      K = (sKernel3 *) calloc(1,sizeof(sKernel3));
      if (K == NULL){
	gft::Error((char *)MSG1,
		   (char *)"Kernel3::Create");
      }
      K->val = gft::AllocFloatArray(A->n);
      K->adj = AdjRel3::Create(A->n);
      
      for (i=0;i<A->n;i++) {
	max.c.x = MAX(A->d[i].axis.x,max.c.x);
	max.c.y = MAX(A->d[i].axis.y,max.c.y);
	max.c.z = MAX(A->d[i].axis.z,max.c.z);
	min.c.x = MIN(A->d[i].axis.x,min.c.x);
	min.c.y = MIN(A->d[i].axis.y,min.c.y);
	min.c.z = MIN(A->d[i].axis.z,min.c.z);
	K->adj->d[i].axis.x = A->d[i].axis.x;
	K->adj->d[i].axis.y = A->d[i].axis.y;
	K->adj->d[i].axis.z = A->d[i].axis.z;
      }
      K->xsize = max.c.x - min.c.x + 1;
      K->ysize = max.c.y - min.c.y + 1;
      K->zsize = max.c.z - min.c.z + 1;
      return(K);
    }


    sKernel3 *Clone(sKernel3 *K){
      sKernel3 *C;
      
      C = (sKernel3 *) calloc(1,sizeof(sKernel3));
      if(C == NULL)
	gft::Error((char *)MSG1,
		   (char *)"Kernel3::Clone");
      
      C->val = gft::AllocFloatArray(K->adj->n);
      memcpy(C->val, K->val,
	     sizeof(float)*K->adj->n);
      C->adj   = AdjRel3::Clone(K->adj);
      C->xsize = K->xsize;
      C->ysize = K->ysize;
      C->zsize = K->zsize;
      return C;
    }

    void Destroy(sKernel3 **K){
      sKernel3 *aux;
      aux = *K;
      if(aux != NULL){
	if (aux->val != NULL) gft::FreeFloatArray(&aux->val); 
	AdjRel3::Destroy(&(aux->adj));
	free(aux);    
	*K = NULL;
      }
    }


    sKernel3 *Normalize(sKernel3 *K){
      sKernel3 *C = Clone(K);
      float wt=0.0;
      int i;
      
      for(i=0; i<K->adj->n; i++)
	wt += K->val[i];
      for(i=0; i<C->adj->n; i++)
	C->val[i] /= wt;
      
      return C;
    }


    sKernel3 *SphericalGaussian(float R, float s, float f){
      float R2,r2;
      sKernel3 *K;
      sAdjRel3 *A;
      int i;
      R2 = R*R;
      A = AdjRel3::Spheric(R);
      K = Create(A);
      for (i=0;i<A->n;i++) {
	r2 = (A->d[i].axis.x * A->d[i].axis.x +
	      A->d[i].axis.y * A->d[i].axis.y +
	      A->d[i].axis.z * A->d[i].axis.z);
	K->val[i] = s*exp (-f*(r2/R2));
      }
      AdjRel3::Destroy(&A);
      return(K);
    }


  } //end Kernel3 namespace

  
  namespace Scene8{


    void ModeFilterLabel(sScene8 *label, float r){
      sAdjRel3 *A = AdjRel3::Spheric(r);
      sAdjVxl  *N;
      int fx,fy,fz,dp;
      int p,q,i,Lmax,l,lmax;
      int *frequency;
      sScene8 *sub=NULL,*mode=NULL,*tmp=NULL;
      struct{
	Voxel v1;
	Voxel v2;
      } box;
      
      MBB(label, &box.v1, &box.v2);
      sub = SubScene(label, box.v1, box.v2);

      AdjRel3::GetFrameSize(A, &fx, &fy, &fz);

      tmp = AddFrame(sub, fx, 0);
      Destroy(&sub);
      sub = tmp;
      mode = Clone(sub);
      
      Lmax = GetMaximumValue(sub);
      frequency = (int *)calloc(Lmax+1,sizeof(int));
      
      N = AdjRel3::AdjVoxels(A, sub);

      dp  = fx*1;
      dp += fy*sub->xsize;
      dp += fz*sub->xsize*sub->ysize;
      
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

      sub = RemFrame(mode, fx);      
      Copy(label, sub, box.v1);
      
      Destroy(&sub);
      Destroy(&mode);
      AdjRel3::Destroy(&A);
      AdjRel3::DestroyAdjVxl(&N);
      free(frequency);
    }
    

  } //end Scene8 namespace



  
  namespace Scene32{


    void ModeFilterLabel(sScene32 *label, float r){
      sAdjRel3 *A = AdjRel3::Spheric(r);
      sAdjVxl  *N;
      int fx,fy,fz,dp;
      int p,q,i,Lmax,l,lmax;
      int *frequency;
      sScene32 *sub=NULL,*mode=NULL,*tmp=NULL;
      struct{
	Voxel v1;
	Voxel v2;
      } box;
      
      MBB(label, &box.v1, &box.v2);
      sub = SubScene(label, box.v1, box.v2);

      AdjRel3::GetFrameSize(A, &fx, &fy, &fz);

      tmp = AddFrame(sub, fx, 0);
      Destroy(&sub);
      sub = tmp;
      mode = Clone(sub);
      
      Lmax = GetMaximumValue(sub);
      frequency = (int *)calloc(Lmax+1,sizeof(int));
      
      N = AdjRel3::AdjVoxels(A, sub);

      dp  = fx*1;
      dp += fy*sub->xsize;
      dp += fz*sub->xsize*sub->ysize;
      
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

      sub = RemFrame(mode, fx);      
      Copy(label, sub, box.v1);
      
      Destroy(&sub);
      Destroy(&mode);
      AdjRel3::Destroy(&A);
      AdjRel3::DestroyAdjVxl(&N);
      free(frequency);
    }
    

    
    sScene32  *AccAbsDiff(sScene32 *scn, float r){
      sScene32 *W;
      sAdjRel3 *A = AdjRel3::Spheric(r);
      int p,q,i,weight;
      Voxel u,v;
      
      W = Create(scn);
      for(p = 0; p < scn->n; p++){
	v.c.x = GetAddressX(scn, p);
	v.c.y = GetAddressY(scn, p);
	v.c.z = GetAddressZ(scn, p);
	
	for(i = 1; i < A->n; i++){
	  u.v = v.v + A->d[i].v;

	  if(gft::Scene32::IsValidVoxel(scn, u)){
	    q = gft::Scene32::GetVoxelAddress(scn, u);
	    weight = abs(scn->data[p] - scn->data[q]);
	    W->data[p] += weight;
	  }
	}
      }
      AdjRel3::Destroy(&A);
      return W;
    }



    int    AreaPercentageHigherThreshold(sCurve *hist,
					 float perc){
      sCurve *nhist;
      double sum=0.0,ratio;
      int i;
      ratio = 1.0 - perc/100.0;
      nhist = Curve::Normalize(hist);
      i = nhist->n-1;
      while(sum<ratio && i>=0){
	sum += nhist->Y[i];
	i--;
      }
      i = MIN(nhist->n-1, i+1);
      Curve::Destroy(&nhist);
      return i;
    }


    void SuppressHighIntensities(sScene32 *scn){
      sCurve *hist;
      int T,p,n;
      hist = Histogram(scn, 1);
      T = AreaPercentageHigherThreshold(hist, 99.95);
      n = scn->xsize*scn->ysize*scn->zsize;
      for(p=0; p<n; p++)
	if(scn->data[p]>T)
	  scn->data[p] = T;
      Curve::Destroy(&hist);
    }


    sScene32 *Convolution(sScene32 *scn, sKernel3 *K){
      sScene32 *cscn;
      Voxel u,v;
      int p,q,i;
      float conv;
      
      cscn = Scene32::Create(scn->xsize,scn->ysize,scn->zsize);
      cscn->dx = scn->dx;
      cscn->dy = scn->dy;
      cscn->dz = scn->dz;
      for(u.c.z=0; u.c.z<scn->zsize; u.c.z++)
	for(u.c.y=0; u.c.y<scn->ysize; u.c.y++)
	  for(u.c.x=0; u.c.x<scn->xsize; u.c.x++){
	    p = GetVoxelAddress(cscn, u);
	    conv=0.0;
	    for (i=0;i<K->adj->n;i++) {
	      v.v = u.v + K->adj->d[i].v;
	      if (Scene32::IsValidVoxel(scn, v)){
		q = GetVoxelAddress(scn, v);
		conv += ((float)scn->data[q])*K->val[i];	   
	      }      
	    }
	    cscn->data[p] = ROUND(conv);
	  }
      return(cscn);
    }


    sScene32 *OptConvolution(sScene32 *scn, sKernel3 *K){
      sScene32 *cscn;
      sAdjVxl *vxl;
      int p,q,i,N,dx,dy,dz,dn;
      float conv;
      N = scn->xsize*scn->ysize*scn->zsize;
      cscn = Scene32::Create(scn->xsize,scn->ysize,scn->zsize);
      cscn->dx = scn->dx;
      cscn->dy = scn->dy;
      cscn->dz = scn->dz;
      vxl = AdjRel3::AdjVoxels(K->adj, scn);
      dx = dy = dz = 0;
      for(i=0;i<K->adj->n;i++){
	dx = MAX(dx, abs(K->adj->d[i].axis.x));
	dy = MAX(dy, abs(K->adj->d[i].axis.y));
	dz = MAX(dz, abs(K->adj->d[i].axis.z));
      }
      dx *= 1;
      dy *= scn->xsize;
      dz *= scn->xsize * scn->ysize;
      dn = dx+dy+dz;
      
      for(p=dn; p<N-dn; p++){
	conv=0.0;
	for(i=0; i<K->adj->n; i++){
	  q = p + vxl->dp[i];
	  conv += ((float)scn->data[q])*K->val[i];	   
	}
	cscn->data[p] = ROUND(conv);
      }
      AdjRel3::DestroyAdjVxl(&vxl);
      return(cscn);
    }


    sScene32   *GaussianBlur(sScene32 *scn){
      sScene32 *blur;
      sKernel3 *K,*NK;
      K  = Kernel3::SphericalGaussian(1.0, 10.0, 1.0); //2.0, 10.0, 1.0
      NK = Kernel3::Normalize(K);
      blur = Convolution(scn, NK);
      blur->dx = scn->dx;
      blur->dy = scn->dy;
      blur->dz = scn->dz;
      Kernel3::Destroy(&NK);
      Kernel3::Destroy(&K);
      return blur;
    }

    /*
    sScene32   *FastGaussianBlur(sScene32 *scn){
      sScene32 *blur;
      sKernel3 *K,*NK;
      K  = Kernel3::SphericalGaussian(1.0, 10.0, 1.0); //2.0, 10.0, 1.0
      NK = Kernel3::Normalize(K);
      blur = FastConvolution(scn, NK);
      blur->dx = scn->dx;
      blur->dy = scn->dy;
      blur->dz = scn->dz;
      Kernel3::Destroy(&NK);
      Kernel3::Destroy(&K);
      return blur;
    }
    */
    
    sScene32   *OptGaussianBlur(sScene32 *scn){
      sScene32 *blur;
      sKernel3 *K,*NK;
      K  = Kernel3::SphericalGaussian(1.0, 10.0, 1.0); //2.0, 10.0, 1.0
      NK = Kernel3::Normalize(K);
      blur = OptConvolution(scn, NK);
      blur->dx = scn->dx;
      blur->dy = scn->dy;
      blur->dz = scn->dz;
      Kernel3::Destroy(&NK);
      Kernel3::Destroy(&K);
      return blur;
    }

    /*
    sScene32   *FastOptGaussianBlur(sScene32 *scn){
      sScene32 *blur;
      sKernel3 *K,*NK;
      K  = Kernel3::SphericalGaussian(1.0, 10.0, 1.0); //2.0, 10.0, 1.0
      NK = Kernel3::Normalize(K);
      blur = FastOptConvolution(scn, NK);
      blur->dx = scn->dx;
      blur->dy = scn->dy;
      blur->dz = scn->dz;
      Kernel3::Destroy(&NK);
      Kernel3::Destroy(&K);
      return blur;
    }
    */

    sScene32   *Subsampling(sScene32 *scn){
      sScene32 *scl;
      Voxel u,v;
      int p,q;
      int xsize,ysize,zsize;
      xsize = (int)(scn->xsize/2) + (scn->xsize%2);
      ysize = (int)(scn->ysize/2) + (scn->ysize%2);
      zsize = (int)(scn->zsize/2) + (scn->zsize%2);
      if(xsize%2==0) xsize++;
      if(ysize%2==0) ysize++;
      if(zsize%2==0) zsize++;
      scl = Scene32::Create(xsize, ysize, zsize);
      scl->dx = 2.0*scn->dx;
      scl->dy = 2.0*scn->dy;
      scl->dz = 2.0*scn->dz;
      
      for(v.c.z=0; v.c.z<scl->zsize; v.c.z++)
	for(v.c.y=0; v.c.y<scl->ysize; v.c.y++)
	  for(v.c.x=0; v.c.x<scl->xsize; v.c.x++){
	    p = Scene32::GetVoxelAddress(scl, v);
	    u.c.x = (v.c.x-scl->xsize/2)*2 + scn->xsize/2;
	    u.c.y = (v.c.y-scl->ysize/2)*2 + scn->ysize/2;
	    u.c.z = (v.c.z-scl->zsize/2)*2 + scn->zsize/2;
	    
	    if(Scene32::IsValidVoxel(scn, u)){
	      q = Scene32::GetVoxelAddress(scn, u);
	      scl->data[p] = scn->data[q];
	    }
	  }
      
      return(scl);
    }


    sScene32 *LaplacianFilter(sScene32 *orig){
      sScene32 *lap;
      int dx,dy,dz,di,m,i,N;
      int acc;
      static int n0[9]={ -1,  0,  1, -1,  0,  1, -1,  0,  1};
      static int n1[9]={ -1, -1, -1,  0,  0,  0,  1,  1,  1};
      static int wt[9]={ -1, -1, -1, -1, -1, -1, -1, -1, -1};
      static int wc[9]={ -1, -1, -1, -1, 26, -1, -1, -1, -1};
      
      lap = Create(orig);
      N = orig->xsize*orig->ysize*orig->zsize;
      dx = 1;
      dy = orig->xsize;
      dz = orig->xsize * orig->ysize;
      di = dx+dy+dz;
      for(i=di; i<N-di; i++){
	// xy plane
	acc = 0;
	for(m=0;m<9;m++){
	  acc += wt[m]*orig->data[ i+dz+dx*n0[m]+dy*n1[m] ];
	  acc += wc[m]*orig->data[ i+ 0+dx*n0[m]+dy*n1[m] ];
	  acc += wt[m]*orig->data[ i-dz+dx*n0[m]+dy*n1[m] ];
	}
	lap->data[i] = acc;
      }
      return lap;
    }
    
    
    
    sScene32 *SobelFilter(sScene32 *scn){
      sScene32 *grad;
      int   dx,dy,dz,di,m,i,N;
      int   acc[3];
      float facc[3];
      static int n0[9]={-1, 0, 1,-1, 0, 1,-1, 0, 1};
      static int n1[9]={-1,-1,-1, 0, 0, 0, 1, 1, 1};
      static int we[9]={ 1, 2, 1, 2, 4, 2, 1, 2, 1};
      
      grad = Create(scn);
      N = scn->xsize*scn->ysize*scn->zsize;
      dx = 1;
      dy = scn->xsize;
      dz = scn->xsize * scn->ysize;
      di = dx+dy+dz;
      
      for(i=di; i<N-di; i++){
	// xy plane
	acc[0] = 0;
	for(m=0; m<9; m++)
	  acc[0] += we[m] * (scn->data[ i+dz+dx*n0[m]+dy*n1[m] ] -
			     scn->data[ i-dz+dx*n0[m]+dy*n1[m] ]);
	// yz plane
	acc[1] = 0;
	for(m=0; m<9; m++)
	  acc[1] += we[m] * (scn->data[ i+dx+dz*n0[m]+dy*n1[m] ] -
			     scn->data[ i-dx+dz*n0[m]+dy*n1[m] ]);
	// xz plane
	acc[2] = 0;
	for(m=0; m<9; m++)
	  acc[2] += we[m] * (scn->data[ i+dy+dx*n0[m]+dz*n1[m] ] -
			     scn->data[ i-dy+dx*n0[m]+dz*n1[m] ]);
	acc[0] >>= 4;
	acc[1] >>= 4;
	acc[2] >>= 4;
	facc[0] = (float)acc[0];
	facc[1] = (float)acc[1];
	facc[2] = (float)acc[2];
	grad->data[i] = ROUND(sqrtf((facc[0]*facc[0] + 
				     facc[1]*facc[1] +
				     facc[2]*facc[2])/3.0));
      }
      return grad;
    }
    

    sScene32 *SphericalGradient(sScene32 *scn, float r){
      sAdjRel3 *A = AdjRel3::Spherical_mm(scn, r);
      sAdjVxl  *N;
      sScene32 *grad=NULL;
      int p,q,i,vp,vq;
      float *mg=NULL;
      float gx,gy,gz,d;
      int fx,fy,fz,dp;
      N  = AdjRel3::AdjVoxels(A, scn);
      mg = AdjRel3::GetDistanceArray_mm(A, scn);
      AdjRel3::GetFrameSize(A, &fx, &fy, &fz);
      dp  = fx*1;
      dp += fy*scn->xsize;
      dp += fz*scn->xsize*scn->ysize;
      grad = Create(scn);
      for(p=dp; p<scn->n-dp; p++){
	vp = scn->data[p];
	gx = gy = gz = 0.0;
	for(i=1; i<N->n; i++){
	  q = p + N->dp[i];
	  vq = scn->data[q];
	  d = (float)(vq-vp);
	  
	  gx  += (d*A->d[i].axis.x)/(mg[i]);
	  gy  += (d*A->d[i].axis.y)/(mg[i]);
	  gz  += (d*A->d[i].axis.z)/(mg[i]);
	}
	gx = ROUND(10.0*gx*scn->dx/fx);
	gy = ROUND(10.0*gy*scn->dy/fy);
	gz = ROUND(10.0*gz*scn->dz/fz);
	grad->data[p] = ROUND(sqrtf(gx*gx + gy*gy + gz*gz));
      }
      //Scene32::ClearAdjFrame(grad, A);
      AdjRel3::Destroy(&A);
      AdjRel3::DestroyAdjVxl(&N);
      free(mg);
      return grad;
    }

    
    
  } //end Scene32 namespace

  

} //end gft namespace





