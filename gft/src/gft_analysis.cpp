
#include "gft_analysis.h"

namespace gft{
  namespace Image32{


    sImage32 *RegMin(sImage32 *img, sAdjRel *A){
      sImage32 *Rmin;
      sStack *Q;
      int r,p,q,px,py,qx,qy,i,val,lb = 0;
      bool isminima;
      int *P;
      Rmin = Create(img);
      P = gft::AllocIntArray(img->n);
      Q = gft::Stack::Create(img->n);
      for(r = 0; r < img->n; r++){
	if(P[r] > 0) continue;
	isminima = true;
	val = img->data[r];
	
	P[r] = 1;
	gft::Stack::Push(Q, r);
	while(!gft::Stack::IsEmpty(Q)){
	  p = gft::Stack::Pop(Q);
	  px = p%img->ncols;
	  py = p/img->ncols;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    if(gft::Image32::IsValidPixel(img, qx,qy)){
	      q = qx + qy*img->ncols;
	      if(img->data[q] == val && P[q] == 0){
		P[q] = 1;
		gft::Stack::Push(Q, q);
	      }
	      else if(img->data[q] < val){
		isminima = false;
	      }
	    }
	  }
	}

	if(!isminima) continue;

	lb++;
	Rmin->data[r] = lb;
	gft::Stack::Push(Q, r);
	while(!gft::Stack::IsEmpty(Q)){
	  p = gft::Stack::Pop(Q);
	  px = p%img->ncols;
	  py = p/img->ncols;
	  for(i=1; i < A->n; i++){
	    qx = px + A->dx[i];
	    qy = py + A->dy[i];
	    if(gft::Image32::IsValidPixel(img, qx,qy)){
	      q = qx + qy*img->ncols;
	      if(img->data[q] == val && Rmin->data[q] == 0){
		Rmin->data[q] = lb;
		gft::Stack::Push(Q, q);
	      }
	    }
	  }
	}
      }

      /*
      printf("lb: %d\n", lb);
      gft::Image32::Write(Rmin, (char *)"regmin.pgm");

      gft::CImage::CImage *ctmp;
      ctmp = gft::CImage::RandomColorize(Rmin);
      gft::CImage::Write(ctmp, (char *)"regmin.ppm");
      gft::CImage::Destroy(&ctmp);
      */
      
      gft::Stack::Destroy(&Q);
      gft::FreeIntArray(&P);
      return Rmin;
    }
    

    sImage32 *LabelBinComp(sImage32 *bin, sAdjRel *A){
      sImage32 *label=NULL;
      int i,j,n,p,q,l=1;
      int px,py,qx,qy;
      int *FIFO=NULL;
      int first=0,last=0;
      
      label = Create(bin->ncols, bin->nrows);
      n  = bin->ncols*bin->nrows;
      FIFO  = gft::AllocIntArray(n);
      for(j=0; j < n; j++){
	if((bin->data[j] != 0)&&(label->data[j]==0)){
	  label->data[j]=l;
	  FIFO[last]=j;
	  last++;
	  while(first != last){
	    p = FIFO[first];
	    px = p%bin->ncols;
	    py = p/bin->ncols;
	    first++;
	    for (i=1; i < A->n; i++){
	      qx = px + A->dx[i];
	      qy = py + A->dy[i];
	      if(gft::Image32::IsValidPixel(bin, qx,qy)){
		q = qx + qy*bin->ncols;
		if ((bin->data[q] != 0)&&(label->data[q] == 0)){
		  label->data[q] = label->data[p];
		  FIFO[last] = q;
		  last++;
		}
	      }
	    }
	  }
	  l++;
	  first=last=0;
	}
      }
      gft::FreeIntArray(&FIFO);
      return(label);
    }


    sImage32   *GetObjBorder(sImage32 *bin){
      sImage32 *border = gft::Image32::Create(bin->ncols,bin->nrows);
      int p,q,i;
      Pixel u,v;
      sAdjRel *A = gft::AdjRel::Circular(1.0);
      
      for (u.y=0; u.y < bin->nrows; u.y++) 
	for (u.x=0; u.x < bin->ncols; u.x++) {
	  p = u.x + bin->ncols * u.y;
	  if (bin->data[p]>0){
	    for (i=1; i < A->n; i++){
	      v.x = u.x + A->dx[i];
	      v.y = u.y + A->dy[i];
	      if (gft::Image32::IsValidPixel(bin,v.x,v.y)){
		q = v.x + bin->ncols * v.y;
		if (bin->data[q]==0){
		  border->data[p]=1;
		  break;
		}
	      }
	    }
	  }
	}
      gft::AdjRel::Destroy(&A);
      return(border);
    }


    sImage32 *GetObjBorders(sImage32 *img, sAdjRel *A){
      sImage32 *himg=NULL;
      int p,q,i;
      Pixel u,v;
      
      himg = gft::Image32::Create(img->ncols, img->nrows);
      for (u.y=0; u.y < himg->nrows; u.y++){
	for (u.x=0; u.x < himg->ncols; u.x++){
	  p = u.x + himg->ncols*u.y;
	  if (img->data[p] != 0) {
	    for (i=1; i < A->n; i++){
	      v.x = u.x + A->dx[i];
	      v.y = u.y + A->dy[i];
	      if (gft::Image32::IsValidPixel(himg,v.x,v.y)){
		q = v.x + himg->ncols*v.y;
		if (img->data[p] != img->data[q]){
		  himg->data[p] = img->data[p];
		  break;
		}
	      } else {
		himg->data[p] = img->data[p];
		break;
	      }
	    }
	  }
	}
      }
      return(himg);
    }

    //------------------------------------



    sImage32 *LabelContour(sImage32 *bin){
      sImage32 *label,*border,*pred;
      sStack *LIFO;
      int *Cor = NULL;
      int l,t,p,q,pp,qq,qe,qd,i,j,closed = 0,exist;
      int qe_x,qe_y,qd_x,qd_y;
      sAdjRel *A8,*A4,*L,*R;
      Pixel u,v,uu,vv;
     
      A4 = gft::AdjRel::Neighborhood_4();
      A8 = gft::AdjRel::Neighborhood_8_clockwise();
      R = gft::AdjRel::RightSide8(A8);
      L = gft::AdjRel::LeftSide8(A8);

      /***********/
      /*      
      for (i=0; i < A8->n; i++){
	printf("i:%02d: ",i);
	printf("A8 (%2d,%2d)\n",A8->dx[i],A8->dy[i]);
      }

      for (i=0; i < L->n; i++){
	printf("i:%02d: ",i);
	printf("L  (%2d,%2d)\n",L->dx[i], L->dy[i]);
      }

      for (i=0; i < R->n; i++){
	printf("i:%02d: ",i);
	printf("R  (%2d,%2d)\n",R->dx[i], R->dy[i]);
      }
      */
      /***********/      
      
      label = gft::Image32::Create(bin->ncols, bin->nrows);
      pred  = gft::Image32::Create(bin->ncols, bin->nrows);
      LIFO = gft::Stack::Create(bin->n);
      Cor  = gft::AllocIntArray(bin->n);
      for(p = 0; p < bin->n; p++){
	label->data[p] = 0;
	Cor[p] = WHITE;
      }
      border = GetObjBorders(bin, A4);
      //gft::Image32::Write(border, (char *)"border.pgm");

      for(p = 0; p < bin->n; p++){
	if(border->data[p] == 0) continue;
	if(Cor[p] != WHITE) continue;
	u.x = p % bin->ncols;
	u.y = p / bin->ncols;
	exist = 0;
	for (i=1; i < A8->n; i++){
	  v.x = u.x + A8->dx[i];
	  v.y = u.y + A8->dy[i];
	  if (gft::Image32::IsValidPixel(bin,v.x,v.y)){
	    q = v.x + bin->ncols*v.y;
	    if(border->data[q] != 0){
	      qe_x = u.x + L->dx[i];
	      qe_y = u.y + L->dy[i];
	      qd_x = u.x + R->dx[i];
	      qd_y = u.y + R->dy[i];
	      if (gft::Image32::IsValidPixel(bin,qe_x,qe_y) &&
		  gft::Image32::IsValidPixel(bin,qd_x,qd_y)){
		qe = qe_x + bin->ncols*qe_y;
		qd = qd_x + bin->ncols*qd_y;
		if(bin->data[qe] != bin->data[qd])
		  exist = 1;
	      }
	    }
	  }
	}
	if(!exist) continue;
	
	closed = 0;
	Cor[p] = GRAY;
	pred->data[p] = NIL;
	gft::Stack::Push(LIFO, p);
	while(!closed && !gft::Stack::IsEmpty(LIFO)){
	  pp = gft::Stack::Pop(LIFO);
	  Cor[pp] = BLACK;
		
	  uu.x = pp % bin->ncols;
	  uu.y = pp / bin->ncols;
	  for (j=1; j < A8->n; j++){
	    vv.x = uu.x + A8->dx[j];
	    vv.y = uu.y + A8->dy[j];
	    if (gft::Image32::IsValidPixel(bin,vv.x,vv.y)){
	      qq = vv.x + bin->ncols*vv.y;
	      if(qq == p && pred->data[pp]!=p){
		closed = 1;
		break;
	      }
	      if(border->data[qq] != 0 && Cor[qq]!=BLACK){
		qe_x = uu.x + L->dx[j];
		qe_y = uu.y + L->dy[j];
		qd_x = uu.x + R->dx[j];
		qd_y = uu.y + R->dy[j];
		if(gft::Image32::IsValidPixel(bin,qe_x,qe_y) &&
		   gft::Image32::IsValidPixel(bin,qd_x,qd_y)){
		  qe = qe_x + bin->ncols*qe_y;
		  qd = qd_x + bin->ncols*qd_y;
		  if(bin->data[qe] != bin->data[qd]){
		    pred->data[qq] = pp;
		    if(Cor[qq] == WHITE){
		      gft::Stack::Push(LIFO, qq);
		      Cor[qq] = GRAY;
		    }
		  }
		}
	      }
	      
	    }
	  }
	}
	l = 1;
	gft::Stack::Clear(LIFO);
	t = pp;
	while(t != NIL){
	  label->data[t] = l;
	  t = pred->data[t];
	  l++;
	}
      }

      gft::Image32::Destroy(&border);
      gft::Image32::Destroy(&pred);
      gft::Stack::Destroy(&LIFO);
      gft::FreeIntArray(&Cor);
      gft::AdjRel::Destroy(&A4);
      gft::AdjRel::Destroy(&A8);
      gft::AdjRel::Destroy(&L);
      gft::AdjRel::Destroy(&R);
      return label;
    }






    sImage32 *LabelContour(sImage32 *bin, sImage32 *contourid){
      sImage32 *label,*border,*pred;
      sStack *LIFO;
      int *Cor = NULL;
      int l,t,p,q,pp,qq,qe,qd,i,j,closed = 0,exist,id=1;
      int qe_x,qe_y,qd_x,qd_y;
      sAdjRel *A8,*A4,*L,*R;
      Pixel u,v,uu,vv;
     
      A4 = gft::AdjRel::Neighborhood_4();
      A8 = gft::AdjRel::Neighborhood_8_clockwise();
      R = gft::AdjRel::RightSide8(A8);
      L = gft::AdjRel::LeftSide8(A8);

      label = gft::Image32::Create(bin->ncols, bin->nrows);
      pred  = gft::Image32::Create(bin->ncols, bin->nrows);
      LIFO = gft::Stack::Create(bin->n);
      Cor  = gft::AllocIntArray(bin->n);
      for(p = 0; p < bin->n; p++){
	label->data[p] = 0;
	Cor[p] = WHITE;
      }
      border = GetObjBorders(bin, A4);
      gft::Image32::Write(border, (char *)"border.pgm");

      for(p = 0; p < bin->n; p++){
	if(border->data[p] == 0) continue;
	if(Cor[p] != WHITE) continue;
	u.x = p % bin->ncols;
	u.y = p / bin->ncols;
	exist = 0;
	for (i=1; i < A8->n; i++){
	  v.x = u.x + A8->dx[i];
	  v.y = u.y + A8->dy[i];
	  if (gft::Image32::IsValidPixel(bin,v.x,v.y)){
	    q = v.x + bin->ncols*v.y;
	    if(border->data[q] != 0){
	      qe_x = u.x + L->dx[i];
	      qe_y = u.y + L->dy[i];
	      qd_x = u.x + R->dx[i];
	      qd_y = u.y + R->dy[i];
	      if (gft::Image32::IsValidPixel(bin,qe_x,qe_y) &&
		  gft::Image32::IsValidPixel(bin,qd_x,qd_y)){
		qe = qe_x + bin->ncols*qe_y;
		qd = qd_x + bin->ncols*qd_y;
		if(bin->data[qe] != bin->data[qd])
		  exist = 1;
	      }
	    }
	  }
	}
	if(!exist) continue;
	
	closed = 0;
	Cor[p] = GRAY;
	pred->data[p] = NIL;
	gft::Stack::Push(LIFO, p);
	while(!closed && !gft::Stack::IsEmpty(LIFO)){
	  pp = gft::Stack::Pop(LIFO);
	  Cor[pp] = BLACK;
		
	  uu.x = pp % bin->ncols;
	  uu.y = pp / bin->ncols;
	  for (j=1; j < A8->n; j++){
	    vv.x = uu.x + A8->dx[j];
	    vv.y = uu.y + A8->dy[j];
	    if (gft::Image32::IsValidPixel(bin,vv.x,vv.y)){
	      qq = vv.x + bin->ncols*vv.y;
	      if(qq == p && pred->data[pp]!=p){
		closed = 1;
		break;
	      }
	      if(border->data[qq] != 0 && Cor[qq]!=BLACK){
		qe_x = uu.x + L->dx[j];
		qe_y = uu.y + L->dy[j];
		qd_x = uu.x + R->dx[j];
		qd_y = uu.y + R->dy[j];
		if(gft::Image32::IsValidPixel(bin,qe_x,qe_y) &&
		   gft::Image32::IsValidPixel(bin,qd_x,qd_y)){
		  qe = qe_x + bin->ncols*qe_y;
		  qd = qd_x + bin->ncols*qd_y;
		  if(bin->data[qe] != bin->data[qd]){
		    pred->data[qq] = pp;
		    if(Cor[qq] == WHITE){
		      gft::Stack::Push(LIFO, qq);
		      Cor[qq] = GRAY;
		    }
		  }
		}
	      }
	      
	    }
	  }
	}
	l = 1;
	gft::Stack::Clear(LIFO);
	t = pp;
	while(t != NIL){
	  label->data[t] = l;
	  contourid->data[t] = id;
	  t = pred->data[t];
	  l++;
	}
	id++;
      }

      gft::Image32::Destroy(&border);
      gft::Image32::Destroy(&pred);
      gft::Stack::Destroy(&LIFO);
      gft::FreeIntArray(&Cor);
      gft::AdjRel::Destroy(&A4);
      gft::AdjRel::Destroy(&A8);
      gft::AdjRel::Destroy(&L);
      gft::AdjRel::Destroy(&R);
      return label;
    }

    
    
    /*    
    Image32 *DistTrans(Image32 *bin, gft::AdjRel::AdjRel *A, char side){
      Image32 *Dx=NULL,*Dy=NULL,*cost;
      Queue *Q=NULL;
      int i,p,q,n,sz;
      gft::Pixel u,v;
      int *sq=NULL,tmp=INT_MAX,dx,dy;
      AdjPxl *N;
      
      n  = MAX(bin->ncols,bin->nrows);
      sq = AllocIntArray(n);
      for (i=0; i < n; i++) 
	sq[i]=i*i;
      
      sz = FrameSize(A);  
      fbin = AddFrame(bin,sz,0);
      fcont = ObjectBorder(fbin);
      fcost = AddFrame(bin,sz,INT_MIN);
      Dx = CreateImage(fcost->ncols,fcost->nrows);
      Dy = CreateImage(fcost->ncols,fcost->nrows);  
      N  = AdjPixels(fcost,A);
      n  = fcost->ncols*fcost->nrows;
      Q = CreateQueue(2*sz*(sz+bin->ncols+bin->nrows),n);
      
      switch (side) {
      case INTERIOR:
	for(p = 0; p < n; p++){
	  if (fbin->val[p] != 0){
	    if (fcont->val[p]>0){
	      fcost->val[p]=0;    
	      InsertQueue(Q,fcost->val[p]%Q->C.nbuckets,p);
	    }
	    else
	      fcost->val[p] = INT_MAX;	  
	  }
	  else{
	    if (fcost->val[p]!=INT_MIN)
	      fcost->val[p] = 0;
	  }
	}
	break;
      case EXTERIOR:
	for(p = 0; p < n; p++){
	  if (fbin->val[p] == 0){
	    if (fcost->val[p]!=INT_MIN)
	      fcost->val[p] = INT_MAX;	  
	  }
	  else{
	    if (fcont->val[p]>0){
	      fcost->val[p]=0;    
	      InsertQueue(Q,fcost->val[p]%Q->C.nbuckets,p);
	    }
	    else
	      fcost->val[p] = 0;
	  }
	}
	break;
      case BOTH:
      default:    
	for(p = 0; p < n; p++){
	  if (fcont->val[p] > 0){
	    fcost->val[p]=0;    
	    InsertQueue(Q,fcost->val[p]%Q->C.nbuckets,p);
	  }
	  else{ 
	    if (fcost->val[p]!=INT_MIN)
	      fcost->val[p]=INT_MAX;    
	  }
	}
      }
      
      DestroyImage(&fcont);
      DestroyImage(&fbin);
      
      while(!EmptyQueue(Q)) {
	p=RemoveQueue(Q);
	for (i=1; i < N->n; i++){
	  q = p + N->dp[i];
	  if (fcost->val[p] < fcost->val[q]){
	    u.x = p%fcost->ncols;
	    u.y = p/fcost->ncols;
	    v.x = u.x + A->dx[i];
	    v.y = u.y + A->dy[i];
	    dx  = Dx->val[p] + abs(v.x-u.x);
	    dy  = Dy->val[p] + abs(v.y-u.y);
	    tmp = sq[dx] + sq[dy];
	    if (tmp < fcost->val[q]){
	      if (fcost->val[q] == INT_MAX)
		InsertQueue(Q,tmp%Q->C.nbuckets,q);
	      else
		UpdateQueue(Q,q,fcost->val[q]%Q->C.nbuckets,tmp%Q->C.nbuckets);
	      fcost->val[q]  = tmp;
	      Dx->val[q] = dx;
	      Dy->val[q] = dy;
	    }
	  }
	}
      }
      
      DestroyQueue(&Q);
      DestroyAdjPxl(&N);
      cost = RemFrame(fcost,sz);
      
      free(sq);
      DestroyImage(&Dx);
      DestroyImage(&Dy);
      DestroyImage(&fcost);
      
      return(cost);
    }

    Image *SignedDistTrans(Image *bin, AdjRel *A, char side){
      Image *Dx=NULL,*Dy=NULL,*fbin,*fcont,*fcost,*cost;
      Queue *Q=NULL;
      int i,p,q,n,sz;
      Pixel u,v;
      int *sq=NULL,tmp=INT_MAX,dx,dy;
      AdjPxl *N;
      
      n  = MAX(bin->ncols,bin->nrows);
      sq = AllocIntArray(n);
      for (i=0; i < n; i++) 
	sq[i]=i*i;
      
      sz = FrameSize(A);  
      fbin = AddFrame(bin,sz,0);
      fcont = ObjectBorder(fbin);
      fcost = AddFrame(bin,sz,INT_MIN);
      Dx = CreateImage(fcost->ncols,fcost->nrows);
      Dy = CreateImage(fcost->ncols,fcost->nrows);  
      N  = AdjPixels(fcost,A);
      n  = fcost->ncols*fcost->nrows;
      Q = CreateQueue(2*sz*(sz+bin->ncols+bin->nrows),n);
      
      switch (side) {
      case INTERIOR:
	for(p = 0; p < n; p++){
	  if (fbin->val[p] != 0){
	    if (fcont->val[p]>0){
	      fcost->val[p]=0;    
	      InsertQueue(Q,fcost->val[p]%Q->C.nbuckets,p);
	    }
	    else
	      fcost->val[p] = INT_MAX;	  
	  }
	  else{
	    if (fcost->val[p]!=INT_MIN)
	      fcost->val[p] = 0;
	  }
	}
	break;
      case EXTERIOR:
	for(p = 0; p < n; p++){
	  if (fbin->val[p] == 0){
	    if (fcost->val[p]!=INT_MIN)
	      fcost->val[p] = INT_MAX;	  
	  }
	  else{
	    if (fcont->val[p]>0){
	      fcost->val[p]=0;    
	      InsertQueue(Q,fcost->val[p]%Q->C.nbuckets,p);
	    }
	    else
	      fcost->val[p] = 0;
	  }
	}
	break;
      case BOTH:
      default:    
	for(p = 0; p < n; p++){
	  if (fcont->val[p] > 0){
	    fcost->val[p]=0;    
	    InsertQueue(Q,fcost->val[p]%Q->C.nbuckets,p);
	  }
	  else{ 
	    if (fcost->val[p]!=INT_MIN)
	      fcost->val[p]=INT_MAX;    
	  }
	}
      }
      
      DestroyImage(&fcont);
      DestroyImage(&fbin);
      
      while(!EmptyQueue(Q)) {
	p=RemoveQueue(Q);
	for (i=1; i < N->n; i++){
	  q = p + N->dp[i];
	  if (fcost->val[p] < fcost->val[q]){
	    u.x = p%fcost->ncols;
	    u.y = p/fcost->ncols;
	    v.x = u.x + A->dx[i];
	    v.y = u.y + A->dy[i];
	    dx  = Dx->val[p] + abs(v.x-u.x);
	    dy  = Dy->val[p] + abs(v.y-u.y);
	    tmp = sq[dx] + sq[dy];
	    if (tmp < fcost->val[q]){
	      if (fcost->val[q] == INT_MAX)
		InsertQueue(Q,tmp%Q->C.nbuckets,q);
	      else
		UpdateQueue(Q,q,fcost->val[q]%Q->C.nbuckets,tmp%Q->C.nbuckets);
	      fcost->val[q]  = tmp;
	      Dx->val[q] = dx;
	      Dy->val[q] = dy;
	    }
	  }
	}
      }
      
      DestroyQueue(&Q);
      DestroyAdjPxl(&N);
      cost = RemFrame(fcost,sz);
      // sign image
      n  = cost->ncols*cost->nrows;  
      
      if (side != INTERIOR)
	for (i=0; i<n; i++) {
	  if (bin->val[i] == 0) {
	    cost->val[i] = -cost->val[i];
	  }
	}
      free(sq);
      DestroyImage(&Dx);
      DestroyImage(&Dy);
      DestroyImage(&fcost);
      
      return(cost);
    }
    */


    sImage32 *Mask2EDT(sImage32 *bin, sAdjRel *A,
		       char side, int limit, char sign){
      sImage32 *Dx=NULL,*Dy=NULL,*cost,*cont;
      sPQueue32 *Q=NULL;
      int i,p,q,n;
      Pixel u,v;
      int *sq=NULL,tmp=INT_MAX,dx,dy;
      sAdjRel *A4 = gft::AdjRel::Circular(1.0);
      
      n  = MAX(bin->ncols, bin->nrows);
      sq = gft::AllocIntArray(n);
      for (i=0; i < n; i++) 
	sq[i]=i*i;

      cost = Create(bin->ncols, bin->nrows);
      cont = GetObjBorders(bin, A4); 
      Dx = Create(cost->ncols, cost->nrows);
      Dy = Create(cost->ncols, cost->nrows);

      n  = cost->ncols*cost->nrows;

      Q = gft::PQueue32::Create(2*(bin->ncols+bin->nrows), n, cost->data);
      
      switch (side) {
      case INTERIOR:
	for(p = 0; p < n; p++){
	  if (bin->data[p] != 0){
	    if (cont->data[p] > 0){
	      cost->data[p] = 0;
	      gft::PQueue32::InsertElem(&Q, p);
	    } else
	      cost->data[p] = INT_MAX;	  
	  }else{
	    if (cost->data[p] != INT_MIN)
	      cost->data[p] = 0;
	  }
	}
	break;
      case EXTERIOR:
	for(p = 0; p < n; p++){
	  if (bin->data[p] == 0){
	    if (cost->data[p]!=INT_MIN)
	      cost->data[p] = INT_MAX;	  
	  }else{
	    if (cont->data[p]>0){
	      cost->data[p]=0;    
	      gft::PQueue32::InsertElem(&Q, p);
	    }else
	      cost->data[p] = 0;
	  }
	}
	break;
      case BOTH:
      default:    
	for(p = 0; p < n; p++){
	  if (cont->data[p] > 0){
	    cost->data[p]=0;    
	    gft::PQueue32::InsertElem(&Q, p);
	  }else{ 
	    if (cost->data[p]!=INT_MIN)
	      cost->data[p]=INT_MAX;    
	  }
	}
      }
      Destroy(&cont);
      
      while(!gft::PQueue32::IsEmpty(Q)) {
	p = gft::PQueue32::RemoveMinFIFO(Q);
	u.x = p % cost->ncols;
	u.y = p / cost->ncols;
	
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];

	  if(!IsValidPixel(cost, v.x, v.y))
	    continue;
	  q = v.x + v.y * cost->ncols;

	  if (cost->data[p] < cost->data[q]){
	    dx  = Dx->data[p] + abs(A->dx[i]);
	    dy  = Dy->data[p] + abs(A->dy[i]);
	    tmp = sq[dx] + sq[dy];
	    if (tmp < cost->data[q] && tmp <= limit){
	      if (cost->data[q] == INT_MAX){
		cost->data[q]  = tmp;
		gft::PQueue32::InsertElem(&Q, q);
	      }
	      else
		gft::PQueue32::UpdateElem(&Q, q, tmp);
	      Dx->data[q] = dx;
	      Dy->data[q] = dy;
	    }
	  }
	}
      }
      gft::PQueue32::Destroy(&Q);
      gft::AdjRel::Destroy(&A4);
      free(sq);
      Destroy(&Dx);
      Destroy(&Dy);
      
      // Eliminate infinite values */
      n = cost->ncols * cost->nrows;
      for (i=0; i<n; i++) {
	if (cost->data[i]==INT_MAX)
	  cost->data[i] = limit;
      }
      
      // sign scene
      if (sign != 0){
	n  = cost->ncols * cost->nrows;
	if (side != INTERIOR)
	  for (i=0; i<n; i++) {
	    if (bin->data[i] == 0) {
	      cost->data[i] = -cost->data[i];
	    }
	  }
      }
      return(cost);
    }
    

    
    void Mask2EDT(sImage32 *bin, sAdjRel *A,
		  char side, int limit, char sign,
		  sImage32 *cost, sImage32 *root){
      sImage32 *Dx=NULL,*Dy=NULL,*cont;
      sPQueue32 *Q=NULL;
      int i,p,q,n;
      Pixel u,v;
      int *sq=NULL,tmp=INT_MAX,dx,dy;
      sAdjRel *A4 = gft::AdjRel::Circular(1.0);
      
      n  = MAX(bin->ncols, bin->nrows);
      sq = gft::AllocIntArray(n);
      for (i=0; i < n; i++) 
	sq[i]=i*i;

      cont = GetObjBorders(bin, A4); 
      Dx = Create(cost->ncols, cost->nrows);
      Dy = Create(cost->ncols, cost->nrows);

      n  = cost->ncols*cost->nrows;

      Q = gft::PQueue32::Create(2*(bin->ncols+bin->nrows),
				n, cost->data);
      
      switch (side) {
      case INTERIOR:
	for(p = 0; p < n; p++){
	  root->data[p] = p;
	  if (bin->data[p] != 0){
	    if (cont->data[p] > 0){
	      cost->data[p] = 0;
	      gft::PQueue32::InsertElem(&Q, p);
	    } else
	      cost->data[p] = INT_MAX;	  
	  }else{
	    if (cost->data[p] != INT_MIN)
	      cost->data[p] = 0;
	  }
	}
	break;
      case EXTERIOR:
	for(p = 0; p < n; p++){
	  root->data[p] = p;
	  if (bin->data[p] == 0){
	    if (cost->data[p]!=INT_MIN)
	      cost->data[p] = INT_MAX;	  
	  }else{
	    if (cont->data[p]>0){
	      cost->data[p]=0;    
	      gft::PQueue32::InsertElem(&Q, p);
	    }else
	      cost->data[p] = 0;
	  }
	}
	break;
      case BOTH:
      default:    
	for(p = 0; p < n; p++){
	  root->data[p] = p;
	  if (cont->data[p] > 0){
	    cost->data[p]=0;    
	    gft::PQueue32::InsertElem(&Q, p);
	  }else{ 
	    if (cost->data[p]!=INT_MIN)
	      cost->data[p]=INT_MAX;    
	  }
	}
      }
      Destroy(&cont);
      
      while(!gft::PQueue32::IsEmpty(Q)) {
	p = gft::PQueue32::RemoveMinFIFO(Q);
	u.x = p % cost->ncols;
	u.y = p / cost->ncols;
	
	for (i=1; i < A->n; i++){
	  v.x = u.x + A->dx[i];
	  v.y = u.y + A->dy[i];

	  if(!IsValidPixel(cost, v.x, v.y))
	    continue;
	  q = v.x + v.y * cost->ncols;

	  if (cost->data[p] < cost->data[q]){
	    dx  = Dx->data[p] + abs(A->dx[i]);
	    dy  = Dy->data[p] + abs(A->dy[i]);
	    tmp = sq[dx] + sq[dy];
	    if (tmp < cost->data[q] && tmp <= limit){
	      if (cost->data[q] == INT_MAX){
		cost->data[q]  = tmp;
		gft::PQueue32::InsertElem(&Q, q);
	      }
	      else
		gft::PQueue32::UpdateElem(&Q, q, tmp);
	      root->data[q] = root->data[p];
	      Dx->data[q] = dx;
	      Dy->data[q] = dy;
	    }
	  }
	}
      }
      gft::PQueue32::Destroy(&Q);
      gft::AdjRel::Destroy(&A4);
      free(sq);
      Destroy(&Dx);
      Destroy(&Dy);
      
      // Eliminate infinite values 
      n = cost->ncols * cost->nrows;
      for (i=0; i<n; i++) {
	if (cost->data[i]==INT_MAX)
	  cost->data[i] = limit;
      }
      
      // sign scene
      if (sign != 0){
	n  = cost->ncols * cost->nrows;
	if (side != INTERIOR)
	  for (i=0; i<n; i++) {
	    if (bin->data[i] == 0) {
	      cost->data[i] = -cost->data[i];
	    }
	  }
      }
    }
    

    sImage32 *Multiscaleskeletons(sImage32 *bin){
      sAdjRel *A8 = NULL, *A4 = NULL;
      sImage32 *cost, *root, *label, *labelc, *D, *contourid;
      int *N = NULL;
      int p, q, i, Lmax, delta_pq, d;
      Pixel u,v;
      A4 = gft::AdjRel::Neighborhood_4();
      A8 = gft::AdjRel::Neighborhood_8();
      labelc = gft::Image32::LabelBinComp(bin, A8);
      cost = gft::Image32::Create(bin->ncols, bin->nrows);
      root = gft::Image32::Create(bin->ncols, bin->nrows);
      Mask2EDT(bin, A8, INTERIOR, INT_MAX, 0, cost, root);
      gft::Image32::Write(cost, (char *)"edt.pgm");
      //label = LabelContour(bin);
      contourid = gft::Image32::Create(bin->ncols, bin->nrows);
      label = LabelContour(bin, contourid);
      gft::Image32::Write(label, (char *)"contour.pgm");
      for(p = 0; p < bin->n; p++){
	label->data[p] = label->data[root->data[p]];
	contourid->data[p] = contourid->data[root->data[p]];
      }
      gft::Image32::Write(label, (char *)"label.pgm");

      Lmax = gft::Image32::GetMaxVal(labelc);
      N = gft::AllocIntArray(Lmax+1);
      for(p = 0; p < bin->n; p++)
	if(label->data[p] > N[labelc->data[p]])
	  N[labelc->data[p]] = label->data[p];

      D = gft::Image32::Create(bin->ncols, bin->nrows);
      for(p = 0; p < bin->n; p++){
	//if(bin->data[p] == 0) continue;
	if(label->data[p] == 0) continue;
	u.x = p % bin->ncols;
	u.y = p / bin->ncols;
	for (i=1; i < A4->n; i++){
	  v.x = u.x + A4->dx[i];
	  v.y = u.y + A4->dy[i];
	  if(!gft::Image32::IsValidPixel(bin, v.x, v.y))
	    continue;
	  q = v.x + v.y * bin->ncols;
	  //if(bin->data[q] == 0) continue;
	  if(label->data[q] == 0) continue;
	  delta_pq = label->data[q] - label->data[p];
	  d = MIN(delta_pq, N[labelc->data[p]] - delta_pq);
	  //-------------------
	  if(d > 0 && contourid->data[p] != contourid->data[q])
	    d = N[labelc->data[p]];
	  //-------------------
	  if(d > D->data[p])
	    D->data[p] = d;
	}
      }
      gft::FreeIntArray(&N);
      gft::Image32::Destroy(&labelc);
      gft::Image32::Destroy(&label);
      gft::Image32::Destroy(&contourid);
      gft::Image32::Destroy(&cost);
      gft::Image32::Destroy(&root);
      gft::AdjRel::Destroy(&A4);
      gft::AdjRel::Destroy(&A8);
      return D;
    }    



    float PerimeterLength(sImage32 *bin){
      sImage32 *lc;
      int *C;
      int p,q,l,px,py,qx,qy,Lmax;
      float length = 0.0;
      lc = LabelContour(bin);
      Lmax = gft::Image32::GetMaxVal(lc);

      C = gft::AllocIntArray(Lmax+1);
      
      for(p = 0; p < bin->n; p++){
	if(lc->data[p] != 0)
	  C[lc->data[p]] = p;
      }
      C[0] = C[Lmax];
      for(l = 1; l <= Lmax; l++){
	p = C[l-1];
	q = C[l];
	px = p%bin->ncols;
	py = p/bin->ncols;
	qx = q%bin->ncols;
	qy = q/bin->ncols;
	if(abs(px-qx) == 1 && abs(py-qy) == 1)
	  length += 1.4;
	else
	  length += 1.0;
      }
     
      gft::FreeIntArray(&C);
      gft::Image32::Destroy(&lc);
      return length;
    }
    
    
  } /*end Image32 namespace*/
} /*end gft namespace*/

