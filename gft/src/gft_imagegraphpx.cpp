
#include "gft_imagegraphpx.h"

namespace gft{
  namespace ImageGraphPx{


    sImageGraphPx *ByAbsDiff(sImage32 *img, float r){
      sImageGraphPx *sg=NULL;
      sAdjRel *A=AdjRel::Circular(r);
      sg = Create(img->ncols, 
		  img->nrows, A);
      sg->type = DISSIMILARITY;
      sg->Wmax = 0;
      sg->W = NULL;
      sg->img = img;
      sg->cimg = NULL;
      return sg;
    }


    sImageGraphPx *ByAbsDiff(sCImage *cimg, float r){
      sImageGraphPx *sg=NULL;
      sAdjRel *A=AdjRel::Circular(r);
      sg = Create(cimg->C[0]->ncols,
		  cimg->C[0]->nrows, A);
      sg->type = DISSIMILARITY;
      sg->Wmax = 0;
      sg->W = NULL;
      sg->img = NULL;
      sg->cimg = cimg;
      return sg;
    }

    
    sImageGraphPx *ByWeightImage(sImage32 *W, float r){
      sImageGraphPx *sg=NULL;
      sAdjRel *A=AdjRel::Circular(r);
      sg = Create(W->ncols, 
		  W->nrows, A);
      sg->type = DISSIMILARITY;
      sg->Wmax = 0;
      sg->W = W;
      sg->img = NULL;
      sg->cimg = NULL;
      return sg;
    }

    

    sImageGraphPx *Create(int ncols, int nrows, sAdjRel *A){
      sImageGraphPx *sg;
      sg = (sImageGraphPx *)calloc(1, sizeof(sImageGraphPx));
      if(sg == NULL)
	gft::Error((char *)MSG1,(char *)"ImageGraphPx::Create");
      sg->A = A;
      sg->type = DISSIMILARITY;
      sg->ncols = ncols;
      sg->nrows = nrows;
      sg->Wmax = 0;
      sg->W = NULL;
      sg->img = NULL;
      sg->cimg = NULL;      
      return sg;
    }


    void    Destroy(sImageGraphPx **sg){
      sImageGraphPx *aux;  
      aux = *sg;
      if(aux != NULL){
	AdjRel::Destroy(&(aux->A));
	free(*sg);
	*sg = NULL;
      }
    }
    

    sImageGraphPx   *Clone(sImageGraphPx *sg){
      sImageGraphPx *clone = NULL;
      clone = Create(sg->ncols, 
		     sg->nrows, 
		     AdjRel::Clone(sg->A));
      clone->Wmax = sg->Wmax;
      clone->type = sg->type;
      clone->W = sg->W;
      clone->img = sg->img;
      clone->cimg = sg->cimg;
      return clone;
    }



    /*

    void Orient2Digraph(sImageGraph *sg, 
			sImage32 *img, float per){
      sAdjRel *A;
      int p,q,n,i,val,max = 0;
      int u_x,u_y,v_x,v_y;
      float alpha;
      alpha = per/100.0;
      n = sg->ncols*sg->nrows;
      A = sg->A;
      
      for(p=0; p<n; p++){
	u_x = p%img->ncols;
	u_y = p/img->ncols;
	
	for(i=1; i<A->n; i++){
	  if((sg->n_link[p])[i] == INT_MAX ||
	     (sg->n_link[p])[i] == 0)
	    continue;
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(!IsValidPixel(img,v_x,v_y)) continue;
	  q = v_x + img->ncols*v_y;
	  
	  val = (sg->n_link[p])[i];
	  if(img->data[p] > img->data[q]){
	    val = MAX(ROUND(val*(1.0 + alpha)), 1);
	  }
	  else if(img->data[p] < img->data[q]){
	    val = MAX(ROUND(val*(1.0 - alpha)), 1);
	  }
	  (sg->n_link[p])[i] = val;
	  if(val > max) max = val;
	}
      }
      sg->Wmax = max;
    }
    

    


    void Orient2DigraphInner(sImageGraph *sg, sImage32 *P_sum){
      sAdjRel *A;
      int p,q,n,i;//val,max = 0;
      int u_x,u_y,v_x,v_y;
      int new_value;

      if(sg->type == DISSIMILARITY)
	new_value = 0;
      else
	new_value = INT_MAX;
      
      n = sg->ncols*sg->nrows;
      A = sg->A;
      
      for(p=0; p<n; p++){
	u_x = p%sg->ncols;
	u_y = p/sg->ncols;
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(!IsValidPixel(P_sum,v_x,v_y)) continue;
	  q = v_x + sg->ncols*v_y;
	  
	  if(P_sum->data[q] == p)
	    sg->n_link[p][i] = new_value;  //w(P_sum(q),q)=0
	  
	  //val = (sg->n_link[p])[i]; 
	  //if(val > max && val != INT_MAX) max = val;
	}
      }
      //sg->Wmax = max;
    }

    
    void Orient2DigraphOuter(sImageGraph *sg, sImage32 *P_sum){
      sAdjRel *A;
      int p,q,n,i;//val,max = 0;
      int u_x,u_y,v_x,v_y;
      int new_value;

      if(sg->type == DISSIMILARITY)
	new_value = 0;
      else
	new_value = INT_MAX;
      
      n = sg->ncols*sg->nrows;
      A = sg->A;

      for(p=0; p<n; p++){
	u_x = p%sg->ncols;
	u_y = p/sg->ncols;
	
	for(i=1; i<A->n; i++){
	  v_x = u_x + A->dx[i];
	  v_y = u_y + A->dy[i];
	  if(!Image32::IsValidPixel(P_sum,v_x,v_y)) continue;
	  q = v_x + sg->ncols*v_y;
	  
	  if (P_sum->data[p] == q) {
	    //j = get_edge_index(q, p, sg);
	    sg->n_link[p][i] = new_value;  // w(q,P_sum(q)) = 0
	  }
	  //val = (sg->n_link[p])[i]; 
	  //if(val > max && val != INT_MAX) max = val;
	}
      }
      //sg->Wmax = max;
    }
    


    void ChangeType(sImageGraph *sg, 
		    int type){
      int i,p,n;
      sAdjRel *A;
      
      if(sg->type==type)
	return;
      
      A = sg->A;
      n = sg->ncols*sg->nrows;
      for(p=0; p<n; p++){
	for(i=1; i<A->n; i++){
	  if((sg->n_link[p])[i] == INT_MAX)
	    (sg->n_link[p])[i] = 0;
	  else if((sg->n_link[p])[i] == 0)
	    (sg->n_link[p])[i] = INT_MAX;
	  else
	    (sg->n_link[p])[i] = sg->Wmax - (sg->n_link[p])[i] + 1;
	}
      }
      if(sg->type==DISSIMILARITY)
	sg->type = CAPACITY;
      else
	sg->type = DISSIMILARITY;
    }



    // Increasing transformation.
    void Pow(sImageGraph *sg, 
	     int power, int max){
      sAdjRel *A;
      int p,n,i;
      long double v,m;
      
      if(sg->Wmax<=1) return;
      if(power<=0)    return;
      
      n = sg->ncols*sg->nrows;
      A = sg->A;
      
      for(p=0; p<n; p++){
	for(i=1; i<A->n; i++){
	  if((sg->n_link[p])[i] == INT_MAX ||
	     (sg->n_link[p])[i] == 0) continue;
	  v = (double)((sg->n_link[p])[i]);
	  m = powl((v/sg->Wmax),power);
	  (sg->n_link[p])[i] = MAX(ROUND(m*max), 1);
	}
      }
      sg->Wmax = max;
    }
    

*/

    
  } /*end ImageGraph namespace*/
} /*end gft namespace*/


