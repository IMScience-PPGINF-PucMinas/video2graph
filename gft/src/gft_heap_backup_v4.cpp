#include "gft_heap.h"

namespace gft{
  namespace Heap{


    void GoUp_MaxPolicy(Heap *H, int i){
      int j = HEAP_DAD(i);
      int p = H->pixel[i];
      float c = H->cost[p];
      while((i > 1)&&(H->cost[H->pixel[j]] < c)){
	H->pixel[i] = H->pixel[j];
	H->pos[H->pixel[j]] = i;
	i = j;
	j = HEAP_DAD(i);
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }


    void GoUp_MinPolicy(Heap *H, int i){
      int j = HEAP_DAD(i);
      int p = H->pixel[i];
      float c = H->cost[p];
      while((i > 1)&&(H->cost[H->pixel[j]] > c)){
	H->pixel[i] = H->pixel[j];
	H->pos[H->pixel[j]] = i;
	i = j;
	j = HEAP_DAD(i);
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }    

    
    /*    
    void GoDown_MaxPolicy(Heap *H, int i) {
      int j, left = HEAP_LEFTSON(i), right = HEAP_RIGHTSON(i);
      
      j = i;
      if ((left <= H->last) && 
	  (H->cost[H->pixel[left]] > H->cost[H->pixel[i]]))
	j = left;
      if ((right <= H->last) && 
	  (H->cost[H->pixel[right]] > H->cost[H->pixel[j]]))
	j = right;
      
      if(j != i) {
	gft::SwapInt(&H->pixel[j], &H->pixel[i]);
	H->pos[H->pixel[i]] = i;
	H->pos[H->pixel[j]] = j;
	GoDown_MaxPolicy(H, j);
      }
    }
    */

    void GoDown_MaxPolicy(Heap *H, int i) {
      int j, left, right, p = H->pixel[i];
      float c = H->cost[p];
      j = i;
      while(1){
	i = j;
	left = HEAP_LEFTSON(i);
	right = left + 1;
	if(right <= H->last){
	  if(H->cost[H->pixel[right]] > H->cost[H->pixel[left]])
	    j = right;
	  else
	    j = left;
	}
	else if(left <= H->last)
	  j = left;
	else
	  break;
	
	if(H->cost[H->pixel[j]] > c) {
	  H->pixel[i] = H->pixel[j];
	  H->pos[H->pixel[j]] = i;
	}
	else break;
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }
    

    void GoDown_MinPolicy(Heap *H, int i) {
      int j, left, right, p = H->pixel[i];
      float c = H->cost[p];
      j = i;
      while(1){
	i = j;
	left = HEAP_LEFTSON(i);
	right = left + 1;
	if(right <= H->last){
	  if(H->cost[H->pixel[right]] < H->cost[H->pixel[left]])
	    j = right;
	  else
	    j = left;
	}
	else if(left <= H->last)
	  j = left;
	else
	  break;
	
	if(H->cost[H->pixel[j]] < c) {
	  H->pixel[i] = H->pixel[j];
	  H->pos[H->pixel[j]] = i;
	}
	else break;
      }
      H->pixel[i] = p;
      H->pos[p] = i;
    }

    
    /*
    void GoDown_MinPolicy(Heap *H, int i) {
      int j, left = HEAP_LEFTSON(i), right = HEAP_RIGHTSON(i);
      
      j = i;
      if ((left <= H->last) && 
	  (H->cost[H->pixel[left]] < H->cost[H->pixel[i]]))
	j = left;
      if ((right <= H->last) && 
	  (H->cost[H->pixel[right]] < H->cost[H->pixel[j]]))
	j = right;

      if(j != i) {
	gft::SwapInt(&H->pixel[j], &H->pixel[i]);
	H->pos[H->pixel[i]] = i;
	H->pos[H->pixel[j]] = j;
	GoDown_MinPolicy(H, j);
      }
    }
    */
    
    
    char IsFull(Heap *H) {
      if (H->last == H->n)
	return 1;
      else
	return 0;
    }
    
    char IsEmpty(Heap *H) {
      if (H->last == 0)
	return 1;
      else
	return 0;
    }


    Heap *Create_MaxPolicy(int n, float **cost) {
      Heap *H = Create_MinPolicy(n, cost);
      H->cost[n] = FLT_MAX;
      return H;
    }
    
    
    Heap *Create_MinPolicy(int n, float **cost) {
      Heap *H = NULL;
      int i;
      
      if (cost == NULL) {
	fprintf(stdout,"Cannot create heap without cost map in Heap::Create");
	return NULL;
      }
      
      H = (Heap *) malloc(sizeof(Heap));
      if (H != NULL) {
	H->n       = n;
	H->cost    = gft::AllocFloatArray(n+1);
	H->color   = (char *) malloc(sizeof(char) * n);
	H->pixel   = (int *) malloc(sizeof(int) * (n+1));
	H->pos     = (int *) malloc(sizeof(int) * n);
	H->last    = 0;
	if (H->color == NULL || H->pos == NULL || H->pixel == NULL)
	  gft::Error((char *)MSG1,(char *)"Heap::Create");
	for (i = 0; i < H->n; i++) {
	  H->color[i] = WHITE;
	  H->pos[i]   = -1;
	  H->pixel[i] = -1;
	}
	H->pixel[0] = n;
	H->pixel[n] = -1;
	H->cost[n] = -FLT_MAX;
	*cost = H->cost;
      }
      else
	gft::Error((char *)MSG1,(char *)"Heap::Create");

      return H;
    }

    
    void Destroy(Heap **H) {
      Heap *aux = *H;
      if (aux != NULL) {
	if (aux->pixel != NULL) free(aux->pixel);
	if (aux->color != NULL) free(aux->color);
	if (aux->pos != NULL)   free(aux->pos);
	free(aux);
	*H = NULL;
      }
    }
    
    void Insert_MaxPolicy(Heap *H, int pixel) {
      H->last++;
      H->pixel[H->last] = pixel;
      H->color[pixel]   = GRAY;
      H->pos[pixel]     = H->last;
      GoUp_MaxPolicy(H, H->last);
    }

    
    void Insert_MinPolicy(Heap *H, int pixel) {
      H->last++;
      H->pixel[H->last] = pixel;
      H->color[pixel]   = GRAY;
      H->pos[pixel]     = H->last;
      GoUp_MinPolicy(H, H->last);
    }

    
    void Remove_MaxPolicy(Heap *H, int *pixel) {
      *pixel = H->pixel[1];
      H->pos[*pixel]   = -1;
      H->color[*pixel] = BLACK;
      if(H->last == 1){
	H->pixel[1] = -1; //Linha nao obrigatoria.
	H->last = 0;
      }
      else{ //if(H->last > 1){
	H->pixel[1]      = H->pixel[H->last];
	H->pos[H->pixel[1]] = 1;
	H->pixel[H->last] = -1;
	H->last--;
	GoDown_MaxPolicy(H, 1);
      }
    }


    void Remove_MinPolicy(Heap *H, int *pixel) {
      *pixel = H->pixel[1];
      H->pos[*pixel]   = -1;
      H->color[*pixel] = BLACK;
      if(H->last == 1){
	H->pixel[1] = -1; //Linha nao obrigatoria.
	H->last = 0;
      }
      else{ //if(H->last > 1){
	H->pixel[1]      = H->pixel[H->last];
	H->pos[H->pixel[1]] = 1;
	H->pixel[H->last] = -1;
	H->last--;
	GoDown_MinPolicy(H, 1);
      }
    }


    void Update_MaxPolicy(Heap *H, int p, float value){
      H->cost[p] = value;
      
      //if (H->color[p] == BLACK) printf("ferrou\n");
      
      if(H->color[p] == WHITE)
	Insert_MaxPolicy(H, p);
      else
	GoUp_MaxPolicy(H, H->pos[p]);
    }


    void Update_MinPolicy(Heap *H, int p, float value){
      H->cost[p] = value;
      
      //if (H->color[p] == BLACK) printf("ferrou\n");
      
      if(H->color[p] == WHITE)
	Insert_MinPolicy(H, p);
      else
	GoUp_MinPolicy(H, H->pos[p]);
    }
    
    
    void Reset(Heap *H){
      int i;
      
      for (i=0; i < H->n; i++) {
	H->color[i] = WHITE;
	H->pos[i]   = -1;
	H->pixel[i] = -1;
      }
      H->pixel[H->n] = -1;
      H->last = 0;
    }


    /*
    char Delete_MaxPolicy(Heap *H, int pixel) {
      int i, t = 0;
      if (!IsEmpty(H)) {
	Update_MaxPolicy(H,pixel,FLT_MAX);
	Remove_MaxPolicy(H, &i);
	H->color[pixel] = WHITE;
	return 1;
      }
      else
	return 0;
    }
    */

    void Delete_MaxPolicy(Heap *H, int pixel) {
      int q;
      //if (IsEmpty(H)) printf("error: delete pixel%d\n", pixel);
      q = H->pixel[H->last];
      H->pixel[H->pos[pixel]] = q;
      H->pos[q] = H->pos[pixel];
      H->pixel[H->last] = -1;
      H->last--;
      H->pos[pixel] = -1;
      H->color[pixel] = WHITE;
      if(H->cost[pixel] > H->cost[q])
	GoDown_MaxPolicy(H, H->pos[q]);
      else 
	GoUp_MaxPolicy(H, H->pos[q]);
    }
    

    void Delete_MinPolicy(Heap *H, int pixel) {
      int q;
      //if (IsEmpty(H)) printf("error: delete pixel%d\n", pixel);
      q = H->pixel[H->last];
      H->pixel[H->pos[pixel]] = q;
      H->pos[q] = H->pos[pixel];
      H->pixel[H->last] = -1;
      H->last--;
      H->pos[pixel] = -1;
      H->color[pixel] = WHITE;
      if(H->cost[pixel] < H->cost[q])
	GoDown_MinPolicy(H, H->pos[q]);
      else 
	GoUp_MinPolicy(H, H->pos[q]);
    }

    
    /*
    char Delete_MinPolicy(Heap *H, int pixel) {
      int i, t = 0;
      if (!IsEmpty(H)) {
	Update_MinPolicy(H,pixel,-FLT_MAX);
	Remove_MinPolicy(H, &i);
	H->color[pixel] = WHITE;
	return 1;
      }
      else
	return 0;
    }
    */
    

    
  } /*end Heap namespace*/
} /*end gft namespace*/


