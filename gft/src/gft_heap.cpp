#include "gft_heap.h"

namespace gft{
  namespace Heap{


    void GoUp_MaxPolicy(sHeap *H, int i){
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


    void GoUp_MinPolicy(sHeap *H, int i){
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
    void GoDown_MaxPolicy(sHeap *H, int i) {
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

    void GoDown_MaxPolicy(sHeap *H, int i) {
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
    

    void GoDown_MinPolicy(sHeap *H, int i) {
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
    void GoDown_MinPolicy(sHeap *H, int i) {
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
    
    
    char IsFull(sHeap *H) {
      if (H->last == H->n)
	return 1;
      else
	return 0;
    }
    
    char IsEmpty(sHeap *H) {
      if (H->last == 0)
	return 1;
      else
	return 0;
    }

    
    sHeap *Create(int n, float *cost) {
      sHeap *H = NULL;
      int i;
      
      if (cost == NULL) {
	fprintf(stdout,"Cannot create heap without cost map in Heap::Create");
	return NULL;
      }
      
      H = (sHeap *) malloc(sizeof(sHeap));
      if (H != NULL) {
	H->n       = n;
	H->cost    = cost;
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
	H->pixel[n] = -1;
      }
      else
	gft::Error((char *)MSG1,(char *)"Heap::Create");

      return H;
    }

    void Destroy(sHeap **H) {
      sHeap *aux = *H;
      if (aux != NULL) {
	if (aux->pixel != NULL) free(aux->pixel);
	if (aux->color != NULL) free(aux->color);
	if (aux->pos != NULL)   free(aux->pos);
	free(aux);
	*H = NULL;
      }
    }
    
    void Insert_MaxPolicy(sHeap *H, int pixel) {
      H->last++;
      H->pixel[H->last] = pixel;
      H->color[pixel]   = GRAY;
      H->pos[pixel]     = H->last;
      GoUp_MaxPolicy(H, H->last);
    }

    
    void Insert_MinPolicy(sHeap *H, int pixel) {
      H->last++;
      H->pixel[H->last] = pixel;
      H->color[pixel]   = GRAY;
      H->pos[pixel]     = H->last;
      GoUp_MinPolicy(H, H->last);
    }

    
    void Remove_MaxPolicy(sHeap *H, int *pixel) {
      *pixel = H->pixel[1];
      H->pos[*pixel]   = -1;
      H->color[*pixel] = BLACK;
      if(H->last == 1){
	H->pixel[1] = -1; //Linha nao obrigatoria.
	H->last = 0;
      }
      else if(H->last > 1){
	H->pixel[1]      = H->pixel[H->last];
	H->pos[H->pixel[1]] = 1;
	H->pixel[H->last] = -1;
	H->last--;
	GoDown_MaxPolicy(H, 1);
      }
    }


    void Remove_MinPolicy(sHeap *H, int *pixel) {
      *pixel = H->pixel[1];
      H->pos[*pixel]   = -1;
      H->color[*pixel] = BLACK;
      if(H->last == 1){
	H->pixel[1] = -1; //Linha nao obrigatoria.
	H->last = 0;
      }
      else if(H->last > 1){
	H->pixel[1]      = H->pixel[H->last];
	H->pos[H->pixel[1]] = 1;
	H->pixel[H->last] = -1;
	H->last--;
	GoDown_MinPolicy(H, 1);
      }
    }


    void Get_MaxPolicy(sHeap *H, int *pixel){
      *pixel = H->pixel[1];
    }

    
    void Get_MinPolicy(sHeap *H, int *pixel){
      *pixel = H->pixel[1];
    }
    

    void Update_MaxPolicy(sHeap *H, int p, float value){
      H->cost[p] = value;
      
      //if (H->color[p] == BLACK) printf("ferrou\n");
      
      if(H->color[p] == WHITE)
	Insert_MaxPolicy(H, p);
      else
	GoUp_MaxPolicy(H, H->pos[p]);
    }


    void Update_MinPolicy(sHeap *H, int p, float value){
      H->cost[p] = value;
      
      //if (H->color[p] == BLACK) printf("ferrou\n");
      
      if(H->color[p] == WHITE)
	Insert_MinPolicy(H, p);
      else
	GoUp_MinPolicy(H, H->pos[p]);
    }
    
    
    void Reset(sHeap *H){
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
    char Delete_MaxPolicy(sHeap *H, int pixel) {
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

    void Delete_MaxPolicy(sHeap *H, int pixel) {
      int q;
      //if (IsEmpty(H)) printf("error: delete pixel%d\n", pixel);
      H->color[pixel] = WHITE;
      q = H->pixel[H->last];
      if(pixel != q){
	H->pixel[H->pos[pixel]] = q;
	H->pos[q] = H->pos[pixel];
	H->pixel[H->last] = -1;
	H->last--;
	H->pos[pixel] = -1;
	if(H->cost[pixel] > H->cost[q])
	  GoDown_MaxPolicy(H, H->pos[q]);
	else 
	  GoUp_MaxPolicy(H, H->pos[q]);
      }
      else{
	H->pixel[H->last] = -1;
	H->last--;
	H->pos[pixel] = -1;
      }
    }
    

    void Delete_MinPolicy(sHeap *H, int pixel) {
      int q;
      //if (IsEmpty(H)) printf("error: delete pixel%d\n", pixel);
      H->color[pixel] = WHITE;
      q = H->pixel[H->last];
      if(pixel != q){
	H->pixel[H->pos[pixel]] = q;
	H->pos[q] = H->pos[pixel];
	H->pixel[H->last] = -1;
	H->last--;
	H->pos[pixel] = -1;
	if(H->cost[pixel] < H->cost[q])
	  GoDown_MinPolicy(H, H->pos[q]);
	else 
	  GoUp_MinPolicy(H, H->pos[q]);
      }
      else{
	H->pixel[H->last] = -1;
	H->last--;
	H->pos[pixel] = -1;
      }
    }

    
    /*
    char Delete_MinPolicy(sHeap *H, int pixel) {
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
    

    int Debug(sHeap *H, int p){
      int i, inside = false, index = -1;
      //printf("pos: %d\n", H->pos[p]);
      for(i = 1; i <= H->last; i++){
	if(H->pixel[i] == p){
	  //printf("Inside Queue at index: %d cost: %f\n", i, H->cost[p]);
	  if(inside == true)
	    printf("Dentro varias vezes\n");
	  index = i;
	  inside = true;
	}
      }
      if(!inside && H->pos[p] == -1)
	return 1;
	//printf("Not inside\n");
      else if(inside && H->pos[p] == index)
	return 2;
      else
	return 0;
    }

    
    
  } /*end Heap namespace*/
} /*end gft namespace*/


