
#include "gft_set.h"

namespace gft{
  namespace Set{

    sSet *Create(){
      return NULL;
    }


    sSet *Create(sImage32 *bin,
		 sAdjRel *A){
      sImage32 *border=NULL;
      sSet *S=NULL;
      int p,n = bin->ncols*bin->nrows;
      
      if(A==NULL){
	for(p = 0; p < n; p++)
	  if(bin->data[p]>0)
	    Insert(&S, p);
      }
      else{
	border = gft::Image32::GetBoundaries(bin, A);
	S = Create(border, NULL);
	gft::Image32::Destroy(&border);
      }
      return S;
    }


    //Image Border.
    sSet *Create(sImage32 *img){
      sSet *B=NULL;
      int x,y,p;
      
      for (x=0, y=0; x < img->ncols; x++){
	p = x + img->ncols*y;
	Insert(&B,p);
      }
      for (x=img->ncols-1, y=1; y < img->nrows; y++){
	p = x + img->ncols*y;
	Insert(&B,p);
      }
      for (x=img->ncols-2, y=img->nrows-1; x >= 0; x--){
	p = x + img->ncols*y;
	Insert(&B,p);
      }
      for (x=0, y=img->nrows-2; y > 0; y--){
	p = x + img->ncols*y;
	Insert(&B,p);
      }
      return(B);
    }

    
    void Destroy(sSet **S){
      sSet *p;
      while(*S != NULL){
	p = *S;
	*S = p->next;
	free(p);
      }
    }

    //sSet *Clone(sSet *S){
    //  sSet *tmp = NULL;
    //  sSet *C = NULL;
    //  int p;
    //
    //  tmp = S;
    //  while(tmp!=NULL){
    //    p = tmp->elem;
    //    InsertSet(&C, p);
    //    tmp = tmp->next;
    //  }
    //  return C;
    //}

    // Copies the set in an ordered fashion
    sSet *Clone(sSet *S){
      sSet *tmp = NULL;
      sSet *C = NULL;
      sSet **tail = NULL;
      int p;

      tmp = S;

      if(tmp != NULL){
  	p = tmp->elem;
	C = (sSet *) calloc(1,sizeof(sSet));
	C->elem = p;
	C->next = NULL;
	tail = &(C->next);
	tmp = tmp->next;
      }

      while(tmp!=NULL){
  	p = tmp->elem;
	*tail = (sSet *) calloc(1,sizeof(sSet));
	(*tail)->elem = p;
	(*tail)->next = NULL;
	tail = &((*tail)->next);
	tmp = tmp->next;
      }
      return C;
    }


    void Insert(sSet **S, int elem){
      sSet *p=NULL;
      
      p = (sSet *) calloc(1,sizeof(sSet));
      if(p == NULL) 
	gft::Error((char *)MSG1,
		   (char *)"Set::Insert");
      if(*S == NULL){
	p->elem  = elem;
	p->next  = NULL;
      }else{
	p->elem  = elem;
	p->next  = *S;
      }
      *S = p;
    }


    int Remove(sSet **S){
      sSet *p;
      int elem = NIL;

      if(*S != NULL){
	p    =  *S;
	elem = p->elem;
	*S   = p->next;
	//printf("RemoveSet before free");
	free(p);
	//printf(" RemoveSet after free: elem is %d\n",elem);
	//if(*S != NULL) printf(" *S->elem is %d\n",(*S)->elem);
      }
      return(elem);
    }


    void RemoveElem(sSet **S, int elem){
      sSet *tmp = NULL, *remove;
      
      tmp = *S;
      if ( tmp->elem == elem ) {
	*S = tmp->next;
	free( tmp );
      }
      else {
	while( tmp->next->elem != elem ) tmp = tmp->next;
	remove = tmp->next;
	tmp->next = remove->next;
	free( remove );
      }
    }


    bool IsInSet(sSet *S, int elem){
      bool flag = false;
      
      while(S!=NULL){
	if(S->elem == elem){
	  flag = true;
	  break;
	}
	S = S->next;
      }
      return flag;
    }


    int     MinimumValue(sSet *S){
      sSet *aux;
      int p,min=INT_MAX;

      aux = S;
      while(aux != NULL){
	p = aux->elem;
	if(p<min) min = p;
	aux = aux->next;
      }
      return min;
    }


    int     MaximumValue(sSet *S){
      sSet *aux;
      int p,max=INT_MIN;
      
      aux = S;
      while(aux != NULL){
	p = aux->elem;
	if(p>max) max = p;
	aux = aux->next;
      }
      return max;
    }


    void    Convert2DisjointSets(sSet **S1,
				 sSet **S2){
      sBMap *bmap;
      sSet **S[2];
      sSet *Cur,*Prev;
      int elem,i,max,min,p,n;
      
      if(*S1==NULL || *S2==NULL)
	return;
      
      min = MIN(MinimumValue(*S1),
		MinimumValue(*S2));
      
      max = MAX(MaximumValue(*S1),
		MaximumValue(*S2));

      n = max-min+1;
      bmap = BMap::Create(n);
      S[0] = S1;
      S[1] = S2;

      for(i=0; i<2; i++){
	Prev = Cur = *(S[i]);
	while(Cur != NULL){
	  elem = Cur->elem;
	  p = elem - min;
	  if( BMap::Get(bmap, p) ){
	    if(Prev==Cur){
	      *(S[i]) = Cur->next;
	      free(Cur);
	      Prev = Cur = *(S[i]);
	    }
	    else{
	      Prev->next = Cur->next;
	      free(Cur);
	      Cur = Prev->next;
	    }
	  }
	  else{
	    BMap::Set1(bmap, p);
	    Prev = Cur;
	    Cur = Cur->next;
	  }
	}
      }
      BMap::Destroy(&bmap);
    }


    int  GetNElems(sSet *S){
      sSet *aux;
      int size=0;

      aux = S;
      while(aux != NULL){
	size++;
	aux = aux->next;
      }
      return size;
    }


    void Merge(sSet **S, sSet **T) {
      sSet *aux;
      if( *S == NULL ) *S = *T;
      else{
	for( aux = *S; aux->next != NULL; aux = aux->next );
	aux->next = *T;
      }
    }


  } //end Set namespace
} //end gft namespace





namespace gft{
  namespace Image32{


    void DrawSet(sImage32 *img,
		 sSet *S,
		 int value){
      sSet *aux;
      int p;
      aux = S;
      while(aux != NULL){
	p = aux->elem;
	img->data[p] = value;
	aux = aux->next;
      }
    }


  } //end Image32 namespace
} //end gft namespace




namespace gft{
  namespace CImage{


    void DrawSet(sCImage *cimg,
		 sSet *S,
		 int color){
      sSet *aux;
      int p;
      aux = S;
      while(aux != NULL){
	p = aux->elem;
	cimg->C[0]->data[p] = gft::Color::Channel0(color);
	cimg->C[1]->data[p] = gft::Color::Channel1(color);
	cimg->C[2]->data[p] = gft::Color::Channel2(color);
	aux = aux->next;
      }
    }


  } //end CImage namespace
} //end gft namespace

