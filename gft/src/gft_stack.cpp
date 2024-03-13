
#include "gft_stack.h"

namespace gft{
  namespace Stack{

    sStack *Create(int n) {
      sStack *S;
      S = (sStack *) malloc(sizeof(sStack));
      if(S==NULL) gft::Error((char *)MSG1,
			     (char *)"Stack::Create");
      S->n   = n;
      S->top = -1;
      S->data = gft::AllocIntArray(n);
      return S;
    }

    void    Destroy(sStack **S) {
      sStack *aux = *S;
      if(aux) {
	if(aux->data) gft::FreeIntArray(&aux->data);
	free(aux);
	*S = NULL;
      }
    }

    void Clear(sStack *S){
      S->top = -1;
    }
    
    void    Push(sStack *S, int p) {
      S->data[++(S->top)] = p;
    }

    int     Pop(sStack *S) {
      if(S->top == -1) return -1;
      return(S->data[S->top--]);
    }

    
  } //end Stack namespace
} //end gft namespace

