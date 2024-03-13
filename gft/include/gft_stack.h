
#ifndef _GFT_STACK_H
#define _GFT_STACK_H 1

#include "gft_common.h"

namespace gft{
  namespace Stack{

    struct sStack {
      int *data;
      int top;
      int n;
    };
    
    sStack *Create(int n);
    void    Destroy(sStack **S);
    void    Push(sStack *S, int p);
    /**
     * @return Returns NIL if empty.
     */
    int    Pop(sStack *S);

    void Clear(sStack *S);
    
    inline int IsEmpty(sStack *S){ return(S->top == -1); }

  } //end Stack namespace

  typedef Stack::sStack sStack;

} //end gft namespace

#endif
