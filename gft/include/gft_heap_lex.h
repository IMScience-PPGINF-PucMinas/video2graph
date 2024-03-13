#ifndef _GFT_HEAP_LEX_H_
#define _GFT_HEAP_LEX_H_

#include "gft_common.h"
#include "gft_heap.h"
#include "gft_gpqueue_by_Falcao.h"

namespace gft{
  namespace Heap_lex{

    struct sHeap_lex {
      float *cost1,*cost2;
      char *color;
      int *pixel;
      int *pos;
      int last;
      int n;
      char removal_policy; /* 0 is MINVALUE and 1 is MAXVALUE */
    };
   

    void SetRemovalPolicy(sHeap_lex *H, char policy);
    char IsFull(sHeap_lex *H);
    char IsEmpty(sHeap_lex *H);
    sHeap_lex *Create(int n, float *cost1, float *cost2);
    void Destroy(sHeap_lex **H);
    char Insert(sHeap_lex *H, int pixel);
    char Remove(sHeap_lex *H, int *pixel);
    void Update(sHeap_lex *H, int p, float value1, float value2);
    void GoUp(sHeap_lex *H, int i);
    void GoDown(sHeap_lex *H, int i);
    void Reset(sHeap_lex *H);

  } //end Heap_lex namespace

  typedef Heap_lex::sHeap_lex sHeap_lex;

} //end gft namespace

#endif



