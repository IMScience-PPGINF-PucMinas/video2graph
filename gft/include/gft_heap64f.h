#ifndef _GFT_HEAP64F_H_
#define _GFT_HEAP64F_H_

#include "gft_common.h"
#include "gft_gpqueue_by_Falcao.h"

#include "gft_heap.h"

namespace gft{
  namespace Heap64f{

    struct sHeap64f {
      double *cost;
      char *color;
      int *pixel;
      int *pos;
      int last;
      int n;
    };

    /* Auxiliary Functions */

    char IsFull(sHeap64f *H);
    char IsEmpty(sHeap64f *H);
    sHeap64f *Create(int n, double *cost);
    void Destroy(sHeap64f **H);

    void Insert_MaxPolicy(sHeap64f *H, int pixel);
    void Insert_MinPolicy(sHeap64f *H, int pixel);
    
    void Remove_MaxPolicy(sHeap64f *H, int *pixel);
    void Remove_MinPolicy(sHeap64f *H, int *pixel);

    void Update_MaxPolicy(sHeap64f *H, int p, double value);
    void Update_MinPolicy(sHeap64f *H, int p, double value);

    void GoUp_MaxPolicy(sHeap64f *H, int i);
    void GoUp_MinPolicy(sHeap64f *H, int i);

    void GoDown_MaxPolicy(sHeap64f *H, int i);
    void GoDown_MinPolicy(sHeap64f *H, int i);

    void Reset(sHeap64f *H);

    void Delete_MaxPolicy(sHeap64f *H, int pixel);
    void Delete_MinPolicy(sHeap64f *H, int pixel);
    

  } //end Heap64f namespace

  typedef Heap64f::sHeap64f sHeap64f;

} //end gft namespace

#endif

