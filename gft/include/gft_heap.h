#ifndef _GFT_HEAP_H_
#define _GFT_HEAP_H_

#include "gft_common.h"
#include "gft_gpqueue_by_Falcao.h"


namespace gft{
  namespace Heap{

    struct sHeap {
      float *cost;
      char *color;
      int *pixel;
      int *pos;
      int last;
      int n;
    };

    /* Auxiliary Functions */

    //#define HEAP_DAD(i) ((i - 1) / 2)
    //#define HEAP_LEFTSON(i) (2 * i + 1)
    //#define HEAP_RIGHTSON(i) (2 * i + 2)

#define HEAP_DAD(i) (i/2)
#define HEAP_LEFTSON(i) (2 * i)
#define HEAP_RIGHTSON(i) (2 * i + 1)
    
    char IsFull(sHeap *H);
    char IsEmpty(sHeap *H);
    sHeap *Create(int n, float *cost);
    void Destroy(sHeap **H);

    void Insert_MaxPolicy(sHeap *H, int pixel);
    void Insert_MinPolicy(sHeap *H, int pixel);
    
    void Remove_MaxPolicy(sHeap *H, int *pixel);
    void Remove_MinPolicy(sHeap *H, int *pixel);

    void Get_MaxPolicy(sHeap *H, int *pixel);
    void Get_MinPolicy(sHeap *H, int *pixel);
    
    void Update_MaxPolicy(sHeap *H, int p, float value);
    void Update_MinPolicy(sHeap *H, int p, float value);

    void GoUp_MaxPolicy(sHeap *H, int i);
    void GoUp_MinPolicy(sHeap *H, int i);

    void GoDown_MaxPolicy(sHeap *H, int i);
    void GoDown_MinPolicy(sHeap *H, int i);

    void Reset(sHeap *H);

    void Delete_MaxPolicy(sHeap *H, int pixel);
    void Delete_MinPolicy(sHeap *H, int pixel);

    int Debug(sHeap *H, int p);
    
  } //end Heap namespace

  typedef Heap::sHeap sHeap;

} //end gft namespace

#endif

