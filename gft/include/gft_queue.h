
#ifndef _GFT_QUEUE_H_
#define _GFT_QUEUE_H_

#include "gft_common.h"

namespace gft{
  namespace Queue{

    /**
     * \brief FIFO Queue with circular and growing features.
     */
    struct sQueue {
      int *data;
      int get, put;
      int nbuckets;
      int nadded;   //!< Number of elements added.
    };

    sQueue *Create(int nbuckets);
    void    Destroy(sQueue **Q);
    
    void    Push(sQueue *Q, int p);

    /**
     * @return Returns NIL if empty.
     */
    int         Pop(sQueue *Q);
    
    void        Reset(sQueue *Q);
    bool        IsEmpty(sQueue *Q);
    bool        IsFull(sQueue *Q);

  } //end Queue namespace

  typedef Queue::sQueue sQueue;

} //end gft namespace

#endif


