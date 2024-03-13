#ifndef _GFT_PQUEUE32_H_
#define _GFT_PQUEUE32_H_

#include "gft_common.h"

namespace gft{
  namespace PQueue32{

    struct sPQNode { 
      int  next;  //!< Next node.
      int  prev;  //!< Previous node.
      char color; //!< WHITE=0, GRAY=1, BLACK=2.
    };
    
    struct sPQDoublyLinkedLists {
      sPQNode *elem; //!< All possible doubly-linked lists of the circular queue.
      int nelems;   //!< Total number of elements.
      int *value;   //!< The value of the nodes in the graph.
    }; 

    struct sPQCircularQueue { 
      int  *first;   //!< List of the first elements of each doubly-linked list.
      int  *last;    //!< List of the last  elements of each doubly-linked list.
      int  nbuckets; //!< Number of buckets in the circular queue.
      int  minvalue; //!< Minimum value of a node in queue.
      int  maxvalue; //!< Maximum value of a node in queue.
    };


    /**
     * \brief Priority queue by Dial implemented as proposed by 
     * A.X. Falcao with circular and growing features.
     */
    struct sPQueue32 { 
      sPQCircularQueue C;
      sPQDoublyLinkedLists L;
      int nadded;      //!< Number of elements added.
    };


    sPQueue32  *Create(int nbuckets, int nelems, int *value);
    void        Destroy(sPQueue32 **Q);
    sPQueue32  *Grow(sPQueue32 **Q, int nbuckets);
    void        Reset(sPQueue32 *Q);
    inline bool IsEmpty(sPQueue32 *Q){ 
      return (Q->nadded==0); 
    }
    inline bool IsFull(sPQueue32 *Q){ 
      return (Q->nadded==(Q->L).nelems); 
    }
    

    /**
     * Generic version with circular and growing features.
     */
    void   InsertElem(sPQueue32 **Q, int elem);
    /**
     * Generic version with circular and growing features.
     */
    void   RemoveElem(sPQueue32 *Q, int elem);
    /**
     * Generic version with circular and growing features.
     */
    void   UpdateElem(sPQueue32 **Q, int elem, int newvalue);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMinFIFO(sPQueue32 *Q);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMinLIFO(sPQueue32 *Q);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMaxFIFO(sPQueue32 *Q);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMaxLIFO(sPQueue32 *Q);


    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    void   FastInsertElem(sPQueue32 *Q, int elem);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    void   FastRemoveElem(sPQueue32 *Q, int elem);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    void   FastUpdateElem(sPQueue32 *Q, int elem, int newvalue);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMinFIFO(sPQueue32 *Q);

    int    FastGetMinFIFO(sPQueue32 *Q);

    int    FastGetMinVal(sPQueue32 *Q);
    int    FastGetMaxVal(sPQueue32 *Q);
    
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMinLIFO(sPQueue32 *Q);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMaxFIFO(sPQueue32 *Q);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMaxLIFO(sPQueue32 *Q);

    void   FastInsertElemAsFirst(sPQueue32 *Q, int elem);
    
  } //end PQueue32 namespace

  typedef PQueue32::sPQueue32 sPQueue32;  

} //end gft namespace

#endif


