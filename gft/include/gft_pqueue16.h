
#ifndef _GFT_PQUEUE16_H_
#define _GFT_PQUEUE16_H_

#include "gft_common.h"

namespace gft{
  namespace PQueue16{

    struct sPQNode { 
      int  next;  //!< Next node.
      int  prev;  //!< Previous node.
      char color; //!< WHITE=0, GRAY=1, BLACK=2.
    };
    
    struct sPQDoublyLinkedLists {
      sPQNode *elem;  //!< All possible doubly-linked lists of the circular queue.
      int nelems;    //!< Total number of elements.
      ushort *value; //!< The value of the nodes in the graph.
    };

    struct sPQCircularQueue { 
      int  *first;   //!< List of the first elements of each doubly-linked list.
      int  *last;    //!< List of the last  elements of each doubly-linked list.
      int  nbuckets; //!< Number of buckets in the circular queue.
      ushort minvalue; //!< Minimum value of a node in queue.
      ushort maxvalue; //!< Maximum value of a node in queue.
    };


    /**
     * \brief Priority queue by Dial implemented as proposed by 
     * A.X. Falcao with circular and growing features.
     */
    struct sPQueue16 { 
      sPQCircularQueue C;
      sPQDoublyLinkedLists L;
      int nadded;      //!< Number of elements added.
    };


    sPQueue16  *Create(int nbuckets, int nelems, ushort *value);
    void        Destroy(sPQueue16 **Q);
    sPQueue16  *Grow(sPQueue16 **Q, int nbuckets);
    void        Reset(sPQueue16 *Q);
    inline bool IsEmpty(sPQueue16 *Q){ 
      return (Q->nadded==0); 
    }
    inline bool IsFull(sPQueue16 *Q){ 
      return (Q->nadded==(Q->L).nelems); 
    }
    
    /**
     * Generic version with circular and growing features.
     */
    void   InsertElem(sPQueue16 **Q, int elem);
    /**
     * Generic version with circular and growing features.
     */
    void   RemoveElem(sPQueue16 *Q, int elem);
    /**
     * Generic version with circular and growing features.
     */
    void   UpdateElem(sPQueue16 **Q, int elem, ushort newvalue);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMinFIFO(sPQueue16 *Q);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMinLIFO(sPQueue16 *Q);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMaxFIFO(sPQueue16 *Q);
    /**
     * Generic version with circular and growing features.
     */
    int    RemoveMaxLIFO(sPQueue16 *Q);


    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    void   FastInsertElem(sPQueue16 *Q, int elem);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    void   FastRemoveElem(sPQueue16 *Q, int elem);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    void   FastUpdateElem(sPQueue16 *Q, int elem, ushort newvalue);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMinFIFO(sPQueue16 *Q);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMinLIFO(sPQueue16 *Q);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMaxFIFO(sPQueue16 *Q);
    /**
     * Faster version to be used when values 
     * are in fixed range [0, nbuckets-1] (e.g., watershed transform).
     */
    int    FastRemoveMaxLIFO(sPQueue16 *Q);

    
  } //end PQueue16 namespace

  typedef PQueue16::sPQueue16 sPQueue16;

} //end gft namespace

#endif


