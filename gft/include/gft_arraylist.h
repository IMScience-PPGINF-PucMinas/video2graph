// version 00.00.02

#ifndef _GFT_ARRAYLIST_H_
#define _GFT_ARRAYLIST_H_

#include "gft_common.h"

namespace gft{
  namespace ArrayList{

    struct sArrayList {
      void **array;
      int cap; //Capacity.
      int n;   //Number of objects added.
      void (*clean)(void**); //Clean function
    };

    sArrayList *Create(int cap);
    void        Destroy(sArrayList **A);
    
    void        SetCleanFunc(sArrayList *A,
			     void (*clean)(void**));

    void       AddElement(sArrayList *A, 
			  void *elem);
    void      *GetElement(sArrayList *A, 
			  int index);
    void       DelElement(sArrayList *A, 
			  int index);
    void       DelElement(sArrayList *A,
			  void **elem);
    
    void       Resize(sArrayList *A, int n);
    
    //Trims the capacity of this ArrayList instance 
    //to be the list's current size.
    void       Trim2Size(sArrayList *A);
    
  } /*end ArrayList namespace*/

  typedef ArrayList::sArrayList sArrayList;
  
} /*end gft namespace*/


#endif

