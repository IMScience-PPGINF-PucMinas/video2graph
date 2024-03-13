
#ifndef _GFT_SET_H_
#define _GFT_SET_H_

#include "gft_common.h"
#include "gft_bmap.h"

#include "gft_image32.h"
#include "gft_adjrel.h"

namespace gft{
  namespace Set{

    struct sSet {
      int elem;
      struct sSet *next;
    };

    
    sSet *Create();
    sSet *Create(sImage32 *bin,
		 sAdjRel *A);
    sSet *Create(sImage32 *img);
    
    void  Destroy(sSet **S);
    sSet *Clone(sSet *S);

    void Insert(sSet **S, int elem);
    int  Remove(sSet **S);
    void RemoveElem(sSet **S, int elem);
    bool IsInSet(sSet *S, int elem);
    int  MinimumValue(sSet *S);
    int  MaximumValue(sSet *S);
    void Convert2DisjointSets(sSet **S1,
			      sSet **S2);
    int  GetNElems(sSet *S);
    
    /**
     * \brief Merge two sets. 
     *
     * The next field of the last element of set S 
     * points to the first element of set T. 
     * T does not change.
     */
    void Merge(sSet **S, sSet **T);

  } //end Set namespace

  typedef Set::sSet sSet;

} //end gft namespace



namespace gft{
  namespace Image32{

    void DrawSet(sImage32 *img,
		 sSet *S, 
		 int value);
    

  } //end Image32 namespace
} //end gft namespace


#include "gft_cimage.h"

namespace gft{
  namespace CImage{
    
    void DrawSet(sCImage *img,
		 sSet *S, 
		 int color);

    
  } //end CImage namespace
} //end gft namespace



#endif

