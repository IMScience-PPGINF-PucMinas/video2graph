
#ifndef _GFT_BMAP_H_
#define _GFT_BMAP_H_

#include "gft_common.h"


namespace gft{
  /**
   * \brief Common definitions and functions to manipulate a vector of booleans.
   */
  namespace BMap{

    /**
     * \brief Vector of booleans. 
     *
     * It uses one bit per boolean (i.e., size = ceil (n / 8)).
     */
    struct sBMap {
      char *data;
      int N, VN;
    };
    
    sBMap *Create(int n);
    void   Destroy(sBMap **b);
    void   Copy(sBMap *dest, sBMap *src);
    void   Fill(sBMap *b, int value);

    inline int   Get(sBMap *b, int p){
      return ((b->data[p>>3]&(1<<(p&0x07)))!=0);
    }

    inline void  Set(sBMap *b, int p, int value){
      if(value) b->data[p>>3]|=(1<<(p&0x07));
      else      b->data[p>>3]&=((~0)^(1<<(p&0x07)));
    }

    inline void Set0(sBMap *b, int p){
      b->data[p>>3]&=((~0)^(1<<(p&0x07)));
    }

    inline void Set1(sBMap *b, int p){
      b->data[p>>3]|=(1<<(p&0x07));
    }

    inline void  Toggle(sBMap *b, int p){
      b->data[p>>3]^=(1<<(p&0x07));
    }


  } //end BMap namespace

  typedef BMap::sBMap sBMap;

} //end gft namespace

#endif

