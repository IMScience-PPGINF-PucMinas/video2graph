
#ifndef _GFT_ATTRIBUTELIST_H_
#define _GFT_ATTRIBUTELIST_H_

#include "gft_common.h"

namespace gft{
  namespace AttributeList{
    
    struct sAttribute{
      char value[128];
      char name[128];
    };
    
    struct sAttributeList{
      sAttribute *list;
      int n;
    };
    
    int GetAttributeValue(sAttributeList *al, char name[], char value[]);

    sAttributeList *Create(int n);
    sAttributeList *Create(char buffer[]);
    void Destroy(sAttributeList **al);

    void ToString(sAttributeList *al, char sep, char str[]);
    void ReplaceCharacterInString(char str[], char oldc, char newc);
    
  } //end AttributeList namespace

  typedef AttributeList::sAttributeList sAttributeList;

} //end gft namespace


#endif

