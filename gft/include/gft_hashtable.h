#ifndef _GFT_HASHTABLE_H_
#define _GFT_HASHTABLE_H_

#include "gft_common.h"

namespace gft{
  namespace HashTable{

    /**
     * \brief A hashing table node.
     */
    struct sHashNode {
      char *key;               //!< Search key.
      void *value;             //!< Associated value.
      struct sHashNode *next;  //!< Pointer to next node.
    };

    struct sHashTable {
      sHashNode **data;
      int size;
    };

    /**
     * \brief A constructor.
     */ 
    sHashTable *Create(int size);

    /**
     * \brief A destructor.
     */
    void       Destroy(sHashTable **ht);
    
    /**
     * \brief Insert the pair (key, value) in the hash table.
     */
    void Insert(sHashTable *ht, char *key, void *value);

    /**
     * \brief Search for the value associated with a key.
     */
    void *Search(sHashTable *ht, char *key);


  } //end HashTable namespace

  typedef HashTable::sHashTable sHashTable;

} //end gft namespace

#endif

