
#include "gft_hashtable.h"

namespace gft{
  namespace HashTable{


    sHashTable *Create(int size){
      sHashTable *ht;
      int i;
      if(size<=0)
	gft::Error((char *)"Invalid size value",
		   (char *)"CreateHashTable");
      ht = (sHashTable *)malloc(sizeof(sHashTable));
      if(ht==NULL)
	gft::Error((char *)MSG1,
		   (char *)"CreateHashTable");
      ht->size = size;
      ht->data = (sHashNode **)malloc(size*sizeof(sHashNode *));
      if(ht->data==NULL)
	gft::Error((char *)MSG1,
		   (char *)"CreateHashTable");
      for(i=0; i<size; i++)
	ht->data[i]=NULL;
      
      return ht;
    }


    void Destroy(sHashTable **ht){
      sHashNode *node,*tmp;
      int i;
      
      if(*ht!=NULL){
	for(i=0; i<(*ht)->size; i++){
	  node = (*ht)->data[i];
	  while(node!=NULL){
	    tmp = node;
	    node = node->next;
	    free(tmp);
	  }
	}
	free((*ht)->data);
	free(*ht);
	*ht = NULL;
      }
    }


    // Function to calculate the position in the hash table
    int HashPosition(sHashTable *ht, char *key){
      int h = 0, a = 127;
      
      for (; *key != '\0'; key++)
	h = (a * h + *key) % ht->size;
      
      return h;
    }


    void Insert(sHashTable *ht, char *key, void *value){
      sHashNode *node;
      int i;
      
      i = HashPosition(ht, key);
      node = (sHashNode *)malloc(sizeof(sHashNode));
      node->key = key;
      node->value = value;
      node->next = ht->data[i];
      ht->data[i] = node;
    }


    void *Search(sHashTable *ht, char *key){
      sHashNode *node;
      int i;
      
      i = HashPosition(ht, key);
      for(node=ht->data[i]; node!=NULL; node=node->next)
	if(strcmp(key, node->key)==0)
	  return node->value;
      
      return NULL;
    }


  } //end HashTable namespace
} //end gft namespace





