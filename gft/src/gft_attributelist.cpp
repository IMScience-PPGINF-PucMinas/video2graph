
#include "gft_attributelist.h"

namespace gft{
  namespace AttributeList{

    int GetAttributeValue(sAttributeList *al, char name[], char value[]){
      int i;
      value[0] = '\0';
      for(i = 0; i < al->n; i++){
	if(strcmp(al->list[i].name, name) == 0){
	  strcpy(value, al->list[i].value);
	  return 1;
	}
      }
      return 0;
    }


    sAttributeList *Create(int n){
      sAttributeList *al;
      al = (sAttributeList *)malloc(sizeof(sAttributeList));
      if(al == NULL){
	printf("Error: AttributeList::Create\n");
	exit(1);
      }
      al->list = (sAttribute *)malloc(n*sizeof(sAttribute));
      al->n = n;
      return al;
    }


    
    sAttributeList *Create(char buffer[]){
      sAttributeList *al;
      int i,n,c,cc;
      cc = n = i = 0;
      while(buffer[i] != '\0'){
	if(buffer[i] == '='){
	  buffer[i] = ' ';
	  n++;
	}
	i++;
      }
      al = Create(n);
      for(i = 0; i < n; i++){
	sscanf(buffer+cc, "%s%n", al->list[i].name, &c);
	cc += c;
	sscanf(buffer+cc, "%s%n", al->list[i].value,  &c);
	cc += c;
      }
      return al;
    }
    

    
    void Destroy(sAttributeList **al){
      sAttributeList *tmp;
      if(al != NULL){
	tmp = *al;
	if(tmp != NULL){
	  if(tmp->list != NULL)
	    free(tmp->list);
	  free(tmp);
	  *al = NULL;
	}
      }
    }
    

    void ToString(sAttributeList *al, char sep, char str[]){
      char s[2];
      int i;
      str[0] = '\0';
      s[0] = sep;
      s[1] = '\0';
      for(i = 0; i < al->n; i++){
	strcat(str, al->list[i].name);
	strcat(str, (char *)"=");
	strcat(str, al->list[i].value);
	if(i < al->n - 1)
	  strcat(str, s);
      }
    }
    
    void ReplaceCharacterInString(char str[], char oldc, char newc){
      int i = 0;
      while(str[i] != '\0'){
	if(str[i] == oldc)
	  str[i] = newc;
	i++;
      }
    }


  } //end AttributeList namespace
} //end gft namespace

