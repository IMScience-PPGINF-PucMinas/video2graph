
#include "gft_configfile.h"

namespace gft{
  namespace ConfigFile{


    sConfigFile *Create(){
      sConfigFile *cf;
      cf = (sConfigFile *)malloc(sizeof(sConfigFile));
      if(cf == NULL){
	printf("Error: ConfigFile::Create\n");
	exit(1);
      }
      cf->lines = NULL;
      cf->nlines = 0;
      return cf;
    }


    void Destroy(sConfigFile **cf){
      sConfigFile *tmp;
      int i;
      if(cf != NULL){
	tmp = *cf;
	if(tmp != NULL){
	  if(tmp->lines != NULL){
	    for(i = 0; i < tmp->nlines; i++){
	      gft::AttributeList::Destroy(&tmp->lines[i]);
	    }
	    free(tmp->lines);
	  }
	  free(tmp);
	  *cf = NULL;
	}
      }
    }
    
    
    sConfigFile *Read(char filename[]){
      sConfigFile *cf;
      char buffer[512];
      int r,nlines,i;
      FILE *fp;
      cf = Create();
      fp = fopen(filename, "r");
      if(fp == NULL){
	printf("Error: ConfigFile::Read\n");
	exit(1);
      }
      nlines = 0;
      while(1){
	r = fscanf(fp," %[^\n]", buffer);
	if(r == EOF)
	  break;
	nlines++;
      }
      rewind(fp);
      cf->nlines = nlines;
      cf->lines = (sAttributeList **)malloc(nlines*sizeof(sAttributeList *));
      if(cf->lines == NULL){
	printf("Error: ConfigFile::Read\n");
	exit(1);
      }
      i = 0;
      while(1){
	r = fscanf(fp," %[^\n]", buffer);
	if(r == EOF)
	  break;
	cf->lines[i] = gft::AttributeList::Create(buffer);
	i++;
      }
      fclose(fp);
      return cf;
    }



  } //end ConfigFile namespace
} //end gft namespace

