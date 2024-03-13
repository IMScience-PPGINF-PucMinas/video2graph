
#ifndef _GFT_CONFIGFILE_H_
#define _GFT_CONFIGFILE_H_

#include "gft_attributelist.h"

namespace gft{
  namespace ConfigFile{

    struct sConfigFile{
      sAttributeList **lines;
      int nlines;
    };
    
    sConfigFile *Create();
    void Destroy(sConfigFile **cf);
    
    sConfigFile *Read(char filename[]);

  } //end ConfigFile namespace

  typedef ConfigFile::sConfigFile sConfigFile;

} //end gft namespace

 
#endif

