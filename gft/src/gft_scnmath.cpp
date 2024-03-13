
#include "gft_scnmath.h"

namespace gft{
  namespace Scene32{

    sScene32 *Sub(sScene32 *scn1, sScene32 *scn2){
      sScene32 *diff=NULL;
      v4si *ptr1,*ptr2,*ptr3;
      int p;
      diff = Create(scn1);
      for(p=0; p<scn1->n; p+=4){
	//diff->data[p] = scn1->data[p] - scn2->data[p];
	ptr1 = (v4si *)(scn1->data + p);
	ptr2 = (v4si *)(scn2->data + p);
	ptr3 = (v4si *)(diff->data + p);
	*ptr3 = *ptr1 - *ptr2;
      }
      return diff;
    }


    //inplace version.
    void   Subinplace(sScene32 *scn1, sScene32 *scn2){
      v4si *ptr1,*ptr2;
      int p;
      for(p=0; p<scn1->n; p+=4){
	//scn1->data[p] = scn1->data[p] - scn2->data[p];
	ptr1 = (v4si *)(scn1->data + p);
	ptr2 = (v4si *)(scn2->data + p);
	*ptr1 = *ptr1 - *ptr2;
      }
    } 


    sScene32 *Add(sScene32 *scn1, sScene32 *scn2){
      sScene32 *sum=NULL;
      v4si *ptr1,*ptr2,*ptr3;
      int p;
      sum = Create(scn1);
      for(p=0; p<scn1->n; p+=4){
	//sum->data[p] = scn1->data[p] + scn2->data[p];
	ptr1 = (v4si *)(scn1->data + p);
	ptr2 = (v4si *)(scn2->data + p);
	ptr3 = (v4si *)(sum->data  + p);
	*ptr3 = *ptr1 + *ptr2;
      } 
      return sum;
    }


    //inplace version.
    void   Addinplace(sScene32 *scn1, sScene32 *scn2){
      v4si *ptr1,*ptr2;
      int p;
      for(p=0; p<scn1->n; p+=4){
	//scn1->data[p] = scn1->data[p] + scn2->data[p];
	ptr1 = (v4si *)(scn1->data + p);
	ptr2 = (v4si *)(scn2->data + p);
	*ptr1 = *ptr1 + *ptr2;
      }
    }


    sScene32 *Add(sScene32 *scn, int value){
      sScene32 *sum=NULL;
      v4si *ptr1,*ptr2;
      v4si v;
      int p;
      sum = Create(scn);
      ((int *)(&v))[0] = value;
      ((int *)(&v))[1] = value;
      ((int *)(&v))[2] = value;
      ((int *)(&v))[3] = value;
      for(p=0; p<scn->n; p+=4){
	//sum->data[p] = scn->data[p] + value;
	ptr1 = (v4si *)(scn->data + p);
	ptr2 = (v4si *)(sum->data + p);
	*ptr2 = *ptr1 + v;
      }
      sum->maxval = scn->maxval + value;
      return sum;
    }


    sScene32 *Mult(sScene32 *scn1, sScene32 *scn2){
      sScene32 *mscn=NULL;
      v4si *ptr1,*ptr2;
      int p;
      mscn = Clone(scn2);
      for(p=0; p<scn1->n; p+=4){
	//mscn->data[p]*=scn1->data[p];
	ptr1 = (v4si *)(scn1->data + p);
	ptr2 = (v4si *)(mscn->data + p);
	(*ptr2) *= (*ptr1);
      }
      return(mscn);
    }


    sScene32 *Or(sScene32 *scn1, sScene32 *scn2){
      sScene32 *scnor=NULL;
      int p;
      scnor = Create(scn1);
      for(p = 0; p < scn1->n; p++)
	scnor->data[p] = MAX(scn1->data[p],scn2->data[p]);

      return scnor;
    }


    //inplace version.
    void   Orinplace(sScene32 *scn1, sScene32 *scn2){
      int p;
      for(p=0; p<scn1->n; p++)
	scn1->data[p] = MAX(scn1->data[p],scn2->data[p]);
    }


    sScene32 *And(sScene32 *scn1, sScene32 *scn2){
      sScene32 *scnand=NULL;
      int p;
      scnand = Create(scn1);
      for(p=0; p<scn1->n; p++)
	scnand->data[p] = MIN(scn1->data[p],scn2->data[p]);
      
      return scnand;
    }


    sScene32 *XOr(sScene32 *scn1, sScene32 *scn2){
      sScene32 *scnxor=NULL;
      int p;
      scnxor = Create(scn1);
      for(p = 0; p < scn1->n; p++){
	//scnxor->data[p] = (MAX(scn1->data[p],scn2->data[p])-
	//		     MIN(scn1->data[p],scn2->data[p]));
	scnxor->data[p] = ((scn1->data[p]>scn2->data[p])?
			   (scn1->data[p]-scn2->data[p]):
			   (scn2->data[p]-scn1->data[p]));
      }
      return scnxor;
    }


    sScene32 *Complement(sScene32 *scn){
      sScene32 *cscn=NULL;
      v4si *ptr1,*ptr2;
      v4si v;
      int p,Imax;
      cscn = Create(scn);
      Imax = GetMaximumValue(scn);
      ((int *)(&v))[0] = Imax;
      ((int *)(&v))[1] = Imax;
      ((int *)(&v))[2] = Imax;
      ((int *)(&v))[3] = Imax;
      for(p=0; p < scn->n; p+=4){
	//cscn->data[p] = Imax - scn->data[p];
	ptr1 = (v4si *)(scn->data + p);
	ptr2 = (v4si *)(cscn->data + p);
	*ptr2 = v - *ptr1;
      }
      return(cscn);
    }


    sScene32 *Abs(sScene32 *scn){
      sScene32 *absscn=NULL;
      int p;
      absscn = Create(scn);
      for(p=0; p<scn->n; p++)
	absscn->data[p] = abs(scn->data[p]);

      return(absscn);
    }


    //inplace version.
    void     Negateinplace(sScene32 *scn){
      v4si *ptr;
      int p;
      for(p=0; p<scn->n; p+=4){
	//scn->data[p] = -scn->data[p];
	ptr = (v4si *)(scn->data + p);
	*ptr = -(*ptr);
      }
    } 


  } //end Scene32 namespace



  namespace Scene16{

    sScene16 *Sub(sScene16 *scn1, sScene16 *scn2){
      sScene16 *diff=NULL;
      v8hi *ptr1,*ptr2,*ptr3;
      int p;
      diff = Create(scn1);
      for(p=0; p<scn1->n; p+=8){
	//diff->data[p] = scn1->data[p] - scn2->data[p];
	ptr1 = (v8hi *)(scn1->data + p);
	ptr2 = (v8hi *)(scn2->data + p);
	ptr3 = (v8hi *)(diff->data + p);
	*ptr3 = *ptr1 - *ptr2;
      }
      return diff;
    }


    //inplace version.
    void     Subinplace(sScene16 *scn1, sScene16 *scn2){
      v8hi *ptr1,*ptr2;
      int p;
      for(p=0; p<scn1->n; p+=8){
	//scn1->data[p] = scn1->data[p] - scn2->data[p];
	ptr1 = (v8hi *)(scn1->data + p);
	ptr2 = (v8hi *)(scn2->data + p);
	*ptr1 = *ptr1 - *ptr2;
      }
    }


    sScene16 *Add(sScene16 *scn1, sScene16 *scn2){
      sScene16 *sum=NULL;
      v8hi *ptr1,*ptr2,*ptr3;
      int p;
      sum = Create(scn1);
      for(p=0; p<scn1->n; p+=8){
	//sum->data[p] = scn1->data[p] + scn2->data[p];
	ptr1 = (v8hi *)(scn1->data + p);
	ptr2 = (v8hi *)(scn2->data + p);
	ptr3 = (v8hi *)(sum->data  + p);
	*ptr3 = *ptr1 + *ptr2;
      } 
      return sum;
    }

    sScene16 *Add(sScene16 *scn, ushort value){
      sScene16 *sum=NULL;
      v8hi *ptr1,*ptr2;
      v8hi v;
      int p,i;
      sum = Create(scn);
      for(i=0; i<8; i++)
	((ushort *)(&v))[i] = value;

      for(p=0; p<scn->n; p+=8){
	//sum->data[p] = scn->data[p] + value;
	ptr1 = (v8hi *)(scn->data + p);
	ptr2 = (v8hi *)(sum->data + p);
	*ptr2 = *ptr1 + v;
      }
      sum->maxval = scn->maxval + value;
      return sum;
    }

    //inplace version.
    void     Addinplace(sScene16 *scn1, sScene16 *scn2){
      v8hi *ptr1,*ptr2;
      int p;
      for(p=0; p<scn1->n; p+=8){
	//scn1->data[p] = scn1->data[p] + scn2->data[p];
	ptr1 = (v8hi *)(scn1->data + p);
	ptr2 = (v8hi *)(scn2->data + p);
	*ptr1 = *ptr1 + *ptr2;
      }
    }

    sScene16 *Mult(sScene16 *scn1, sScene16 *scn2){
      sScene16 *mscn=NULL;
      v8hi *ptr1,*ptr2;
      int p;
      mscn = Clone(scn2);
      for(p=0; p<scn1->n; p+=8){
	//mscn->data[p]*=scn1->data[p];
	ptr1 = (v8hi *)(scn1->data + p);
	ptr2 = (v8hi *)(mscn->data + p);
	(*ptr2) *= (*ptr1);
      }
      return(mscn);
    }


    sScene16 *Or(sScene16 *scn1, sScene16 *scn2){
      sScene16 *scnor=NULL;
      int p;
      scnor = Create(scn1);
      for(p=0; p<scn1->n; p++)
	scnor->data[p] = MAX(scn1->data[p],scn2->data[p]);

      return scnor;
    }

    //inplace version.
    void     Orinplace(sScene16 *scn1, sScene16 *scn2){
      int p;
      for(p=0; p<scn1->n; p++)
	scn1->data[p] = MAX(scn1->data[p],scn2->data[p]);
    }

    sScene16 *And(sScene16 *scn1, sScene16 *scn2){
      sScene16 *scnand=NULL;
      int p;
      scnand = Create(scn1);
      for(p=0; p<scn1->n; p++)
	scnand->data[p] = MIN(scn1->data[p],scn2->data[p]);
      
      return scnand;
    }

    sScene16 *XOr(sScene16 *scn1, sScene16 *scn2){
      sScene16 *scnxor=NULL;
      int p;
      scnxor = Create(scn1);
      for(p=0; p<scn1->n; p++){
	//scnxor->data[p] = (MAX(scn1->data[p],scn2->data[p])-
	//		     MIN(scn1->data[p],scn2->data[p]));
	scnxor->data[p] = ((scn1->data[p]>scn2->data[p])?
			   (scn1->data[p]-scn2->data[p]):
			   (scn2->data[p]-scn1->data[p]));
      }
      return scnxor;
    }

    sScene16 *Complement(sScene16 *scn){
      sScene16 *cscn=NULL;
      v8hi *ptr1,*ptr2;
      v8hi v;
      int p,i;
      ushort Imax;
      cscn = Create(scn);
      Imax = GetMaximumValue(scn);

      for(i=0; i<8; i++)
	((ushort *)(&v))[i] = Imax;

      for(p=0; p<scn->n; p+=8){
	//cscn->data[p] = Imax - scn->data[p];
	ptr1 = (v8hi *)(scn->data  + p);
	ptr2 = (v8hi *)(cscn->data + p);
	*ptr2 = v - *ptr1;
      }
      return(cscn);
    }

  } //end Scene16 namespace


  namespace Scene8{

    sScene8 *Sub(sScene8 *scn1, sScene8 *scn2){
      sScene8 *diff=NULL;
      v16qi *ptr1,*ptr2,*ptr3;
      int p;
      diff = Create(scn1);
      for(p=0; p<scn1->n; p+=16){
	//diff->data[p] = scn1->data[p] - scn2->data[p];
	ptr1 = (v16qi *)(scn1->data + p);
	ptr2 = (v16qi *)(scn2->data + p);
	ptr3 = (v16qi *)(diff->data + p);
	*ptr3 = *ptr1 - *ptr2;
      }
      return diff;
    }

    //inplace version.
    void    Subinplace(sScene8 *scn1, sScene8 *scn2){
      v16qi *ptr1,*ptr2;
      int p;
      for(p=0; p<scn1->n; p+=16){
	//scn1->data[p] = scn1->data[p] - scn2->data[p];
	ptr1 = (v16qi *)(scn1->data + p);
	ptr2 = (v16qi *)(scn2->data + p);
	*ptr1 = *ptr1 - *ptr2;
      }
    }

    sScene8 *Add(sScene8 *scn1, sScene8 *scn2){
      sScene8 *sum=NULL;
      v16qi *ptr1,*ptr2,*ptr3;
      int p;
      sum = Create(scn1);
      for(p=0; p<scn1->n; p+=16){
	//sum->data[p] = scn1->data[p] + scn2->data[p];
	ptr1 = (v16qi *)(scn1->data + p);
	ptr2 = (v16qi *)(scn2->data + p);
	ptr3 = (v16qi *)(sum->data  + p);
	*ptr3 = *ptr1 + *ptr2;
      } 
      return sum;
    }

    sScene8 *Add(sScene8 *scn, uchar value){
      sScene8 *sum=NULL;
      v16qi *ptr1,*ptr2;
      v16qi v;
      int p,i;
      sum = Create(scn);
      for(i=0; i<16; i++)
	((uchar *)(&v))[i] = value;

      for(p=0; p<scn->n; p+=16){
	//sum->data[p] = scn->data[p] + value;
	ptr1 = (v16qi *)(scn->data + p);
	ptr2 = (v16qi *)(sum->data + p);
	*ptr2 = *ptr1 + v;
      }
      sum->maxval = scn->maxval + value;
      return sum;
    }

    //inplace version.
    void    Addinplace(sScene8 *scn1, sScene8 *scn2){
      v16qi *ptr1,*ptr2;
      int p;
      for(p=0; p<scn1->n; p+=16){
	//scn1->data[p] = scn1->data[p] + scn2->data[p];
	ptr1 = (v16qi *)(scn1->data + p);
	ptr2 = (v16qi *)(scn2->data + p);
	*ptr1 = *ptr1 + *ptr2;
      }
    }

    sScene8 *Mult(sScene8 *scn1, sScene8 *scn2){
      sScene8 *mscn=NULL;
      v16qi *ptr1,*ptr2;
      int p;
      mscn = Clone(scn2);
      for(p=0; p<scn1->n; p+=16){
	//mscn->data[p]*=scn1->data[p];
	ptr1 = (v16qi *)(scn1->data + p);
	ptr2 = (v16qi *)(mscn->data + p);
	(*ptr2) *= (*ptr1);
      }
      return(mscn);
    }


    sScene8 *Or(sScene8 *scn1, sScene8 *scn2){
      sScene8 *scnor=NULL;
      int p;
      scnor = Create(scn1);
      for(p=0; p<scn1->n; p++)
	scnor->data[p] = MAX(scn1->data[p],scn2->data[p]);

      return scnor;
    }

    //inplace version.
    void     Orinplace(sScene8 *scn1, sScene8 *scn2){
      int p;
      for(p=0; p<scn1->n; p++)
	scn1->data[p] = MAX(scn1->data[p],scn2->data[p]);
    }

    sScene8 *And(sScene8 *scn1, sScene8 *scn2){
      sScene8 *scnand=NULL;
      int p;
      scnand = Create(scn1);
      for(p=0; p<scn1->n; p++)
	scnand->data[p] = MIN(scn1->data[p],scn2->data[p]);
      
      return scnand;
    }

    sScene8 *XOr(sScene8 *scn1, sScene8 *scn2){
      sScene8 *scnxor=NULL;
      int p;
      scnxor = Create(scn1);
      for(p=0; p<scn1->n; p++){
	//scnxor->data[p] = (MAX(scn1->data[p],scn2->data[p])-
	//		     MIN(scn1->data[p],scn2->data[p]));
	scnxor->data[p] = ((scn1->data[p]>scn2->data[p])?
			   (scn1->data[p]-scn2->data[p]):
			   (scn2->data[p]-scn1->data[p]));
      }
      return scnxor;
    }


    /*
    sScene8 *Or(sScene8 *scn1, sScene8 *scn2){
      sScene8 *scnor=NULL;
      v16qi *ptr1,*ptr2,*ptr3;
      int p;
      scnor = Create(scn1);
      for(p=0; p<scn1->n; p+=16){
	//scnor->data[p] = MAX(scn1->data[p],scn2->data[p]);
	ptr1 = (v16qi *)(scn1->data  + p);
	ptr2 = (v16qi *)(scn2->data  + p);
	ptr3 = (v16qi *)(scnor->data + p);
	*ptr3 = __builtin_ia32_pmaxub128(*ptr1, *ptr2);
      }
      return scnor;
    }

    //inplace version.
    void    Orinplace(sScene8 *scn1, sScene8 *scn2){
      v16qi *ptr1,*ptr2;
      int p;
      for(p=0; p<scn1->n; p+=16){
	//scn1->data[p] = MAX(scn1->data[p],scn2->data[p]);
	ptr1 = (v16qi *)(scn1->data + p);
	ptr2 = (v16qi *)(scn2->data + p);
	*ptr1 = __builtin_ia32_pmaxub128(*ptr1, *ptr2);
      }
    }

    sScene8 *And(sScene8 *scn1, sScene8 *scn2){
      sScene8 *scnand=NULL;
      v16qi *ptr1,*ptr2,*ptr3;
      int p;
      scnand = Create(scn1);
      for(p=0; p<scn1->n; p+=16){
	//scnand->data[p] = MIN(scn1->data[p],scn2->data[p]);
	ptr1 = (v16qi *)(scn1->data   + p);
	ptr2 = (v16qi *)(scn2->data   + p);
	ptr3 = (v16qi *)(scnand->data + p);
	*ptr3 = __builtin_ia32_pminub128(*ptr1, *ptr2);	
      }
      return scnand;
    }

    sScene8 *XOr(sScene8 *scn1, sScene8 *scn2){
      sScene8 *scnxor=NULL;
      v16qi *ptr1,*ptr2,*ptr3;
      v16qi min,max;
      int p;
      scnxor = Create(scn1);
      for(p=0; p<scn1->n; p+=16){
	//scnxor->data[p] = (MAX(scn1->data[p],scn2->data[p])-
	//		     MIN(scn1->data[p],scn2->data[p]));
	ptr1 = (v16qi *)(scn1->data   + p);
	ptr2 = (v16qi *)(scn2->data   + p);
	ptr3 = (v16qi *)(scnxor->data + p);
	max  = __builtin_ia32_pmaxub128(*ptr1, *ptr2);
	min  = __builtin_ia32_pminub128(*ptr1, *ptr2);
	*ptr3 = max - min;
      }
      return scnxor;
    }
    */

    sScene8 *Complement(sScene8 *scn){
      sScene8 *cscn=NULL;
      v16qi *ptr1,*ptr2;
      v16qi v;
      int p,i;
      uchar Imax;
      cscn = Create(scn);
      Imax = GetMaximumValue(scn);

      for(i=0; i<16; i++)
	((uchar *)(&v))[i] = Imax;

      for(p=0; p<scn->n; p+=16){
	//cscn->data[p] = Imax - scn->data[p];
	ptr1 = (v16qi *)(scn->data  + p);
	ptr2 = (v16qi *)(cscn->data + p);
	*ptr2 = v - *ptr1;
      }
      return(cscn);
    }

  } //end Scene8 namespace


  namespace Scene{

    sScene *Sub(sScene *scn1, sScene *scn2){
      sScene *diff;

      if(scn1->nbits!=scn2->nbits)
	gft::Error((char *)"Incompatible number of bits",
		   (char *)"Scene::Sub");

      diff = (sScene *) calloc(1,sizeof(sScene));
      if(diff == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::Sub");

      diff->nbits = scn1->nbits;
      switch(scn1->nbits){
      case  8: 
	diff->ptr.scn8  = Scene8::Sub(scn1->ptr.scn8,
				      scn2->ptr.scn8);  
	break;
      case 16: 
	diff->ptr.scn16 = Scene16::Sub(scn1->ptr.scn16,
				       scn2->ptr.scn16); 
	break;
      case 32: 
	diff->ptr.scn32 = Scene32::Sub(scn1->ptr.scn32,
				       scn2->ptr.scn32);
	break;
      }
      return diff;
    }

    //inplace version.
    void   Subinplace(sScene *scn1, sScene *scn2){
      if(scn1->nbits!=scn2->nbits)
	gft::Error((char *)"Incompatible number of bits",
		   (char *)"Scene::Subinplace");

      switch(scn1->nbits){
      case  8: 
	Scene8::Subinplace(scn1->ptr.scn8,
			   scn2->ptr.scn8);  
	break;
      case 16: 
	Scene16::Subinplace(scn1->ptr.scn16,
			    scn2->ptr.scn16); 
	break;
      case 32: 
	Scene32::Subinplace(scn1->ptr.scn32,
			    scn2->ptr.scn32);
	break;
      }
    }

    sScene *Add(sScene *scn1, sScene *scn2){
      sScene *sum;

      if(scn1->nbits!=scn2->nbits)
	gft::Error((char *)"Incompatible number of bits",
		   (char *)"Scene::Add");

      sum = (sScene *) calloc(1,sizeof(sScene));
      if(sum == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::Add");

      sum->nbits = scn1->nbits;
      switch(scn1->nbits){
      case  8: 
	sum->ptr.scn8  = Scene8::Add(scn1->ptr.scn8,
				     scn2->ptr.scn8);  
	break;
      case 16: 
	sum->ptr.scn16 = Scene16::Add(scn1->ptr.scn16,
				      scn2->ptr.scn16); 
	break;
      case 32: 
	sum->ptr.scn32 = Scene32::Add(scn1->ptr.scn32,
				      scn2->ptr.scn32);
	break;
      }
      return sum;
    }


    sScene *Add(sScene *scn, int value){
      sScene *sum;

      sum = (sScene *) calloc(1,sizeof(sScene));
      if(sum == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::Add");

      sum->nbits = scn->nbits;
      switch(scn->nbits){
      case  8: 
	sum->ptr.scn8  = Scene8::Add(scn->ptr.scn8, value);
	break;
      case 16: 
	sum->ptr.scn16 = Scene16::Add(scn->ptr.scn16, value); 
	break;
      case 32: 
	sum->ptr.scn32 = Scene32::Add(scn->ptr.scn32, value);
	break;
      }
      return sum;
    }

    //inplace version.
    void   Addinplace(sScene *scn1, sScene *scn2){
      if(scn1->nbits!=scn2->nbits)
	gft::Error((char *)"Incompatible number of bits",
		   (char *)"Scene::Addinplace");

      switch(scn1->nbits){
      case  8: 
	Scene8::Addinplace(scn1->ptr.scn8,
			   scn2->ptr.scn8);  
	break;
      case 16: 
	Scene16::Addinplace(scn1->ptr.scn16,
			    scn2->ptr.scn16); 
	break;
      case 32: 
	Scene32::Addinplace(scn1->ptr.scn32,
			    scn2->ptr.scn32);
	break;
      }
    }


    sScene *Mult(sScene *scn1, sScene *scn2){
      sScene *mscn;

      if(scn1->nbits!=scn2->nbits)
	gft::Error((char *)"Incompatible number of bits",
		   (char *)"Scene::Mult");

      mscn = (sScene *) calloc(1,sizeof(sScene));
      if(mscn == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::Mult");

      mscn->nbits = scn1->nbits;
      switch(scn1->nbits){
      case  8: 
	mscn->ptr.scn8  = Scene8::Mult(scn1->ptr.scn8,
				       scn2->ptr.scn8);  
	break;
      case 16: 
	mscn->ptr.scn16 = Scene16::Mult(scn1->ptr.scn16,
					scn2->ptr.scn16); 
	break;
      case 32: 
	mscn->ptr.scn32 = Scene32::Mult(scn1->ptr.scn32,
					scn2->ptr.scn32);
	break;
      }
      return mscn;
    }


    sScene *Or(sScene *scn1, sScene *scn2){
      sScene *scnor=NULL;

      if(scn1->nbits!=scn2->nbits)
	gft::Error((char *)"Incompatible number of bits",
		   (char *)"Scene::Or");

      scnor = (sScene *) calloc(1,sizeof(sScene));
      if(scnor == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::Or");

      scnor->nbits = scn1->nbits;
      switch(scn1->nbits){
      case  8: 
	scnor->ptr.scn8  = Scene8::Or(scn1->ptr.scn8,
				      scn2->ptr.scn8);  
	break;
      case 16: 
	scnor->ptr.scn16 = Scene16::Or(scn1->ptr.scn16,
				       scn2->ptr.scn16); 
	break;
      case 32: 
	scnor->ptr.scn32 = Scene32::Or(scn1->ptr.scn32,
				       scn2->ptr.scn32);
	break;
      }
      return scnor;
    }


    //inplace version.
    void   Orinplace(sScene *scn1, sScene *scn2){
      if(scn1->nbits!=scn2->nbits)
	gft::Error((char *)"Incompatible number of bits",
		   (char *)"Scene::Orinplace");

      switch(scn1->nbits){
      case  8: 
	Scene8::Orinplace(scn1->ptr.scn8,
			  scn2->ptr.scn8);  
	break;
      case 16: 
	Scene16::Orinplace(scn1->ptr.scn16,
			   scn2->ptr.scn16); 
	break;
      case 32: 
	Scene32::Orinplace(scn1->ptr.scn32,
			   scn2->ptr.scn32);
	break;
      }
    }


    sScene *And(sScene *scn1, sScene *scn2){
      sScene *scnand=NULL;

      if(scn1->nbits!=scn2->nbits)
	gft::Error((char *)"Incompatible number of bits",
		   (char *)"Scene::And");

      scnand = (sScene *) calloc(1,sizeof(sScene));
      if(scnand == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::And");

      scnand->nbits = scn1->nbits;
      switch(scn1->nbits){
      case  8: 
	scnand->ptr.scn8  = Scene8::And(scn1->ptr.scn8,
					scn2->ptr.scn8);  
	break;
      case 16: 
	scnand->ptr.scn16 = Scene16::And(scn1->ptr.scn16,
					 scn2->ptr.scn16); 
	break;
      case 32: 
	scnand->ptr.scn32 = Scene32::And(scn1->ptr.scn32,
					 scn2->ptr.scn32);
	break;
      }
      return scnand;
    }


    sScene *XOr(sScene *scn1, sScene *scn2){
      sScene *scnxor=NULL;

      if(scn1->nbits!=scn2->nbits)
	gft::Error((char *)"Incompatible number of bits",
		   (char *)"Scene::XOr");

      scnxor = (sScene *) calloc(1,sizeof(sScene));
      if(scnxor == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::XOr");

      scnxor->nbits = scn1->nbits;
      switch(scn1->nbits){
      case  8: 
	scnxor->ptr.scn8  = Scene8::XOr(scn1->ptr.scn8,
					scn2->ptr.scn8);  
	break;
      case 16: 
	scnxor->ptr.scn16 = Scene16::XOr(scn1->ptr.scn16,
					 scn2->ptr.scn16); 
	break;
      case 32: 
	scnxor->ptr.scn32 = Scene32::XOr(scn1->ptr.scn32,
					 scn2->ptr.scn32);
	break;
      }
      return scnxor;
    }


    sScene *Complement(sScene *scn){
      sScene *cscn=NULL;

      cscn = (sScene *) calloc(1,sizeof(sScene));
      if(cscn == NULL) 
	gft::Error((char *)MSG1,(char *)"Scene::Complement");

      cscn->nbits = scn->nbits;
      switch(scn->nbits){
      case  8: 
	cscn->ptr.scn8  = Scene8::Complement(scn->ptr.scn8);
	break;
      case 16: 
	cscn->ptr.scn16 = Scene16::Complement(scn->ptr.scn16);
	break;
      case 32: 
	cscn->ptr.scn32 = Scene32::Complement(scn->ptr.scn32);
	break;
      }
      return cscn;
    }


  } //end Scene namespace


} //end gft namespace

