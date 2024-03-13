
#include "gft.h"

#include <string>

#define IMG_COLOR     1
#define IMG_GRAYSCALE 0

gft::sCImage *ReadAnyCImage(char *file){
  gft::sCImage *cimg;
  char command[512];
  int s;

  s = strlen(file);
  if(strcasecmp(&file[s-3], "ppm") == 0){
    cimg = gft::CImage::Read(file);
  }
  else{
    sprintf(command, "convert %s cimage_tmp.ppm", file);
    system(command);
    cimg = gft::CImage::Read("cimage_tmp.ppm");
    system("rm cimage_tmp.ppm");
  }
  return cimg;
}


gft::sImage32 *ReadAnyImage(char *file){
  gft::sImage32 *img;
  char command[512];
  int s;

  s = strlen(file);
  if(strcasecmp(&file[s-3], "pgm") == 0){
    img = gft::Image32::Read(file);
  }
  else{
    sprintf(command, "convert %s image_tmp.pgm", file);
    system(command);
    img = gft::Image32::Read("image_tmp.pgm");
    system("rm image_tmp.pgm");
  }
  return img;
}


int main(int argc, char **argv){
  gft::sCImage *cimg=NULL, *spixels, *tmp;
  gft::sImage32 *img=NULL, *label;
  char filename[512];
  int k, color = 0xFF0000, black=0x000000, type;
  double alpha;
  //clock_t end, start;
  //double totaltime;

  if(argc < 5){
    fprintf(stdout,"usage:\n");
    fprintf(stdout,"superpixels <type> <image> <k> <alpha>\n");
    fprintf(stdout,"\t type.... 0: grayscale\n");
    fprintf(stdout,"\t          1: color\n");
    fprintf(stdout,"\t k....... The desired number of superpixels\n");
    fprintf(stdout,"\t alpha... The relative importance between color similarity and spatial proximity.\n");
    fprintf(stdout,"\t          It can be in the range [0.01, 0.10]\n");
    exit(0);
  }

  type = atoi(argv[1]);
  strcpy(filename, argv[2]);
  k = atoi(argv[3]);
  alpha = atof(argv[4]);

  if(type == IMG_GRAYSCALE){
    img  = ReadAnyImage(filename);
    cimg = gft::CImage::Clone(img);
  }
  else
    cimg = ReadAnyCImage(filename);


  //-------------------------------------------------
  // IFT-SLIC:
  //start = clock();
  if(type == IMG_GRAYSCALE)  
    label = gft::Superpixels::IFT_SLIC(img, k, alpha, //12.0,
				       5.0, 2.0, 10);
  else
    label = gft::Superpixels::IFT_SLIC(cimg, k, alpha, //12.0,
				       5.0, 2.0, 10);
  //end = clock();
  //totaltime = ((double)(end-start))/CLOCKS_PER_SEC;
  //printf("IFT-SLIC Time: %f sec\n", totaltime);

  tmp = gft::Highlight::CWideLabels(cimg, label, 2.0, &black, 0.0,
				    true, true);
  spixels = gft::Highlight::CWideLabels(tmp, label, 1.0, &color, 0.0,
					true, true);
  gft::CImage::Destroy(&tmp);

  printf("N = %d\n", gft::Superpixels::GetNumberOfSuperPixels(label));

  gft::Image32::Write(label, (char *)"./out/label_ift.pgm");
  gft::CImage::Write(spixels, (char *)"./out/spixels_ift.ppm");

  gft::Image32::Destroy(&label);
  gft::CImage::Destroy(&spixels);
  if(cimg != NULL) gft::CImage::Destroy(&cimg);
  if(img != NULL)  gft::Image32::Destroy(&img);
  return 0;
}



