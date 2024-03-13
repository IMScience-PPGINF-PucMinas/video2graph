
#include "gft.h"

#include <string>

int main(int argc, char **argv){
  gft::sScene32 *scn, *label;
  char filename[512];
  int k;
  double alpha;
  clock_t end, start;
  double totaltime;

  if(argc < 4){
    fprintf(stdout,"usage:\n");
    fprintf(stdout,"supervoxels <image> <k> <alpha>\n");
    fprintf(stdout,"\t k....... The desired number of superpixels\n");
    fprintf(stdout,"\t alpha... The relative importance between color similarity and spatial proximity.\n");
    fprintf(stdout,"\t          It can be in the range [0.01, 0.10]\n");
    exit(0);
  }

  strcpy(filename, argv[1]);
  k = atoi(argv[2]);
  alpha = atof(argv[3]);

  scn = gft::Scene32::Read(filename);

  //-------------------------------------------------
  // IFT-SLIC:
  start = clock();

  label = gft::Superpixels::IFT_SLIC(scn, k, alpha, 5.0, 2.0, 10);
  
  end = clock();
  totaltime = ((double)(end-start))/CLOCKS_PER_SEC;
  printf("IFT-SLIC Time: %f sec\n", totaltime);

  gft::Scene32::Write(label, (char *)"./out/label_ift.scn.bz2");

  gft::Scene32::Destroy(&label);
  gft::Scene32::Destroy(&scn);
  return 0;
}



