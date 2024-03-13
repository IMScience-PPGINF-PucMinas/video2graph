
#include "gft.h"


int main(){
  gft::sImage32 *I,*AC;
  gft::sAdjRel *A;
  char filename[512];
  int T;

  printf("Enter the filename: ");  
  scanf("%s", filename);
  printf("Area threshold T: ");
  scanf("%d",&T);
  I = gft::Image32::Read(filename); //(char *)"./dat/cheese.pgm");
  A = gft::AdjRel::Circular(1.0);
  AC = gft::Image32::AreaClosing(A, I, T);

  gft::Image32::Write(AC, (char *)"AC.pgm");

  gft::Image32::Destroy(&I);
  gft::Image32::Destroy(&AC);
  gft::AdjRel::Destroy(&A);
  return 0;
}

