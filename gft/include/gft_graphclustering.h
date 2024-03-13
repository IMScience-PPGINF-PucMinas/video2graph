
#ifndef _GFT_GRAPHCLUSTERING_H_
#define _GFT_GRAPHCLUSTERING_H_

#include "gft_common.h"
#include "gft_graph.h"
#include "gft_pqueue32.h"
#include "gft_queue.h"
#include "gft_ift.h"

namespace gft{
  namespace Graph{

    struct sClusteringTree {
      sGraph *G;
      int right;
      int left;
      int seed;
    };

    struct sMergeHistory {
      int p;
      int q;
      int energy;
    };
    
    //Divisive nearest-neighbor (single-linkage) algorithm:
    //ClusteringTree *DivisiveClusteringByMST(Graph *G, int c);
    int *DivisiveClusteringByMST(sGraph *G, int c);

    //Felzenszwalb and Huttenlocher 2004:
    int *ClusteringByMST2(sGraph *G, float k);
    
    int *DivisiveClusteringByOIFT(sGraph *G, int c);
    sMergeHistory *DivisiveClusteringByOIFT(sGraph *G);

    gft::sImage32 *ComputeEnegyMap_UOIFT(gft::sImage32 *spixels, sGraph *G);
    
    //Tiebreaking by DCC size:
    int *DivisiveClusteringByOIFT_2(sGraph *G, int c);
    sMergeHistory *DivisiveClusteringByOIFT_2(sGraph *G);

    //Tiebreaking by DCC size embedded in the cost:
    int *DivisiveClusteringByOIFT_3(sGraph *G, int c);
    
    void WriteMergeHistory(sMergeHistory *MH, int n, char *file);
    
  } //end Graph namespace

  
} //end gft namespace

    
#endif

