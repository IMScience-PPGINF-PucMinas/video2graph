
#ifndef _GFT_CTREE_H_
#define _GFT_CTREE_H_

#include "gft_common.h"
#include "gft_graph.h"
#include "gft_pqueue32.h"

namespace gft{

  typedef enum {height, area, volume} AttributeType;
  
  namespace CTree{

    struct sCTreeNode {
      int level;
      //Attributes:
      int height;
      int area;
      int volume;
      sCTreeNode *parent;
      int nsons;
      sCTreeNode **sons;
    };
    
    struct sCTree {
      int nnodes;
      int arraysize;
      sCTreeNode *nodes;
      int nleaves;

      int nspels;
      int *mapping; //spel to ctreenode
    };


    sCTree *Create(int size);
    void Destroy(sCTree **ct);

    int GetNodeIndex(sCTree *ct, sCTreeNode *node);
    
    void Print(sCTree *ct);

    sCTree *EdgeBasedMinTree(sGraph *graph, int Lmin);
    sCTree *EdgeBasedMinTree(sImageGraph *graph, int Lmin);
    
    void ComputeHeight(sCTree *ct);
    void ComputeArea(sCTree *ct, int *area);
    void ComputeVolume(sCTree *ct);

    int *ComputeExtinctionValue(sCTree *ct, AttributeType type);
    
  } //end CTree namespace

  typedef CTree::sCTree sCTree;
  typedef CTree::sCTreeNode sCTreeNode;
  
} //end gft namespace
    


    
#endif

