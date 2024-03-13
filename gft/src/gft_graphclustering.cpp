
#include "gft_graphclustering.h"

namespace gft{
  namespace Graph{


    int *DivisiveClusteringByMST(sGraph *G, int c){
      sPQueue32 *Q;
      sQueue *F;
      int *label, *pred, *cost;
      int p,q, Wmax, i;
      int S[2];
      S[0] = 1;
      S[1] = 0;
      label = gft::AllocIntArray(G->nnodes);
      pred  = gft::AllocIntArray(G->nnodes);
      cost  = gft::AllocIntArray(G->nnodes);
      for(p = 0; p < G->nnodes; p++)
	label[p] = NIL;
      label[0] = 0;
      gft::ift::IFT_fw(G, S, label, cost, pred);
      
      Wmax = GetMaximumArc(G);
      Q = PQueue32::Create(Wmax+2, G->nnodes, cost);

      for(p = 1; p < G->nnodes; p++)
	PQueue32::FastInsertElem(Q, p);
      i = 0;
      while(!PQueue32::IsEmpty(Q) && i < c-1) {
	p = PQueue32::FastRemoveMaxFIFO(Q);
	pred[p] = NIL;
	i++;
      }
      PQueue32::Destroy(&Q);

      F = Queue::Create(G->nnodes);
      for(p = 0; p < G->nnodes; p++)
	label[p] = NIL;
      i = 0;
      for(p = 0; p < G->nnodes; p++){
	if(pred[p] == NIL){
	  Queue::Push(F, p);
	  label[p] = i;
	  i++;
	}
      }

      while(!Queue::IsEmpty(F)){
	p = Queue::Pop(F);
	for(i = 0; i < G->nodes[p].outdegree; i++){
	  q = G->nodes[p].adjList[i];
	  if(pred[q] == p){
	    label[q] = label[p];
	    Queue::Push(F, q);
	  }
	}
      }
      Queue::Destroy(&F);
      gft::FreeIntArray(&pred);
      gft::FreeIntArray(&cost);
      return label;
    }


    int GetNodeIndex(sGraph *G, int id){
      int p;
      for(p = 0; p < G->nnodes; p++)
	if(G->nodes[p].id == id)
	  return p;
      return NIL;
    }



    int ComputeDCCsize(sGraph *G, int energy, int p){
      gft::sQueue *Q;
      gft::sBMap *B;
      int size, q,i,w;
      Q = gft::Queue::Create(G->nnodes);
      B = gft::BMap::Create(G->nnodes);
      gft::Queue::Push(Q, p);
      gft::BMap::Set1(B, p);
      size = 1;
      while(!gft::Queue::IsEmpty(Q)){
	p = gft::Queue::Pop(Q);
	for(i = 0; i < G->nodes[p].outdegree; i++){
	  q = G->nodes[p].adjList[i];
	  w = G->nodes[p].Warcs[i];
	  if(w != NIL && w < energy && gft::BMap::Get(B, q)==0){
	    gft::Queue::Push(Q, q);
	    gft::BMap::Set1(B, q);
	    size++;  
	  }
	}
      }
      gft::Queue::Destroy(&Q);
      gft::BMap::Destroy(&B);
      return size;
    }



    gft::sImage32 *ComputeEnegyMap_UOIFT(gft::sImage32 *spixels, sGraph *G){
      gft::sImage32 *energy;
      sGraph *T;
      int S[3];
      int nnodes, p, lb;
      int *label, *cost;
      S[0] = 1;
      S[1] = 0;
      nnodes = G->nnodes;
      energy = gft::Image32::Create(spixels);
      label = gft::AllocIntArray(nnodes);
      cost  = gft::AllocIntArray(nnodes);
      for(p = 0; p < nnodes; p++)
	label[p] = NIL;
      label[0] = 0;
      T = Transpose(G);
      ift::IFT_fw(T, S, label, cost);
      Destroy(&T);

      for(p = 0; p < spixels->n; p++){
	lb = spixels->data[p];
	energy->data[p] = cost[lb];
      }
      gft::FreeIntArray(&label);
      gft::FreeIntArray(&cost);
      return energy;
    }
   

    int *DivisiveClusteringByOIFT(sGraph *G, int c){
      sPQueue32 *Q;
      sClusteringTree *tree;
      int *label, *cost, *oift_label;
      sGraph *T;
      int n,i,Wmax,lb,p,q,t,l,nnodes;
      int S[3];
      int *hist, *mapping;
      tree = (sClusteringTree *)malloc(sizeof(sClusteringTree)*c*2);
      if(tree == NULL){
	printf("Error: DivisiveClusteringByOIFT\n");
	exit(1);
      }
      for(i = 0; i < c*2; i++){
	tree[i].G = NULL;
	tree[i].right = NIL;
	tree[i].left = NIL;
      }
      tree[0].G = gft::Graph::Clone(G);
      tree[0].seed = 0;
      n = 1;

      S[0] = 1;
      S[1] = 0;
      nnodes = G->nnodes;
      label = gft::AllocIntArray(nnodes);
      cost  = gft::AllocIntArray(nnodes);
      for(p = 0; p < nnodes; p++)
	label[p] = NIL;
      label[0] = 0;
      T = Transpose(G);
      ift::IFT_fw(T, S, label, cost);
      Destroy(&T);

      Wmax = GetMaximumArc(G);
      Q = PQueue32::Create(Wmax+2, nnodes, cost);

      for(p = 1; p < G->nnodes; p++)
	PQueue32::FastInsertElem(Q, p);
      i = 0;
      while(!PQueue32::IsEmpty(Q) && i < c-1) {
	p = PQueue32::FastRemoveMaxFIFO(Q);
	
	//printf("2) size: %d, cost: %d\n",ComputeDCCsize(G, cost[p], p), cost[p]);

	lb = label[p];
	oift_label = gft::AllocIntArray(tree[lb].G->nnodes);
	T = Transpose(tree[lb].G);
	S[0] = 2;
	S[1] = tree[lb].seed;
	S[2] = GetNodeIndex(tree[lb].G, p);
	q = tree[lb].G->nodes[S[1]].id;
	for(t = 0; t < tree[lb].G->nnodes; t++)
	  oift_label[t] = NIL;
	oift_label[S[1]] = 0;
	oift_label[S[2]] = 1;
	ift::OIFT(tree[lb].G, T, S, oift_label);
	Destroy(&T);

	tree[n].G = Split(&(tree[lb].G), oift_label, 1);
	tree[n+1].G = tree[lb].G;
	tree[n].seed   = GetNodeIndex(tree[n].G,   p);
	tree[n+1].seed = GetNodeIndex(tree[n+1].G, q);
	tree[lb].G = NULL;
	tree[lb].left = n;
	tree[lb].right = n+1;
	for(t = 0; t < tree[n].G->nnodes; t++)
	  label[tree[n].G->nodes[t].id] = n;
	for(t = 0; t < tree[n+1].G->nnodes; t++)
	  label[tree[n+1].G->nodes[t].id] = n+1;
	n += 2;
	gft::FreeIntArray(&oift_label);

	i++;
      }
      PQueue32::Destroy(&Q);
      for(i = 0; i < c*2; i++){
	if(tree[i].G != NULL)
	  gft::Graph::Destroy(&tree[i].G);
      }
      free(tree);
      gft::FreeIntArray(&cost);

      //Fixing the labels to match the number of partitions:
      hist = gft::AllocIntArray(n);
      for(p = 0; p < nnodes; p++){
	hist[label[p]] += 1;
      }
      mapping = gft::AllocIntArray(n);
      i = 0;
      for(l = 0; l < n; l++){
	if(hist[l] > 0){
	  mapping[l] = i;
	  i++;
	}
      }
      for(p = 0; p < nnodes; p++)
	label[p] = mapping[label[p]];
      gft::FreeIntArray(&mapping);
      gft::FreeIntArray(&hist);
      
      return label;
    }



    int *DivisiveClusteringByOIFT_2(sGraph *G, int c){
      sPQueue32 *Q,*Q_size;
      sClusteringTree *tree;
      int *label, *cost, *oift_label,*DCC_size;
      sGraph *T;
      int n,i,Wmax,lb,p,q,t,l,nnodes,energy;
      int S[3];
      int *hist, *mapping;
      tree = (sClusteringTree *)malloc(sizeof(sClusteringTree)*c*2);
      if(tree == NULL){
	printf("Error: DivisiveClusteringByOIFT_2\n");
	exit(1);
      }
      for(i = 0; i < c*2; i++){
	tree[i].G = NULL;
	tree[i].right = NIL;
	tree[i].left = NIL;
      }
      tree[0].G = gft::Graph::Clone(G);
      tree[0].seed = 0;
      n = 1;

      S[0] = 1;
      S[1] = 0;
      nnodes = G->nnodes;
      label = gft::AllocIntArray(nnodes);
      cost  = gft::AllocIntArray(nnodes);
      DCC_size = gft::AllocIntArray(nnodes);
      for(p = 0; p < nnodes; p++)
	label[p] = NIL;
      label[0] = 0;
      T = Transpose(G);
      ift::IFT_fw(T, S, label, cost);
      Destroy(&T);

      Wmax = GetMaximumArc(G);
      Q = PQueue32::Create(Wmax+2, nnodes, cost);
      Q_size = PQueue32::Create(nnodes+2, nnodes, DCC_size);
      
      for(p = 1; p < G->nnodes; p++)
	PQueue32::FastInsertElem(Q, p);
      i = 0;
      while(!PQueue32::IsEmpty(Q) && i < c-1) {
	energy = PQueue32::FastGetMaxVal(Q);
	
	while(!PQueue32::IsEmpty(Q) &&
	      PQueue32::FastGetMaxVal(Q) == energy){
	  p = PQueue32::FastRemoveMaxFIFO(Q);
	  DCC_size[p] = ComputeDCCsize(G, energy, p);
	  PQueue32::FastInsertElem(Q_size, p);
	  //printf("1) energy: %d, size: %d, cost: %d\n",energy, DCC_size[p], cost[p]);
	}
	
	while(!PQueue32::IsEmpty(Q_size) && i < c-1){
	  p = PQueue32::FastRemoveMaxFIFO(Q_size);
	  //printf("2) energy: %d, size: %d, cost: %d\n",energy, DCC_size[p], cost[p]);

	  lb = label[p];
	  oift_label = gft::AllocIntArray(tree[lb].G->nnodes);
	  T = Transpose(tree[lb].G);
	  S[0] = 2;
	  S[1] = tree[lb].seed;
	  S[2] = GetNodeIndex(tree[lb].G, p);
	  q = tree[lb].G->nodes[S[1]].id;
	  for(t = 0; t < tree[lb].G->nnodes; t++)
	    oift_label[t] = NIL;
	  oift_label[S[1]] = 0;
	  oift_label[S[2]] = 1;
	  ift::OIFT(tree[lb].G, T, S, oift_label);
	  Destroy(&T);
	  
	  tree[n].G = Split(&(tree[lb].G), oift_label, 1);
	  tree[n+1].G = tree[lb].G;
	  tree[n].seed   = GetNodeIndex(tree[n].G,   p);
	  tree[n+1].seed = GetNodeIndex(tree[n+1].G, q);
	  tree[lb].G = NULL;
	  tree[lb].left = n;
	  tree[lb].right = n+1;
	  for(t = 0; t < tree[n].G->nnodes; t++)
	    label[tree[n].G->nodes[t].id] = n;
	  for(t = 0; t < tree[n+1].G->nnodes; t++)
	    label[tree[n+1].G->nodes[t].id] = n+1;
	  n += 2;
	  gft::FreeIntArray(&oift_label);
	  
	  i++;
	}
      }

      PQueue32::Destroy(&Q_size);
      PQueue32::Destroy(&Q);
      for(i = 0; i < c*2; i++){
	if(tree[i].G != NULL)
	  gft::Graph::Destroy(&tree[i].G);
      }
      free(tree);
      gft::FreeIntArray(&cost);
      gft::FreeIntArray(&DCC_size);
      
      //Fixing the labels to match the number of partitions:
      hist = gft::AllocIntArray(n);
      for(p = 0; p < nnodes; p++){
	hist[label[p]] += 1;
      }
      mapping = gft::AllocIntArray(n);
      i = 0;
      for(l = 0; l < n; l++){
	if(hist[l] > 0){
	  mapping[l] = i;
	  i++;
	}
      }
      for(p = 0; p < nnodes; p++)
	label[p] = mapping[label[p]];
      gft::FreeIntArray(&mapping);
      gft::FreeIntArray(&hist);
      
      return label;
    }



    int *IFT_fwT_dcc(sGraph *graph,
		     sGraph *graphT,
		     int *S,
		     int *label,
		     int *cost){
      sPQueue32 *Q;
      int tmp, w, Wmax;
      int n,p,q,i;
      int *dcc_size;
      n = graph->nnodes;
      dcc_size = gft::AllocIntArray(n);
      Wmax = Graph::GetMaximumArc(graph);
      Q = PQueue32::Create((Wmax+1)*n + 2, n, cost);
      
      for(p = 0; p < n; p++){
	if(label[p]==NIL) cost[p] = INT_MAX;
	else              cost[p] = 0;
      }
      
      if(S != NULL){
	for(i = 1; i <= S[0]; i++)
	  PQueue32::FastInsertElem(Q, S[i]);
      }
      else{
	for(p = 0; p < n; p++)
	  if(label[p]!=NIL)
	    PQueue32::FastInsertElem(Q, p);
      }
      
      while(!PQueue32::IsEmpty(Q)){
	p = PQueue32::FastRemoveMinFIFO(Q);
	
	for(i = 0; i < graphT->nodes[p].outdegree; i++){
	  q = graphT->nodes[p].adjList[i];
	  if(Q->L.elem[q].color != BLACK){
	    w = graphT->nodes[p].Warcs[i];
	    tmp = w*n + (n-ComputeDCCsize(graph, w, q));

	    if(tmp < cost[q]){
	      if(Q->L.elem[q].color == GRAY)
		PQueue32::FastRemoveElem(Q, q);
	      cost[q] = tmp;
	      label[q] = label[p];
	      PQueue32::FastInsertElem(Q, q);
	    }
	  }
	}
      }
      PQueue32::Destroy(&Q);
      for(p = 0; p < n; p++){
	dcc_size[p] = n-(cost[p]%n);
	cost[p] = cost[p]/n;
      }
      return dcc_size;
    }
        

    int *DivisiveClusteringByOIFT_3(sGraph *G, int c){
      sPQueue32 *Q,*Q_size;
      sClusteringTree *tree;
      int *label, *cost, *oift_label,*DCC_size;
      sGraph *T;
      int n,i,Wmax,lb,p,q,t,l,nnodes,energy;
      int S[3];
      int *hist, *mapping;
      tree = (sClusteringTree *)malloc(sizeof(sClusteringTree)*c*2);
      if(tree == NULL){
	printf("Error: DivisiveClusteringByOIFT_3\n");
	exit(1);
      }
      for(i = 0; i < c*2; i++){
	tree[i].G = NULL;
	tree[i].right = NIL;
	tree[i].left = NIL;
      }
      tree[0].G = gft::Graph::Clone(G);
      tree[0].seed = 0;
      n = 1;

      S[0] = 1;
      S[1] = 0;
      nnodes = G->nnodes;
      label = gft::AllocIntArray(nnodes);
      cost  = gft::AllocIntArray(nnodes);
      for(p = 0; p < nnodes; p++)
	label[p] = NIL;
      label[0] = 0;
      T = Transpose(G);
      DCC_size = IFT_fwT_dcc(G, T, S, label, cost);
      Destroy(&T);

      Wmax = GetMaximumArc(G);
      Q = PQueue32::Create(Wmax+2, nnodes, cost);
      Q_size = PQueue32::Create(nnodes+2, nnodes, DCC_size);
      
      for(p = 1; p < G->nnodes; p++)
	PQueue32::FastInsertElem(Q, p);
      i = 0;
      while(!PQueue32::IsEmpty(Q) && i < c-1) {
	energy = PQueue32::FastGetMaxVal(Q);
	
	while(!PQueue32::IsEmpty(Q) &&
	      PQueue32::FastGetMaxVal(Q) == energy){
	  p = PQueue32::FastRemoveMaxFIFO(Q);
	  //DCC_size[p] = ComputeDCCsize(G, energy, p);
	  PQueue32::FastInsertElem(Q_size, p);
	  //printf("1) energy: %d, size: %d, cost: %d\n",energy, DCC_size[p], cost[p]);
	}
	
	while(!PQueue32::IsEmpty(Q_size) && i < c-1){
	  p = PQueue32::FastRemoveMaxFIFO(Q_size);
	  //printf("2) energy: %d, size: %d, cost: %d\n",energy, DCC_size[p], cost[p]);

	  lb = label[p];
	  oift_label = gft::AllocIntArray(tree[lb].G->nnodes);
	  T = Transpose(tree[lb].G);
	  S[0] = 2;
	  S[1] = tree[lb].seed;
	  S[2] = GetNodeIndex(tree[lb].G, p);
	  q = tree[lb].G->nodes[S[1]].id;
	  for(t = 0; t < tree[lb].G->nnodes; t++)
	    oift_label[t] = NIL;
	  oift_label[S[1]] = 0;
	  oift_label[S[2]] = 1;
	  ift::OIFT(tree[lb].G, T, S, oift_label);
	  Destroy(&T);
	  
	  tree[n].G = Split(&(tree[lb].G), oift_label, 1);
	  tree[n+1].G = tree[lb].G;
	  tree[n].seed   = GetNodeIndex(tree[n].G,   p);
	  tree[n+1].seed = GetNodeIndex(tree[n+1].G, q);
	  tree[lb].G = NULL;
	  tree[lb].left = n;
	  tree[lb].right = n+1;
	  for(t = 0; t < tree[n].G->nnodes; t++)
	    label[tree[n].G->nodes[t].id] = n;
	  for(t = 0; t < tree[n+1].G->nnodes; t++)
	    label[tree[n+1].G->nodes[t].id] = n+1;
	  n += 2;
	  gft::FreeIntArray(&oift_label);
	  
	  i++;
	}
      }

      PQueue32::Destroy(&Q_size);
      PQueue32::Destroy(&Q);
      for(i = 0; i < c*2; i++){
	if(tree[i].G != NULL)
	  gft::Graph::Destroy(&tree[i].G);
      }
      free(tree);
      gft::FreeIntArray(&cost);
      gft::FreeIntArray(&DCC_size);
      
      //Fixing the labels to match the number of partitions:
      hist = gft::AllocIntArray(n);
      for(p = 0; p < nnodes; p++){
	hist[label[p]] += 1;
      }
      mapping = gft::AllocIntArray(n);
      i = 0;
      for(l = 0; l < n; l++){
	if(hist[l] > 0){
	  mapping[l] = i;
	  i++;
	}
      }
      for(p = 0; p < nnodes; p++)
	label[p] = mapping[label[p]];
      gft::FreeIntArray(&mapping);
      gft::FreeIntArray(&hist);
      
      return label;
    }

    
    void FindMerger(sMergeHistory *MH,
		    int s,
		    sGraph *G,
		    int *oift_label,
		    int energy){
      int i,t,w,smin,tmin,wmin;
      MH->energy = energy;
      for(i = 0; i < G->nodes[s].outdegree; i++){
	t = G->nodes[s].adjList[i];
	w = G->nodes[s].Warcs[i];
	if(oift_label[t] == 0 && w == energy){
	  MH->p = G->nodes[s].id;
	  MH->q = G->nodes[t].id;
	  return;
	}
      }
      smin = tmin = NIL;
      wmin = INT_MAX;
      for(s = 0; s < G->nnodes; s++){
	if(oift_label[s] == 0) continue;
	for(i = 0; i < G->nodes[s].outdegree; i++){
	  t = G->nodes[s].adjList[i];
	  w = G->nodes[s].Warcs[i];
	  if(oift_label[t] == 0 && w < wmin){
	    wmin = w;
	    smin = s;
	    tmin = t;
	  }
	}
      }
      MH->p = G->nodes[smin].id;
      MH->q = G->nodes[tmin].id;
    }


    sMergeHistory *DivisiveClusteringByOIFT_2(sGraph *G){
      sMergeHistory *MH = NULL;
      sPQueue32 *Q,*Q_size;
      sClusteringTree *tree;
      int *label, *cost, *oift_label,*DCC_size;
      sGraph *T;
      int n,i,Wmax,lb,p,q,t,l,nnodes,c,energy;
      int S[3];
      c = G->nnodes;
      MH = (sMergeHistory *)malloc(sizeof(sMergeHistory)*(c-1));
      tree = (sClusteringTree *)malloc(sizeof(sClusteringTree)*c*2);
      if(tree == NULL){
	printf("Error: DivisiveClusteringByOIFT_2\n");
	exit(1);
      }
      for(i = 0; i < c*2; i++){
	tree[i].G = NULL;
	tree[i].right = NIL;
	tree[i].left = NIL;
      }
      tree[0].G = gft::Graph::Clone(G);
      tree[0].seed = 0;
      n = 1;

      S[0] = 1;
      S[1] = 0;
      nnodes = G->nnodes;
      label = gft::AllocIntArray(nnodes);
      cost  = gft::AllocIntArray(nnodes);
      DCC_size = gft::AllocIntArray(nnodes);
      for(p = 0; p < nnodes; p++)
	label[p] = NIL;
      label[0] = 0;
      T = Transpose(G);
      ift::IFT_fw(T, S, label, cost);
      Destroy(&T);

      Wmax = GetMaximumArc(G);
      Q = PQueue32::Create(Wmax+2, nnodes, cost);
      Q_size = PQueue32::Create(nnodes+2, nnodes, DCC_size);
      
      for(p = 1; p < G->nnodes; p++)
	PQueue32::FastInsertElem(Q, p);
      i = 0;
      while(!PQueue32::IsEmpty(Q) && i < c-1) {
	energy = PQueue32::FastGetMaxVal(Q);
	
	while(!PQueue32::IsEmpty(Q) &&
	      PQueue32::FastGetMaxVal(Q) == energy){
	  p = PQueue32::FastRemoveMaxFIFO(Q);
	  DCC_size[p] = ComputeDCCsize(G, energy, p);
	  PQueue32::FastInsertElem(Q_size, p);
	  //printf("1) energy: %d, size: %d, cost: %d\n",energy, DCC_size[p], cost[p]);
	}
	
	while(!PQueue32::IsEmpty(Q_size) && i < c-1){
	  p = PQueue32::FastRemoveMaxFIFO(Q_size);

	  lb = label[p];
	  oift_label = gft::AllocIntArray(tree[lb].G->nnodes);
	  T = Transpose(tree[lb].G);
	  S[0] = 2;
	  S[1] = tree[lb].seed;
	  S[2] = GetNodeIndex(tree[lb].G, p);
	  q = tree[lb].G->nodes[S[1]].id;
	  for(t = 0; t < tree[lb].G->nnodes; t++)
	    oift_label[t] = NIL;
	  oift_label[S[1]] = 0;
	  oift_label[S[2]] = 1;
	  ift::OIFT(tree[lb].G, T, S, oift_label);
	  Destroy(&T);
	  
	  //-----------------------------
	  FindMerger(&(MH[i]), S[2], tree[lb].G, oift_label, cost[p]);
	  //MH[i].p = p;
	  //MH[i].q = q;
	  //MH[i].energy = cost[p];
	  //-----------------------------	
	  
	  tree[n].G = Split(&(tree[lb].G), oift_label, 1);
	  tree[n+1].G = tree[lb].G;
	  tree[n].seed   = GetNodeIndex(tree[n].G,   p);
	  tree[n+1].seed = GetNodeIndex(tree[n+1].G, q);
	  tree[lb].G = NULL;
	  tree[lb].left = n;
	  tree[lb].right = n+1;
	  for(t = 0; t < tree[n].G->nnodes; t++)
	    label[tree[n].G->nodes[t].id] = n;
	  for(t = 0; t < tree[n+1].G->nnodes; t++)
	    label[tree[n+1].G->nodes[t].id] = n+1;
	  n += 2;
	  gft::FreeIntArray(&oift_label);
	  
	  i++;
	}
      }
      PQueue32::Destroy(&Q_size);
      PQueue32::Destroy(&Q);
      for(i = 0; i < c*2; i++){
	if(tree[i].G != NULL)
	  gft::Graph::Destroy(&tree[i].G);
      }
      free(tree);
      gft::FreeIntArray(&cost);
      gft::FreeIntArray(&DCC_size);
      gft::FreeIntArray(&label);
      return MH;
    }

    

    sMergeHistory *DivisiveClusteringByOIFT(sGraph *G){
      sMergeHistory *MH = NULL;
      sPQueue32 *Q;
      sClusteringTree *tree;
      int *label, *cost, *oift_label;
      sGraph *T;
      int n,i,Wmax,lb,p,q,t,l,nnodes,c;
      int S[3];
      c = G->nnodes;
      MH = (sMergeHistory *)malloc(sizeof(sMergeHistory)*(c-1));
      tree = (sClusteringTree *)malloc(sizeof(sClusteringTree)*c*2);
      if(tree == NULL){
	printf("Error: DivisiveClusteringByOIFT\n");
	exit(1);
      }
      for(i = 0; i < c*2; i++){
	tree[i].G = NULL;
	tree[i].right = NIL;
	tree[i].left = NIL;
      }
      tree[0].G = gft::Graph::Clone(G);
      tree[0].seed = 0;
      n = 1;

      S[0] = 1;
      S[1] = 0;
      nnodes = G->nnodes;
      label = gft::AllocIntArray(nnodes);
      cost  = gft::AllocIntArray(nnodes);
      for(p = 0; p < nnodes; p++)
	label[p] = NIL;
      label[0] = 0;
      T = Transpose(G);
      ift::IFT_fw(T, S, label, cost);
      Destroy(&T);

      Wmax = GetMaximumArc(G);
      Q = PQueue32::Create(Wmax+2, nnodes, cost);

      for(p = 1; p < G->nnodes; p++)
	PQueue32::FastInsertElem(Q, p);
      i = 0;
      while(!PQueue32::IsEmpty(Q) && i < c-1) {
	p = PQueue32::FastRemoveMaxFIFO(Q);
	lb = label[p];

	oift_label = gft::AllocIntArray(tree[lb].G->nnodes);
	T = Transpose(tree[lb].G);
	S[0] = 2;
	S[1] = tree[lb].seed;
	S[2] = GetNodeIndex(tree[lb].G, p);
	q = tree[lb].G->nodes[S[1]].id;
	for(t = 0; t < tree[lb].G->nnodes; t++)
	  oift_label[t] = NIL;
	oift_label[S[1]] = 0;
	oift_label[S[2]] = 1;
	ift::OIFT(tree[lb].G, T, S, oift_label);
	Destroy(&T);

	//-----------------------------
	FindMerger(&(MH[i]), S[2], tree[lb].G, oift_label, cost[p]);
	//MH[i].p = p;
	//MH[i].q = q;
	//MH[i].energy = cost[p];
	//-----------------------------	
	
	tree[n].G = Split(&(tree[lb].G), oift_label, 1);
	tree[n+1].G = tree[lb].G;
	tree[n].seed   = GetNodeIndex(tree[n].G,   p);
	tree[n+1].seed = GetNodeIndex(tree[n+1].G, q);
	tree[lb].G = NULL;
	tree[lb].left = n;
	tree[lb].right = n+1;
	for(t = 0; t < tree[n].G->nnodes; t++)
	  label[tree[n].G->nodes[t].id] = n;
	for(t = 0; t < tree[n+1].G->nnodes; t++)
	  label[tree[n+1].G->nodes[t].id] = n+1;
	n += 2;
	gft::FreeIntArray(&oift_label);

	i++;
      }
      PQueue32::Destroy(&Q);
      for(i = 0; i < c*2; i++){
	if(tree[i].G != NULL)
	  gft::Graph::Destroy(&tree[i].G);
      }
      free(tree);
      gft::FreeIntArray(&cost);
      gft::FreeIntArray(&label);
      return MH;
    }


    void WriteMergeHistory(sMergeHistory *MH, int n, char *file){
      int i;
      FILE *fp;
      fp = fopen(file, "w");
      if(fp == NULL){
	printf("Cannot create file %s\n", file);
	return;
      }
      fprintf(fp, "%d\n", n);
      for(i = 0; i < n; i++){
	fprintf(fp, "%d %d %d\n", MH[i].p, MH[i].q, MH[i].energy);
      }
      fclose(fp);
    }
   

    typedef struct _Arc{
      int p;
      int q;
    } Arc;
    
    void Union(int *R, int *rp, int rq, int *N){
      if(N[*rp] >= N[rq]){
	N[*rp] = N[*rp] + N[rq];
	R[rq] = *rp;
      }
      else{
	N[rq] = N[rq] + N[*rp];
	R[*rp] = rq;
	*rp = rq;
      }
    }
    
    int Find(int *R, int p){
      if(R[p] == p)
	return p;
      else{
	R[p] = Find(R, R[p]);
	return R[p];
      }
    }

    
    //Felzenszwalb and Huttenlocher 2004:
    int *ClusteringByMST2(sGraph *G, float k){
      int Gmax, b, p, q, i, w, n, s, rp, rq, l;
      float tau_p, tau_q, mint;
      Arc **bucket;
      int *nbucket; //nelements
      int *sbucket; //maximum size
      //--------------------------------
      int *R; //representative elements
      int *N; //number of elements
      int *I; //internal difference
      int *L;
      R = (int *) calloc(G->nnodes, sizeof(int));
      N = (int *) calloc(G->nnodes, sizeof(int));
      I = (int *) calloc(G->nnodes, sizeof(int));
      L = (int *) calloc(G->nnodes, sizeof(int));
      if(R == NULL || N == NULL || I == NULL || L == NULL)
	Error(MSG1, "gft::Graph::ClusteringByMST2");
      //--------------------------------
      //Data structure for counting sort:
      Gmax = GetMaximumArc(G);
      bucket = (Arc **) calloc(Gmax+1, sizeof(Arc *));
      nbucket = (int *) calloc(Gmax+1, sizeof(int));
      sbucket = (int *) calloc(Gmax+1, sizeof(int));
      if(bucket == NULL || nbucket == NULL || sbucket == NULL)
	Error(MSG1, "gft::Graph::ClusteringByMST2");
      for(b = 0; b <= Gmax; b++){
	sbucket[b] = 4;
	bucket[b] = (Arc *) calloc(4, sizeof(Arc));
	if(bucket[b] == NULL)
	  Error(MSG1, "gft::Graph::ClusteringByMST2");
      }
      //--------------------------------
      //Sorting the edges:
      for(p = 0; p < G->nnodes; p++){
	for(i = 0; i < G->nodes[p].outdegree; i++){
	  q = G->nodes[p].adjList[i];
	  w = G->nodes[p].Warcs[i];
	  //We are assuming an undirected graph:
	  if(p > q){
	    n = nbucket[w];
	    s = sbucket[w];
	    if(n == s){
	      sbucket[w] = s*2;
	      bucket[w] = (Arc *)realloc(bucket[w], sizeof(Arc)*s*2);
	      if(bucket[w] == NULL)
		Error(MSG1, "gft::Graph::ClusteringByMST2");
	    }
	    bucket[w][n].p = p;
	    bucket[w][n].q = q;
	    nbucket[w] += 1;
	  }
	}
      }
      //--------------------------------
      for(p = 0; p < G->nnodes; p++){
	R[p] = p;
	N[p] = 1;
	I[p] = 0;
	L[p] = 0;
      }

      //Main loop:
      b = 0;
      while(b <= Gmax){
	//skipping empty buckets:
	while(b <= Gmax && nbucket[b] == 0){
	  b++;
	}
	if(b <= Gmax){
	  nbucket[b] -= 1;
	  n = nbucket[b];
	  p = bucket[b][n].p;
	  q = bucket[b][n].q;
	  w = b;
	  rp = Find(R, p);
	  rq = Find(R, q);
	  if(rp == rq) continue;

	  tau_p = k/N[rp];
	  tau_q = k/N[rq];
	  mint = MIN(I[rp]+tau_p, I[rq]+tau_q);
	  
	  if(w <= mint){
	    I[rp] = w;
	    I[rq] = w;
	    Union(R, &rp, rq, N);
	  }
	}
      }
      l = 0;
      for(p = 0; p < G->nnodes; p++){
	R[p] = Find(R, p);
	if(R[p] == p){
	  L[p] = l;
	  l++;
	}
      }
      for(p = 0; p < G->nnodes; p++){
	L[p] = L[R[p]];
      }
      //--------------------------------
      for(b = 0; b <= Gmax; b++)
	free(bucket[b]);
      free(bucket);
      free(nbucket);
      free(sbucket);
      free(R);
      free(N);
      free(I);
      return L;
    }
    

  } //end Graph namespace
} //end gft namespace

