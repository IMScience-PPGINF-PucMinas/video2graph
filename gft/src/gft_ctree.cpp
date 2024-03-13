
#include "gft_ctree.h"

namespace gft{
  namespace CTree{

    sCTree *Create(int size){
      sCTree *ct = NULL;
      int i;
      ct = (sCTree *)calloc(1,sizeof(sCTree));
      if(ct == NULL)
	gft::Error((char *)MSG1,(char *)"CTree::Create");
      ct->nodes = (sCTreeNode *)calloc(size, sizeof(sCTreeNode));
      if(ct->nodes == NULL)
	gft::Error((char *)MSG1,(char *)"CTree::Create");
      for(i = 0; i < size; i++){
	ct->nodes[i].parent = NULL;
	ct->nodes[i].nsons = 0;
	ct->nodes[i].sons = NULL;
      }
      ct->arraysize = size;
      ct->nnodes = 0;
      ct->nleaves = 0;
      ct->nspels = 0;
      ct->mapping = NULL;
      return ct;
    }
    
    void Destroy(sCTree **ct){
      sCTree *aux;
      int i;
      if(ct != NULL){
	aux = *ct;
	if(aux != NULL){
	  for(i = 0; i < aux->arraysize; i++)
	    if(aux->nodes[i].sons != NULL)
	      free(aux->nodes[i].sons);
	  if(aux->nodes != NULL)
	    free(aux->nodes);
	  if(aux->mapping != NULL)
	    free(aux->mapping);
	  free(aux);
	  *ct = NULL;
	}
      }
    }


    void Print(sCTree *ct){
      int i,j,b;
      //--------------------------------------------
      sCTreeNode *T;
      int k;
      printf("-------PRINT-----------------\n");
      printf("nnodes: %d\n",ct->nnodes);
      for(i = 0; i < ct->nnodes; i++){
	printf("Node #%d: level: %d, height: %d, area: %d, volume: %d",
	       i,ct->nodes[i].level,
	       ct->nodes[i].height,
	       ct->nodes[i].area,
	       ct->nodes[i].volume);
	b = 0;
	for(j = 0; j < ct->nodes[i].nsons; j++){
	  T = ct->nodes[i].sons[j];
	  k = GetNodeIndex(ct, T);
	  if(b)
	    printf(", %d",  k);
	  else{
	    printf(", sons: %d", k);
	    b = 1;
	  }
	}
	printf("\n");
      }
      printf("------------------------------\n");
      //--------------------------------------------
      /*
      printf("nnodes: %d\n",ct->nnodes);
      for(i = 0; i < ct->nnodes; i++){
	printf("Node #%d: level: %d",i,ct->nodes[i].level);
	b = 0;
	for(j = 0; j < ct->nnodes; j++){
	  if(ct->nodes[j].parent == &(ct->nodes[i])){
	    if(b)
	      printf(", %d",j);
	    else{
	      printf(", sons: %d",j);
	      b = 1;
	    }
	  }
	}
	printf("\n");
      }
      */
    }


    int GetNodeIndex(sCTree *ct, sCTreeNode *node){
      int k;
      k = ((long long)node-(long long)&(ct->nodes[0]))/sizeof(sCTreeNode);
      return k;
    }


    sCTree *EdgeBasedMinTree(sGraph *graph, int Lmin){
      sCTree *ct,*ct_fix;
      sCTreeNode *p_ancestor, *q_ancestor;
      sCTreeNode **stack;
      sPQueue32 *Q;
      int nedges, Wmax, p, q, t, w, i, j, k, l, nsons, nunnecessary;
      bool unnecessary;
      int *arc_w, *arc_p, *arc_q, *first_k;
      ct = CTree::Create(graph->nnodes*2);
      ct->nnodes = graph->nnodes;
      
      //printf("nnodes: %d\n", graph->nnodes);

      for(p = 0; p < graph->nnodes; p++){
	ct->nodes[p].level = Lmin; //NIL;
	ct->nodes[p].parent = NULL;
      }
      
      nedges = Graph::GetNumberOfArcs(graph);
      Wmax = Graph::GetMaximumArc(graph);
      arc_w = (int *)malloc(sizeof(int)*nedges);
      arc_p = (int *)malloc(sizeof(int)*nedges);
      arc_q = (int *)malloc(sizeof(int)*nedges);
      first_k = (int *)malloc(sizeof(int)*graph->nnodes);
      stack = (sCTreeNode **)malloc(sizeof(sCTreeNode *)*graph->nnodes);
      if(arc_w == NULL || arc_p == NULL || arc_q == NULL)
	gft::Error((char *)MSG1,(char *)"Graph::EdgeBasedMinTree");
      Q = PQueue32::Create(Wmax + 2, nedges, arc_w);

      k = 0;
      for(p = 0; p < graph->nnodes; p++){
	first_k[p] = k;
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  q = graph->nodes[p].adjList[i];
	  w = graph->nodes[p].Warcs[i];

	  arc_w[k] = w;
	  arc_p[k] = p;
	  arc_q[k] = q;
	  PQueue32::FastInsertElem(Q, k);
	  k++;
	}
      }

      while(!PQueue32::IsEmpty(Q)){
	k = PQueue32::RemoveMinFIFO(Q);
	w = arc_w[k];
	p = arc_p[k];
	q = arc_q[k];

	p_ancestor = &(ct->nodes[p]);
	while(p_ancestor->parent != NULL){
	  t = 0;
	  stack[t] = p_ancestor;
	  t++;
	  p_ancestor = p_ancestor->parent;
	  
	  l = p_ancestor->level;
	  while(p_ancestor->parent != NULL &&
		p_ancestor->parent->level == l){
	    stack[t] = p_ancestor;
	    t++;
	    p_ancestor = p_ancestor->parent;
	  }
	  //path compression:
	  while(t > 0){
	    t--;
	    stack[t]->parent = p_ancestor;
	  }
	}

	q_ancestor = &(ct->nodes[q]);
	while(q_ancestor->parent != NULL){
	  t = 0;
	  stack[t] = q_ancestor;
	  t++;
	  q_ancestor = q_ancestor->parent;
	  
	  l = q_ancestor->level;
	  while(q_ancestor->parent != NULL &&
		q_ancestor->parent->level == l){
	    stack[t] = q_ancestor;
	    t++;
	    q_ancestor = q_ancestor->parent;
	  }
	  //path compression:
	  while(t > 0){
	    t--;
	    stack[t]->parent = q_ancestor;
	  }
	}

	
	if(p_ancestor != q_ancestor){
	  /*
	  if(p == 28 || q == 28){
	    printf("***********************\n");
	    printf("...p_ancestor_level: %d\n",p_ancestor->level);
	    printf("...q_ancestor_level: %d\n",q_ancestor->level);
	    printf("...w: %d\n",w);
	    printf("...p: %d, q: %d\n",p,q);
	    printf("***********************\n");
	  }
	  */
	  if(w > p_ancestor->level && w > q_ancestor->level){
	    ct->nodes[ct->nnodes].level = w;
	    ct->nodes[ct->nnodes].parent = NULL;
	    p_ancestor->parent = &(ct->nodes[ct->nnodes]);
	    q_ancestor->parent = &(ct->nodes[ct->nnodes]);
	    ct->nnodes++;
	  }
	  else if(w == p_ancestor->level && w > q_ancestor->level)
	    q_ancestor->parent = p_ancestor;
	  else if(w > p_ancestor->level && w == q_ancestor->level)
	    p_ancestor->parent = q_ancestor;
	  else{
	    /*
	    printf("Error: unforeseen case.\n");
	    printf("...p_ancestor_level: %d\n",p_ancestor->level);
	    printf("...q_ancestor_level: %d\n",q_ancestor->level);
	    printf("...w: %d\n",w);
	    printf("...p: %d, q: %d\n",p,q);
	    */
	    if(p_ancestor < q_ancestor)
	      q_ancestor->parent = p_ancestor;
	    else
	      p_ancestor->parent = q_ancestor;	      
	  }
	}
	
	//Vizinhos t de p:
	for(i = 0; i < graph->nodes[p].outdegree; i++){
	  if(w != graph->nodes[p].Warcs[i])
	    continue;
	  t = graph->nodes[p].adjList[i];
	  //Obter k de pt:
	  k = first_k[p] + i;
	  if(Q->L.elem[k].color == GRAY){
	    PQueue32::RemoveElem(Q, k);
	    PQueue32::FastInsertElemAsFirst(Q, k);
	  }
	}

	//Vizinhos t de q:
	for(i = 0; i < graph->nodes[q].outdegree; i++){
	  if(w != graph->nodes[q].Warcs[i])
	    continue;
	  t = graph->nodes[q].adjList[i];
	  //Obter k de qt:
	  k = first_k[q] + i;
	  if(Q->L.elem[k].color == GRAY){
	    PQueue32::RemoveElem(Q, k);
	    PQueue32::FastInsertElemAsFirst(Q, k);
	  }
	}
	
      }

      //Final path compression:
      for(p = 0; p < graph->nnodes; p++){
      
        p_ancestor = &(ct->nodes[p]);
	while(p_ancestor->parent != NULL){
	  t = 0;
	  stack[t] = p_ancestor;
	  t++;
	  p_ancestor = p_ancestor->parent;
	  
	  l = p_ancestor->level;
	  while(p_ancestor->parent != NULL &&
		p_ancestor->parent->level == l){
	    stack[t] = p_ancestor;
	    t++;
	    p_ancestor = p_ancestor->parent;
	  }
	  //path compression:
	  while(t > 0){
	    t--;
	    stack[t]->parent = p_ancestor;
	  }
	}

      }

      //--------New code to remove unnecessary nodes---------
      nunnecessary = 0;
      for(p = 0; p < ct->nnodes; p++){
	l = ct->nodes[p].level;
	p_ancestor = ct->nodes[p].parent;
	if(p_ancestor != NULL)
	  if(p_ancestor->level == l)
	    nunnecessary++;
      }
      printf("nunnecessary: %d\n", nunnecessary);

      ct_fix = CTree::Create(ct->nnodes - nunnecessary);
      ct_fix->mapping = (int *)malloc(sizeof(int)*ct->nnodes);
      ct_fix->nnodes = ct->nnodes - nunnecessary;

      q = 0;
      for(p = 0; p < ct->nnodes; p++){
	ct_fix->mapping[p] = NIL;
	l = ct->nodes[p].level;
	p_ancestor = ct->nodes[p].parent;
	unnecessary = false;
	if(p_ancestor != NULL)
	  if(p_ancestor->level == l)
	    unnecessary = true;
	if(!unnecessary){
	  ct_fix->nodes[q] = ct->nodes[p];
	  ct_fix->mapping[p] = q;
	  q++;
	}
      }

      for(k = 0; k < ct_fix->nnodes; k++){
	if(ct_fix->nodes[k].parent != NULL){
	  p = GetNodeIndex(ct, ct_fix->nodes[k].parent);
	  q = ct_fix->mapping[p];
	  ct_fix->nodes[k].parent = &(ct_fix->nodes[q]);
	}
      }

      
      for(p = 0; p < ct->nnodes; p++){
	if(ct_fix->mapping[p] == NIL){
	  q = GetNodeIndex(ct, ct->nodes[p].parent);
	  t = ct_fix->mapping[q];
	  ct_fix->mapping[p] = t;
	}
      }
     
      
      ct_fix->nspels = graph->nnodes;
      ct_fix->mapping = (int *)realloc(ct_fix->mapping,
				       sizeof(int)*graph->nnodes);
      gft::CTree::Destroy(&ct);
      ct = ct_fix;

      //-----------------------------------------------------
      
      //Compute nsons:
      for(p = 0; p < ct->nnodes; p++){
	if(ct->nodes[p].parent != NULL)
	  (ct->nodes[p].parent)->nsons++;
      }
      //Alloc sons:
      ct->nleaves = 0;
      for(p = 0; p < ct->nnodes; p++){
	nsons = ct->nodes[p].nsons;
	if(nsons > 0)
	  ct->nodes[p].sons = (sCTreeNode **)malloc(nsons*sizeof(sCTreeNode *));
	else{
	  ct->nodes[p].sons = NULL;
	  ct->nleaves++;
	}
      }
      //Set sons:
      for(p = 0; p < ct->nnodes; p++)
	ct->nodes[p].nsons = 0;
      for(p = 0; p < ct->nnodes; p++){
	if(ct->nodes[p].parent != NULL){
	  nsons = (ct->nodes[p].parent)->nsons;
	  (ct->nodes[p].parent)->sons[nsons] = &(ct->nodes[p]);
	  (ct->nodes[p].parent)->nsons++;
	}
      }

      PQueue32::Destroy(&Q);
      free(arc_w);
      free(arc_p);
      free(arc_q);
      free(first_k);
      free(stack);
      return ct;
    }



    sCTree *EdgeBasedMinTree(sImageGraph *graph, int Lmin){
      sCTree *ct;
      sGraph *G;
      G = Graph::Clone(graph);
      ct = EdgeBasedMinTree(G, Lmin);
      Graph::Destroy(&G);
      return ct;
    }

    

    void ComputeHeight(sCTree *ct){
      int i,j,p,q,dl,root,h,hmax;
      sCTreeNode *T;
      gft::sStack *S;
      gft::sQueue *Q;
      Q = gft::Queue::Create(ct->nnodes);
      S = gft::Stack::Create(ct->nnodes);
      for(i = 0; i < ct->nnodes; i++)
	if(ct->nodes[i].parent == NULL)
	  root = i;

      //Reverse Breadth-first search:
      gft::Queue::Push(Q, root);
      while(!gft::Queue::IsEmpty(Q)){
	p = gft::Queue::Pop(Q);
	gft::Stack::Push(S, p);

	for(j = 0; j < ct->nodes[p].nsons; j++){
	  T = ct->nodes[p].sons[j];
	  q = GetNodeIndex(ct, T);
	  gft::Queue::Push(Q, q);
	}
      }
      
      while(!gft::Stack::IsEmpty(S)){
	p = gft::Stack::Pop(S);

	if(ct->nodes[p].parent != NULL)
	  dl = abs(ct->nodes[p].level - (ct->nodes[p].parent)->level);
	else
	  dl = 0;
	  
	hmax = dl;
	for(j = 0; j < ct->nodes[p].nsons; j++){
	  T = ct->nodes[p].sons[j];
	  h = T->height + dl;
	  if(h > hmax)
	    hmax = h;
	}
	ct->nodes[p].height = hmax;
      }
      gft::Stack::Destroy(&S);
      gft::Queue::Destroy(&Q);
    }


    void ComputeArea(sCTree *ct, int *area){
      int i,j,p,q,l,root,area_sum;
      sCTreeNode *T;
      gft::sStack *S;
      gft::sQueue *Q;
      bool flag = false;

      if(area == NULL){
	flag = true;
	area = (int *)calloc(ct->nnodes, sizeof(int));
	for(j = 0; j < ct->nspels; j++){
	  i = ct->mapping[j];
	  //printf("i: %d\n", i);
	  area[i]++;
	}
      }

      Q = gft::Queue::Create(ct->nnodes);
      S = gft::Stack::Create(ct->nnodes);
      for(i = 0; i < ct->nnodes; i++)
	if(ct->nodes[i].parent == NULL)
	  root = i;

      //Reverse Breadth-first search:
      gft::Queue::Push(Q, root);
      while(!gft::Queue::IsEmpty(Q)){
	p = gft::Queue::Pop(Q);
	gft::Stack::Push(S, p);

	for(j = 0; j < ct->nodes[p].nsons; j++){
	  T = ct->nodes[p].sons[j];
	  q = GetNodeIndex(ct, T);
	  gft::Queue::Push(Q, q);
	}
      }
      
      while(!gft::Stack::IsEmpty(S)){
	p = gft::Stack::Pop(S);
	if(ct->nodes[p].nsons == 0)
	  ct->nodes[p].area = area[p];
	else{
	  l = ct->nodes[p].level;
	  area_sum = 0;
	  for(j = 0; j < ct->nodes[p].nsons; j++){
	    T = ct->nodes[p].sons[j];
	    area_sum += T->area;
	  }
	  ct->nodes[p].area = area_sum;
	}
      }
      if(flag) free(area);
      gft::Stack::Destroy(&S);
      gft::Queue::Destroy(&Q);
    }


    //You must call 'ComputeArea' first.
    void ComputeVolume(sCTree *ct){
      int i,j,p,q,root,dl,vsum;
      sCTreeNode *T;
      gft::sStack *S;
      gft::sQueue *Q;
      Q = gft::Queue::Create(ct->nnodes);
      S = gft::Stack::Create(ct->nnodes);
      for(i = 0; i < ct->nnodes; i++)
	if(ct->nodes[i].parent == NULL)
	  root = i;

      //Reverse Breadth-first search:
      gft::Queue::Push(Q, root);
      while(!gft::Queue::IsEmpty(Q)){
	p = gft::Queue::Pop(Q);
	gft::Stack::Push(S, p);

	for(j = 0; j < ct->nodes[p].nsons; j++){
	  T = ct->nodes[p].sons[j];
	  q = GetNodeIndex(ct, T);
	  gft::Queue::Push(Q, q);
	}
      }
      
      while(!gft::Stack::IsEmpty(S)){
	p = gft::Stack::Pop(S);

	if(ct->nodes[p].parent != NULL)
	  dl = abs(ct->nodes[p].level - (ct->nodes[p].parent)->level);
	else
	  dl = 0;
	  
	vsum = dl*ct->nodes[p].area;
	for(j = 0; j < ct->nodes[p].nsons; j++){
	  T = ct->nodes[p].sons[j];
	  vsum += T->volume;
	}
	ct->nodes[p].volume = vsum;
      }
      gft::Stack::Destroy(&S);
      gft::Queue::Destroy(&Q);
    }
    


    int *ComputeExtinctionValue(sCTree *ct, AttributeType type){
      int *ext, *attr;
      int l,a,c,i,k,extinction;
      bool continue_loop;
      bool *visited;
      sCTreeNode *node_p, *node_a, *node_c;
      ext = (int *)malloc(ct->nleaves*sizeof(int));
      visited = (bool *)malloc(ct->nnodes*sizeof(bool));
      attr = (int *)malloc(ct->nnodes*sizeof(int));

      for(i = 0; i < ct->nnodes; i++){
	visited[i] = false;
	switch(type){
	case height:
	  attr[i] = ct->nodes[i].height;
	  break;
	case area:
	  attr[i] = ct->nodes[i].area;
	  break;
	case volume:
	  attr[i] = ct->nodes[i].volume;
	  break;
	}
      }

      for(l = 0; l < ct->nleaves; l++){
	extinction = INT_MAX;
	continue_loop = true;
	node_p = &(ct->nodes[l]);
	while(continue_loop && node_p != NULL){
	  node_a = node_p;
	  a = GetNodeIndex(ct, node_a);
	  node_p = node_a->parent;

	  if(node_p != NULL && node_p->nsons > 1){
	    for(k = 0; k < node_p->nsons && continue_loop; k++){
	      node_c = node_p->sons[k];
	      c = GetNodeIndex(ct, node_c);
	      
	      if((visited[c] == true &&
		  node_c != node_a &&
		  attr[c] == attr[a])||
		 (node_c != node_a && attr[c] > attr[a]))
		continue_loop = false;

	      visited[c] = true;
	    }
	  }
	  if(node_p != NULL)
	    extinction = attr[a];
	  ext[l] = extinction;
	}
      }
      free(attr);
      free(visited);
      return ext;
    }


    
  }
}
