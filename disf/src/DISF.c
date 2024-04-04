#include "DISF.h"
#include "Color.h"
#include "Image.h"
#include "IntList.h"
#include "PrioQueue.h"
#include "Utils.h"
#include <math.h>

//=============================================================================
// Bool Functions
//=============================================================================
inline bool areValidNodeCoords(Graph *graph, NodeCoords coords)
{
    return (coords.x >= 0 && coords.x < graph->num_cols) &&
           (coords.y >= 0 && coords.y < graph->num_rows) &&
           (coords.z >= 0 && coords.z < graph->num_frames);
}

//=============================================================================
// Int Functions
//=============================================================================
inline int getNodeIndex(Graph *graph, NodeCoords coords)
{
    return (coords.y * graph->num_cols + coords.x) + (coords.z * graph->num_rows * graph->num_cols);
}

//=============================================================================
// Double Functions
//=============================================================================
inline double euclDistance(float *feat1, float *feat2, int num_feats)
{
    double dist;

    dist = 0;

    for (int i = 0; i < num_feats; i++)
        dist += (feat1[i] - feat2[i]) * (feat1[i] - feat2[i]);
    dist = sqrtf(dist);

    return dist;
}

inline double taxicabDistance(float *feat1, float *feat2, int num_feats)
{
    double dist;

    dist = 0;

    for (int i = 0; i < num_feats; i++)
        dist += fabs(feat1[i] - feat2[i]);

    return dist;
}

//=============================================================================
// NodeCoords Functions
//=============================================================================
inline NodeCoords getAdjacentNodeCoords(NodeAdj *adj_rel, NodeCoords coords, int id)
{
    NodeCoords adj_coords;

    adj_coords.x = coords.x + adj_rel->dx[id];
    adj_coords.y = coords.y + adj_rel->dy[id];
    adj_coords.z = coords.z + adj_rel->dz[id];

    return adj_coords;
}

inline NodeCoords getNodeCoords(Graph *graph, int index)
{
    NodeCoords coords;

    int np = (graph->num_cols * graph->num_rows);
    coords.z = floor(index / np);
    int tmp = index % np;
    coords.y = floor(tmp / graph->num_cols);
    coords.x = tmp % graph->num_cols;

    return coords;
}

//=============================================================================
// Float* Functions
//=============================================================================
inline float *meanTreeFeatVector(Tree *tree)
{
    float *mean_feat;

    mean_feat = (float *)calloc(tree->num_feats, sizeof(float));

    for (int i = 0; i < tree->num_feats; i++)
        mean_feat[i] = tree->sum_feat[i] / (float)tree->num_nodes;

    return mean_feat;
}

//=============================================================================
// Double* Functions
//=============================================================================
double *computeGradient(Graph *graph)
{
    float max_adj_dist, sum_weight;
    float *dist_weight;
    double *grad;
    NodeAdj *adj_rel;

    grad = (double *)calloc(graph->num_nodes, sizeof(double));
    // adj_rel = create8NeighAdj();
    adj_rel = create6NeighAdj();

    max_adj_dist = sqrtf(2); // Diagonal distance for 8-neighborhood
    dist_weight = (float *)calloc(adj_rel->size, sizeof(float));
    sum_weight = 0;

    // Compute the inverse distance weights (closer --> higher; farther --> lower)
    for (int i = 0; i < adj_rel->size; i++)
    {
        float div;

        // Distance between the adjacent and the center
        div = sqrtf(adj_rel->dx[i] * adj_rel->dx[i] + adj_rel->dy[i] * adj_rel->dy[i] + adj_rel->dz[i] * adj_rel->dz[i]);

        dist_weight[i] = max_adj_dist / div;
        sum_weight += dist_weight[i];
    }

    // Normalize values
    for (int i = 0; i < adj_rel->size; i++)
        dist_weight[i] /= sum_weight;

// Compute the gradients
#pragma omp parallel for
    for (int i = 0; i < graph->num_nodes; i++)
    {
        float *feats;
        NodeCoords coords;

        feats = graph->feats[i];
        coords = getNodeCoords(graph, i);

        // For each adjacent node
        for (int j = 0; j < adj_rel->size; j++)
        {
            float *adj_feats;
            NodeCoords adj_coords;

            adj_coords = getAdjacentNodeCoords(adj_rel, coords, j);

            // Is Valid?
            if (areValidNodeCoords(graph, adj_coords))
            {
                int adj_index;
                double dist;

                adj_index = getNodeIndex(graph, adj_coords);

                adj_feats = graph->feats[adj_index];

                // Compute L1 Norm between center and adjacent
                dist = taxicabDistance(adj_feats, feats, graph->num_feats);

                // Weight by its distance relevance
                grad[i] += dist * dist_weight[j];
            }
        }
    }

    Image **video = createVideo(graph->num_rows, graph->num_cols, 1, graph->num_frames);
    // Insert the lowest gradient positions as seeds, but avoiding
    // repetition
    int k = 0;
    for (int f = 0; f < graph->num_frames; f++)
    {
        char concat_label[256];
        int id = f + 1;
        sprintf(concat_label, "%s%d%s", "GA_", id, ".pgm");
        for (int i = 0; i < graph->num_rows * graph->num_cols; i++)
        {
            video[f]->val[i][0] = (int)grad[k];
            k++;
        }
        // writeImagePGM(video[f], concat_label);
        // freeImage(&video[f]);
    }

    free(dist_weight);
    freeNodeAdj(&adj_rel);

    return grad;
}

//=============================================================================
// NodeAdj* Functions
//=============================================================================
NodeAdj *create4NeighAdj()
{
    NodeAdj *adj_rel;

    adj_rel = (NodeAdj *)calloc(1, sizeof(NodeAdj));

    adj_rel->size = 4;
    adj_rel->dx = (int *)calloc(4, sizeof(int));
    adj_rel->dy = (int *)calloc(4, sizeof(int));

    adj_rel->dx[0] = -1;
    adj_rel->dy[0] = 0; // Left
    adj_rel->dx[1] = 1;
    adj_rel->dy[1] = 0; // Right

    adj_rel->dx[2] = 0;
    adj_rel->dy[2] = -1; // Top
    adj_rel->dx[3] = 0;
    adj_rel->dy[3] = 1; // Bottom

    return adj_rel;
}

NodeAdj *create8NeighAdj()
{
    NodeAdj *adj_rel;

    adj_rel = (NodeAdj *)calloc(1, sizeof(NodeAdj));

    adj_rel->size = 8;
    adj_rel->dx = (int *)calloc(8, sizeof(int));
    adj_rel->dy = (int *)calloc(8, sizeof(int));

    adj_rel->dx[0] = -1;
    adj_rel->dy[0] = 0; // Center-Left
    adj_rel->dx[1] = 1;
    adj_rel->dy[1] = 0; // Center-Right

    adj_rel->dx[2] = 0;
    adj_rel->dy[2] = -1; // Top-Center
    adj_rel->dx[3] = 0;
    adj_rel->dy[3] = 1; // Bottom-Center

    adj_rel->dx[4] = -1;
    adj_rel->dy[4] = 1; // Bottom-Left
    adj_rel->dx[5] = 1;
    adj_rel->dy[5] = -1; // Top-Right

    adj_rel->dx[6] = -1;
    adj_rel->dy[6] = -1; // Top-Left
    adj_rel->dx[7] = 1;
    adj_rel->dy[7] = 1; // Bottom-Right

    return adj_rel;
}

NodeAdj *create6NeighAdj()
{
    NodeAdj *adj_rel;

    adj_rel = (NodeAdj *)calloc(1, sizeof(NodeAdj));

    adj_rel->size = 6;
    adj_rel->dx = (int *)calloc(6, sizeof(int));
    adj_rel->dy = (int *)calloc(6, sizeof(int));
    adj_rel->dz = (int *)calloc(6, sizeof(int));
    // Pixel frame =================================================================
    adj_rel->dx[0] = -1;
    adj_rel->dy[0] = 0;
    adj_rel->dz[0] = 0; // Left
    adj_rel->dx[1] = 1;
    adj_rel->dy[1] = 0;
    adj_rel->dz[0] = 0; // Right

    adj_rel->dx[2] = 0;
    adj_rel->dy[2] = -1;
    adj_rel->dz[0] = 0; // Top
    adj_rel->dx[3] = 0;
    adj_rel->dy[3] = 1;
    adj_rel->dz[0] = 0; // Bottom
                        // Previous frame ==============================================================
    adj_rel->dx[4] = 0;
    adj_rel->dy[4] = 0;
    adj_rel->dz[4] = -1; // Center-back
                         // Posterior frame ==============================================================
    adj_rel->dx[5] = 0;
    adj_rel->dy[5] = 0;
    adj_rel->dz[5] = 1; // Center-front

    return adj_rel;
}

NodeAdj *create18NeighAdj()
{
    NodeAdj *adj_rel;

    adj_rel = (NodeAdj *)calloc(1, sizeof(NodeAdj));

    adj_rel->size = 18;
    adj_rel->dx = (int *)calloc(18, sizeof(int));
    adj_rel->dy = (int *)calloc(18, sizeof(int));
    adj_rel->dz = (int *)calloc(18, sizeof(int));
    // Pixel frame =================================================================
    adj_rel->dx[0] = -1;
    adj_rel->dy[0] = 0;
    adj_rel->dz[0] = 0; // Center-Left
    adj_rel->dx[1] = 1;
    adj_rel->dy[1] = 0;
    adj_rel->dz[1] = 0; // Center-Right

    adj_rel->dx[2] = 0;
    adj_rel->dy[2] = -1;
    adj_rel->dz[2] = 0; // Top-Center
    adj_rel->dx[3] = 0;
    adj_rel->dy[3] = 1;
    adj_rel->dz[3] = 0; // Bottom-Center

    adj_rel->dx[4] = -1;
    adj_rel->dy[4] = 1;
    adj_rel->dz[4] = 0; // Bottom-Left
    adj_rel->dx[5] = 1;
    adj_rel->dy[5] = -1;
    adj_rel->dz[5] = 0; // Top-Right

    adj_rel->dx[6] = -1;
    adj_rel->dy[6] = -1;
    adj_rel->dz[6] = 0; // Top-Left
    adj_rel->dx[7] = 1;
    adj_rel->dy[7] = 1;
    adj_rel->dz[7] = 0; // Bottom-Right

    // Previous frame ==============================================================
    adj_rel->dx[8] = -1;
    adj_rel->dy[8] = 0;
    adj_rel->dz[8] = -1; // Center-Left
    adj_rel->dx[9] = 1;
    adj_rel->dy[9] = 0;
    adj_rel->dz[9] = -1; // Center-Right

    adj_rel->dx[10] = 0;
    adj_rel->dy[10] = -1;
    adj_rel->dz[10] = -1; // Top-Center
    adj_rel->dx[11] = 0;
    adj_rel->dy[11] = 1;
    adj_rel->dz[11] = -1; // Bottom-Center

    adj_rel->dx[12] = 0;
    adj_rel->dy[12] = 0;
    adj_rel->dz[12] = -1; // Center-back

    // Posterior frame ==============================================================
    adj_rel->dx[13] = -1;
    adj_rel->dy[13] = 0;
    adj_rel->dz[13] = 1; // Center-Left
    adj_rel->dx[14] = 1;
    adj_rel->dy[14] = 0;
    adj_rel->dz[14] = 1; // Center-Right

    adj_rel->dx[15] = 0;
    adj_rel->dy[15] = -1;
    adj_rel->dz[15] = 1; // Top-Center
    adj_rel->dx[16] = 0;
    adj_rel->dy[16] = 1;
    adj_rel->dz[16] = 1; // Bottom-Center

    adj_rel->dx[17] = 0;
    adj_rel->dy[17] = 0;
    adj_rel->dz[17] = 1; // Center-front

    return adj_rel;
}

NodeAdj *create26NeighAdj()
{
    NodeAdj *adj_rel;

    adj_rel = (NodeAdj *)calloc(1, sizeof(NodeAdj));

    adj_rel->size = 26;
    adj_rel->dx = (int *)calloc(26, sizeof(int));
    adj_rel->dy = (int *)calloc(26, sizeof(int));
    adj_rel->dz = (int *)calloc(26, sizeof(int));
    // Pixel frame =================================================================
    adj_rel->dx[0] = -1;
    adj_rel->dy[0] = 0;
    adj_rel->dz[0] = 0; // Center-Left
    adj_rel->dx[1] = 1;
    adj_rel->dy[1] = 0;
    adj_rel->dz[1] = 0; // Center-Right

    adj_rel->dx[2] = 0;
    adj_rel->dy[2] = -1;
    adj_rel->dz[2] = 0; // Top-Center
    adj_rel->dx[3] = 0;
    adj_rel->dy[3] = 1;
    adj_rel->dz[3] = 0; // Bottom-Center

    adj_rel->dx[4] = -1;
    adj_rel->dy[4] = 1;
    adj_rel->dz[4] = 0; // Bottom-Left
    adj_rel->dx[5] = 1;
    adj_rel->dy[5] = -1;
    adj_rel->dz[5] = 0; // Top-Right

    adj_rel->dx[6] = -1;
    adj_rel->dy[6] = -1;
    adj_rel->dz[6] = 0; // Top-Left
    adj_rel->dx[7] = 1;
    adj_rel->dy[7] = 1;
    adj_rel->dz[7] = 0; // Bottom-Right

    // Previous frame ==============================================================
    adj_rel->dx[8] = -1;
    adj_rel->dy[8] = 0;
    adj_rel->dz[8] = -1; // Center-Left
    adj_rel->dx[9] = 1;
    adj_rel->dy[9] = 0;
    adj_rel->dz[9] = -1; // Center-Right

    adj_rel->dx[10] = 0;
    adj_rel->dy[10] = -1;
    adj_rel->dz[10] = -1; // Top-Center
    adj_rel->dx[11] = 0;
    adj_rel->dy[11] = 1;
    adj_rel->dz[11] = -1; // Bottom-Center

    adj_rel->dx[12] = -1;
    adj_rel->dy[12] = 1;
    adj_rel->dz[12] = -1; // Bottom-Left
    adj_rel->dx[13] = 1;
    adj_rel->dy[13] = -1;
    adj_rel->dz[13] = -1; // Top-Right

    adj_rel->dx[14] = -1;
    adj_rel->dy[14] = -1;
    adj_rel->dz[14] = -1; // Top-Left
    adj_rel->dx[15] = 1;
    adj_rel->dy[15] = 1;
    adj_rel->dz[15] = -1; // Bottom-Right

    adj_rel->dx[16] = 0;
    adj_rel->dy[16] = 0;
    adj_rel->dz[16] = -1; // Center-back

    // Posterior frame ==============================================================
    adj_rel->dx[17] = -1;
    adj_rel->dy[17] = 0;
    adj_rel->dz[17] = 1; // Center-Left
    adj_rel->dx[18] = 1;
    adj_rel->dy[18] = 0;
    adj_rel->dz[18] = 1; // Center-Right

    adj_rel->dx[19] = 0;
    adj_rel->dy[19] = -1;
    adj_rel->dz[19] = 1; // Top-Center
    adj_rel->dx[20] = 0;
    adj_rel->dy[20] = 1;
    adj_rel->dz[20] = 1; // Bottom-Center

    adj_rel->dx[21] = -1;
    adj_rel->dy[21] = 1;
    adj_rel->dz[21] = 1; // Bottom-Left
    adj_rel->dx[22] = 1;
    adj_rel->dy[22] = -1;
    adj_rel->dz[22] = 1; // Top-Right

    adj_rel->dx[23] = -1;
    adj_rel->dy[23] = -1;
    adj_rel->dz[23] = 1; // Top-Left
    adj_rel->dx[24] = 1;
    adj_rel->dy[24] = 1;
    adj_rel->dz[24] = 1; // Bottom-Right

    adj_rel->dx[25] = 0;
    adj_rel->dy[25] = 0;
    adj_rel->dz[25] = 1; // Center-front

    return adj_rel;
}
//=============================================================================
// Graph* Functions
//=============================================================================
Graph *createGraph(Image **video, int num_frames, int start_frame)
{
    // printf("========Creating Graph from video========\n");
    int normval;
    Graph *graph;
    Image *img = video[start_frame]; // first frame
    normval = getNormValue(img);
    graph = (Graph *)calloc(1, sizeof(Graph));

    graph->num_cols = img->num_cols;
    graph->num_rows = img->num_rows;
    graph->num_frames = num_frames;
    graph->num_feats = 3; // L*a*b cspace
    graph->num_nodes = img->num_pixels * num_frames;

    // printf("Graph information ---\n");
    // printf("cols:%d row:%d\n", graph->num_cols, graph->num_rows);
    // printf("frames:%d nodes:%d\n", graph->num_frames, graph->num_nodes);
    // printf("feats:%d\n", graph->num_feats);

    graph->feats = (float **)calloc(graph->num_nodes, sizeof(float *));
    int count = 0;

    for (int f = 0; f < graph->num_frames; ++f)
    {
        //  #pragma omp parallel for
        for (int i = 0; i < img->num_pixels; ++i)
        {
            if (video[f + start_frame]->num_channels == 1) // Grayscale
                graph->feats[count] = convertGrayToLab(video[f + start_frame]->val[i], normval);
            else // sRGB
                graph->feats[count] = convertsRGBToLab(video[f + start_frame]->val[i], normval);
            count++;
        }
    }
    // printf("Final count index :%d\n", count);
    return graph;
}

//=============================================================================
// Tree* Functions
//=============================================================================
Tree *createTree(int root_index, int num_feats)
{
    Tree *tree;

    tree = (Tree *)calloc(1, sizeof(Tree));

    tree->root_index = root_index;
    tree->num_nodes = 0;
    tree->num_feats = num_feats;

    tree->sum_feat = (float *)calloc(num_feats, sizeof(float));

    return tree;
}

//=============================================================================
// Image* Functions
//=============================================================================
Image **runDISF(Graph *graph, int n_0, int n_f, int adj_op, int sampl_op, int path_op, int rem_op, Image ***border_vid, int initial_label)
{
    // printf("========Running DISF Algorithm========\n");
    bool want_borders;
    int num_rem_seeds, iter;
    int *pred_map, *rnd_indexes;
    double *cost_map;
    NodeAdj *adj_rel;
    IntList *seed_set;
    int num_pixels = graph->num_cols * graph->num_rows;
    // Image *label_img;
    Image **label_video;
    PrioQueue *queue;
    NodeCoords coords;

    // Iniciar estruturas auxiliares
    cost_map = (double *)calloc(graph->num_nodes, sizeof(double));
    pred_map = (int *)calloc(graph->num_nodes, sizeof(int));

    if (adj_op == 1)
        adj_rel = create4NeighAdj();
    else if (adj_op == 2)
        adj_rel = create8NeighAdj();
    else if (adj_op == 3)
        adj_rel = create6NeighAdj();
    else if (adj_op == 4)
        adj_rel = create18NeighAdj();
    else
        adj_rel = create26NeighAdj();

    // label_img = createImage(graph->num_rows, graph->num_cols, 1);
    label_video = createVideo(graph->num_rows, graph->num_cols, 1, graph->num_frames);
    queue = createPrioQueue(graph->num_nodes, cost_map, MINVAL_POLICY);

    want_borders = border_vid != NULL;

    if (sampl_op == 1)
    {
        seed_set = createIntList();
        rnd_indexes = randomSampling(0, graph->num_nodes - 1, n_0);
        for (int i = 0; i < n_0; ++i)
            insertIntListTail(&seed_set, rnd_indexes[i]);
        free(rnd_indexes);
    }
    else
        seed_set = gridSampling(graph, n_0);

    // printf("NUM SEEDS: %d\n", seed_set->size);
    // printf("n_0: %d\n", n_0);

    iter = 1; // Pelo menos uma única iteração é executada
    do
    {
        int seed_label, num_trees, num_maintain, sedex;
        Tree **trees;
        IntList **tree_adj;
        bool **are_trees_adj;
        trees = (Tree **)calloc(seed_set->size, sizeof(Tree *));
        tree_adj = (IntList **)calloc(seed_set->size, sizeof(IntList *));
        are_trees_adj = (bool **)calloc(seed_set->size, sizeof(bool *));

        // Assign initial values for all nodes

        // #pragma omp parallel for
        int c = 0;
        // printf("Numero de frames no for do DISF: %d\n", graph->num_frames);
        // printf("number of pixels %d\n", num_pixels );
        for (int f = 0; f < graph->num_frames; f++)
        {
            for (int i = 0; i < num_pixels; i++)
            {
                cost_map[c] = INFINITY;
                pred_map[c] = -1;

                label_video[f]->val[i][0] = -1;
                c++;

                if (want_borders)
                    (*border_vid)[f]->val[i][0] = 0;
            }
        }
        // printf("Values Assigned for all nodes------>\n");
        // for(int f = 0; f < graph->num_frames; f++){
        //  Atribuir valores iniciais para todas as sementes amostradas
        //  DANI - Cada semente da grid é uma árvore, pois recebe um seed_label diferente
        seed_label = initial_label;
        // printf("-------------seed_label: %d\n", seed_label);

        for (IntCell *ptr = seed_set->head; ptr != NULL; ptr = ptr->next)
        {
            // printf("To aqui 1,5\n");
            int seed_index;

            seed_index = ptr->elem;
            coords = getNodeCoords(graph, seed_index);
            sedex = (coords.y * graph->num_cols + coords.x);

            // printf("Seed Index = %d, z coord= %d\n", seed_index, coords.z);
            cost_map[seed_index] = 0;
            label_video[coords.z]->val[sedex][0] = seed_label;

            trees[seed_label] = createTree(seed_index, graph->num_feats);
            tree_adj[seed_label] = createIntList();
            are_trees_adj[seed_label] = (bool *)calloc(seed_set->size, sizeof(bool));

            seed_label++;
            insertPrioQueue(&queue, seed_index);
        }
        //}
        // printf("Values Assigned for all seed sampled <------\n");

        // for(int f = 0; f < graph->num_frames; f++){
        //  Para cada nó dentro da fila
        while (!isPrioQueueEmpty(queue))
        {
            int nodex, node_index, node_label, addex, root_index;
            NodeCoords node_coords, A_coord;
            float *mean_feat_tree;

            node_index = popPrioQueue(&queue);
            node_coords = getNodeCoords(graph, node_index);
            nodex = (node_coords.y * graph->num_cols + node_coords.x);
            node_label = label_video[node_coords.z]->val[nodex][0];

            // Inserimos os recursos na respectiva árvore neste
            // momento, pois é garantido que este nó não
            // será inserido nunca mais.
            insertNodeInTree(graph, node_index, &(trees[node_label]));
            root_index = trees[node_label]->root_index;

            // propósitos de excesso de velocidade
            mean_feat_tree = meanTreeFeatVector(trees[node_label]);

            // For each adjacent node
            for (int i = 0; i < adj_rel->size; i++)
            {
                NodeCoords adj_coords;

                adj_coords = getAdjacentNodeCoords(adj_rel, node_coords, i);

                // Is valid?
                if (areValidNodeCoords(graph, adj_coords))
                {
                    int adj_index, adj_label;

                    adj_index = getNodeIndex(graph, adj_coords);
                    // printf("ADJ index = %d\n", adj_index );

                    A_coord = getNodeCoords(graph, adj_index);
                    addex = (A_coord.y * graph->num_cols + A_coord.x);

                    adj_label = label_video[A_coord.z]->val[addex][0];

                    // This adjacent was already added to a tree?
                    if (queue->state[adj_index] != BLACK_STATE)
                    {
                        double arc_cost, path_cost;

                        if (path_op == 1)
                        {
                            arc_cost = euclDistance(graph->feats[root_index], graph->feats[adj_index], graph->num_feats);
                        }
                        else
                        {
                            arc_cost = euclDistance(mean_feat_tree, graph->feats[adj_index], graph->num_feats);
                        }
                        // printf("arc_cost: %d\n", arc_cost);

                        path_cost = MAX(cost_map[node_index], arc_cost);

                        // Can this node be conquered by the current tree?
                        if (path_cost < cost_map[adj_index])
                        {
                            cost_map[adj_index] = path_cost;
                            label_video[A_coord.z]->val[addex][0] = node_label;
                            pred_map[adj_index] = node_index;

                            // Update if it is already in the queue
                            if (queue->state[adj_index] == GRAY_STATE)
                                moveIndexUpPrioQueue(&queue, adj_index);
                            else
                                insertPrioQueue(&queue, adj_index);
                        }
                    }
                    else if (node_label != adj_label) // If they differ, their trees are adjacent
                    {

                        if (want_borders) // Both depicts a border between their superpixels
                        {
                            // printf("Want border = %d\n", want_borders );
                            (*border_vid)[A_coord.z]->val[nodex][0] = 255;
                            (*border_vid)[A_coord.z]->val[addex][0] = 255;
                        }

                        // Were they defined as adjacents?
                        if (!are_trees_adj[node_label][adj_label])
                        {
                            insertIntListTail(&(tree_adj[node_label]), adj_label);
                            insertIntListTail(&(tree_adj[adj_label]), node_label);
                            are_trees_adj[adj_label][node_label] = true;
                            are_trees_adj[node_label][adj_label] = true;
                        }
                    }
                }
            }
        }

        //}
        // Compute the number of seeds to be preserved
        num_maintain = MAX(n_0 * exp(-iter), n_f);

        // Auxiliar var
        num_trees = seed_set->size;

        if (rem_op == 1)
        {
            IntList *new_seed_set = createIntList();
            rnd_indexes = randomSampling(0, seed_set->size - 1, num_maintain);
            for (int i = 0; i < num_maintain; ++i)
            {
                int seed_index = getIntListElem(seed_set, rnd_indexes[i]);
                insertIntListTail(&new_seed_set, seed_index);
            }
            freeIntList(&seed_set);
            free(rnd_indexes);
            seed_set = new_seed_set;
        }
        else
        {
            freeIntList(&seed_set);
            // Select the most relevant superpixels
            // Aqui que é calculado o V(x) e é retirado as árvores de menos relevancia da fila de prioridade
            seed_set = selectKMostRelevantSeeds(trees, tree_adj, graph->num_nodes, num_trees, num_maintain);
        }

        // Compute the number of seeds to be removed
        num_rem_seeds = num_trees - seed_set->size;

        iter++;
        resetPrioQueue(&queue); // Clear the queue

        // printf("Queue processedf \n");
        //
    } while (num_rem_seeds > 0);

    // printf("Number of Iterations: %d\n", iter);

    return label_video;
}

//=============================================================================
// IntList* Functions
//=============================================================================
IntList *gridSampling(Graph *graph, int num_seeds)
{
    // printf("Creating grid ------>\n");
    float size, stride, delta_x, delta_y, delta_z;
    double *grad;
    bool *is_seed;
    IntList *seed_set;
    NodeAdj *adj_rel;
    Image **video;

    seed_set = createIntList();
    is_seed = (bool *)calloc(graph->num_nodes, sizeof(bool));
    // printf("gridSampling graph->num_nodes: %d\n", graph->num_nodes);
    //  Compute the approximate superpixel size and stride
    size = 0.5 + (float)(graph->num_nodes / (float)num_seeds);
    stride = pow(size, 1.0 / 3) + 0.5;
    // printf("stride gridSampling: %lf\n", stride);

    delta_x = delta_y = delta_z = stride / 2.0;
    // printf("delta_z gridSampling: %lf\n", delta_z);
    // printf("delta x %f, y %f, z %f\n", delta_x, delta_y, delta_z );

    if (delta_x < 1.0 || delta_y < 1.0 || delta_z < 1.0)
        printError("gridSampling", "The number of samples is too high");

    grad = computeGradient(graph);
    adj_rel = create6NeighAdj();

    // Iterate through the nodes coordinates
    for (int z = (int)delta_z; z < graph->num_frames; z += stride)
    {
        // printf("Entrei no for do gridSampling\n");
        for (int y = (int)delta_y; y < graph->num_rows; y += stride)
        {
            for (int x = (int)delta_x; x < graph->num_cols; x += stride)
            {
                int min_grad_index;
                NodeCoords curr_coords;

                curr_coords.x = x;
                curr_coords.y = y;
                curr_coords.z = z;

                min_grad_index = getNodeIndex(graph, curr_coords);

                // For each adjacent node
                for (int i = 0; i < adj_rel->size; i++)
                {
                    NodeCoords adj_coords;

                    adj_coords = getAdjacentNodeCoords(adj_rel, curr_coords, i);

                    // Is valid?
                    if (areValidNodeCoords(graph, adj_coords))
                    {
                        int adj_index;

                        adj_index = getNodeIndex(graph, adj_coords);

                        // The gradient in the adjacent is minimum?
                        if (grad[adj_index] < grad[min_grad_index])
                            min_grad_index = adj_index;
                    }
                }

                // Select the position with lowest gradient
                is_seed[min_grad_index] = true;
            }
        }
    }

    video = createVideo(graph->num_rows, graph->num_cols, 1, graph->num_frames);
    // Insert the lowest gradient positions as seeds, but avoiding
    // repetition
    int k = 0;
    for (int f = 0; f < graph->num_frames; f++)
    {
        char concat_label[256];
        int id = f + 1;
        sprintf(concat_label, "%s%d%s", "grid_", id, ".pgm");
        for (int i = 0; i < graph->num_rows * graph->num_cols; i++)
        {
            if (is_seed[k])
            {
                // printf("Entrei no if (is_seed[k])\n");
                insertIntListTail(&seed_set, k);
                video[f]->val[i][0] = 255;
            }
            k++;
        }
        //  writeImagePGM(video[f], concat_label);
        //  freeImage(&video[f]);
    }

    // printf("Finished grid <------\n");

    free(grad);
    free(is_seed);
    freeNodeAdj(&adj_rel);

    return seed_set;
}

IntList *selectKMostRelevantSeeds(Tree **trees, IntList **tree_adj, int num_nodes, int num_trees, int num_maintain)
{
    double *tree_prio;
    IntList *rel_seeds;
    PrioQueue *queue;

    tree_prio = (double *)calloc(num_trees, sizeof(double));
    rel_seeds = createIntList();
    queue = createPrioQueue(num_trees, tree_prio, MAXVAL_POLICY);

// For each tree
#pragma omp parallel for
    for (int i = 0; i < num_trees; i++)
    {
        double area_prio, grad_prio;
        float *mean_feat_i;

        // Compute the area relevance
        // |Tx| / |N|
        area_prio = trees[i]->num_nodes / (float)num_nodes;

        // Initial values for the computation of gradient relevance
        grad_prio = INFINITY;
        mean_feat_i = meanTreeFeatVector(trees[i]); // Speeding purposes

        // For each adjacent tree
        for (IntCell *ptr = tree_adj[i]->head; ptr != NULL; ptr = ptr->next)
        {
            int adj_tree_id;
            float *mean_feat_j;
            double dist;

            adj_tree_id = ptr->elem;
            mean_feat_j = meanTreeFeatVector(trees[adj_tree_id]);

            // Compute the L2 norm between trees
            dist = euclDistance(mean_feat_i, mean_feat_j, trees[i]->num_feats);

            // Get the minimum gradient value
            grad_prio = MIN(grad_prio, dist);
        }

        // Compute the superpixel relevance
        tree_prio[i] = area_prio * grad_prio;

#pragma omp critical
        insertPrioQueue(&queue, i);
    }

    // While it is possible to get relevant seeds
    for (int i = 0; i < num_maintain && !isPrioQueueEmpty(queue); i++)
    {
        int tree_id, root_index;

        tree_id = popPrioQueue(&queue);
        root_index = trees[tree_id]->root_index;

        insertIntListTail(&rel_seeds, root_index);
    }
    // The rest is discarted

    freePrioQueue(&queue);
    free(tree_prio);

    return rel_seeds;
}

//=============================================================================
// Void Functions
//=============================================================================
void freeNodeAdj(NodeAdj **adj_rel)
{
    if (*adj_rel != NULL)
    {
        NodeAdj *tmp;

        tmp = *adj_rel;

        free(tmp->dx);
        free(tmp->dy);
        free(tmp);

        *adj_rel = NULL;
    }
}

void freeTree(Tree **tree)
{
    if (*tree != NULL)
    {
        Tree *tmp;

        tmp = *tree;

        free(tmp->sum_feat);
        free(tmp);

        *tree = NULL;
    }
}

void freeGraph(Graph **graph)
{
    if (*graph != NULL)
    {
        Graph *tmp;

        tmp = *graph;

        for (int i = 0; i < tmp->num_nodes; i++)
            free(tmp->feats[i]);
        free(tmp->feats);
        free(tmp);

        *graph = NULL;
    }
}

void insertNodeInTree(Graph *graph, int index, Tree **tree)
{
    (*tree)->num_nodes++;

    for (int i = 0; i < graph->num_feats; i++)
        (*tree)->sum_feat[i] += graph->feats[index][i];
}