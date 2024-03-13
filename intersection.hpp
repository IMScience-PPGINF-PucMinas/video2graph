#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

#include <iostream>
#include "Image.h"

//=============================================================================
// Structures
//=============================================================================
/**
* @brief    Estrutura que armazena informações da interseção e do bloco que terá
*           a propagação das árvores
*/
typedef struct
{
    int amountLabels;                       /* Quantidade total de labels */
    int amountTreeIntersection;             /* Intervalo da quantidade de árvores na imagem de interseção */
    int biggestIntersectionLabel;           /* Maior valor de label da árvore da imagem da interseção */
    int smallestIntersectionLabel;          /* Menor valor de label da árvore da imagem da interseção */
    int* labels;                            /* Labels alterados após propagação */
    int* amountVertices;                    /* Quantidade de vértices nas árvores no bloco que está recebendo a propagação */
    int* amountVerticesIntersection;        /* Quantidade de vértices nas árvores da interseção */
    int** depthTreesBlock;                  /* Profundidade das árvores em cada bloco */
    int** positionLabelsBlock;              /* Matrix que em cada linha é referente a um label e possui a posição de onde os vértices deste label estão na imagem */
    int** positionLabelsIntersection;       /* Matrix que contém a posição dos vértices de cada label da interseção */
} Intersection;

//=============================================================================
// Functions
//=============================================================================
Intersection* newIntersection(int size);
Intersection* changeLabelValue(int maxLabel, Intersection* intersection);
int insertVerticePositionBlock (int line, int position, Intersection* intersection);
int insertVerticePositionIntersection (int line, int position, Intersection* intersection);
int chooseLabelWithMoreVertices(Intersection* intersection, int amountVerticesTree, int max);
int getMaxNumLabel(Intersection* intersection);
void deleteIntersection(Intersection* intersection);
void getDepthTreesBlock (Intersection* intersection, Image **imgBlock, int num_pixels, int numBlock, int num_frames);
void propagateIntersectingTrees (Intersection* intersection, Image **intersectionImage, int propagatedMoreLabels);
void insertIntersection (Image **intersectionImage, Image **intersectionBlock, int num_pixels, Intersection* intersection);
void printPositionVerticesIntersection (Intersection* intersection);
void printPositionVerticesBlock (Intersection* intersection);
void printLabels (Intersection* intersection);

#endif