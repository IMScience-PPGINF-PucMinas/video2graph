
#include "gft.h"
#include "intersection.hpp"

#include <string>
#include <dirent.h>
#include <stdio.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <set>
#include <map>

#include "Image.h"
#include "DISF.h"
#include "Utils.h"

#include <fstream>
#include <string>
#include <random>

#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

gft::sImage32 *ReadAnyImage(char *file)
{
    gft::sImage32 *img;
    char command[512];
    int s;

    s = strlen(file);
    if (strcasecmp(&file[s - 3], "pgm") == 0)
    {
        img = gft::Image32::Read(file);
    }
    else
    {
        sprintf(command, "convert %s image_tmp.pgm", file);
        system(command);
        img = gft::Image32::Read("image_tmp.pgm");
        system("rm image_tmp.pgm");
    }
    return img;
}

float **Alocar_matriz(int m, int n)
{
    float **v;
    int i;
    if (m < 1 || n < 1)
    { /* verifica parametros recebidos */
        printf("** Erro: Parametro invalido**\n");
        return (NULL);
    }
    /* aloca as linhas da matriz */
    v = (float **)calloc(m, sizeof(float *));
    if (v == NULL)
    {
        printf("** Erro: Memoria Insuficiente **");
        return (NULL);
    }
    /* aloca as colunas da matriz */
    for (i = 0; i < m; i++)
    {
        v[i] = (float *)calloc(n, sizeof(float));
        if (v[i] == NULL)
        {
            printf("** Erro: Memoria Insuficiente **");
            return (NULL);
        }
    }
    return (v); /* retorna o ponteiro para a matriz */
}

float **Liberar_matriz(int m, int n, float **v)
{
    int i;
    if (v == NULL)
        return (NULL);
    if (m < 1 || n < 1)
    {
        printf("** Erro: Parametro invalido**\n");
        return (v);
    }
    for (i = 0; i < m; i++)
    {
        if (v[i] != NULL)
        {
            free(v[i]); /* libera as linhas da matriz */
            v[i] = NULL;
        }
    }

    free(v);       /* libera a matriz */
    return (NULL); /* retorna um ponteiro nulo */
}

void Liberar_Video(Image **video, int frames, int num_pixels)
{

    for (int j = 0; j < frames; j++)
    {

        for (int i = 0; i < num_pixels; i++)
        {
            if (video[j]->val[i])
            {
                free(video[j]->val[i]);
            }
        }
        if (video[j]->val)
        {
            free(video[j]->val);
        }
        if (video[j])
        {
            free(video[j]);
        }
    }
    free(video);
}

int GetMaxValVoxel(gft::sImage32 **scnLabel, int num_frames)
{
    int i, max, n, p = 0;
    // printf("scnLabel[p]->ncols = %d \t scnLabel[p]->nrows = %d\n", scnLabel[p]->ncols, scnLabel[p]->nrows);
    // printf("NUM_FRAMES = %d\n", num_frames);
    n = scnLabel[p]->ncols * scnLabel[p]->nrows;
    max = scnLabel[p]->data[0];
    for (p = 0; p < num_frames; p++)
        for (i = 1; i < n; i++)
            if (scnLabel[p]->data[i] > max)
                max = scnLabel[p]->data[i];

    return (max);
}

float **getAverageColor(gft::sImage32 **scnLabel, Image **video, int num_frames, int Imax, int **randowColorLabels, int typeOfColoring)
{
    int nnodes, i, k, l, n;
    // float *value[3];
    int *size;
    float **value;
    n = scnLabel[0]->ncols * scnLabel[0]->nrows;
    value = Alocar_matriz(3, (Imax + 1));
    nnodes = GetMaxValVoxel(scnLabel, num_frames) + 1;
    size = gft::AllocIntArray(nnodes);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 255);

    if (value == NULL || size == NULL)
    {
        printf("Error\n");
        exit(1);
    }

    // printf("num_frames: %d\n", num_frames);
    for (i = 0; i < n; i++)
    {
        for (k = 0; k < num_frames; k++)
        {
            l = scnLabel[k]->data[i];
            size[l] += 1;
            value[0][l] += video[k]->val[i][0];
            value[1][l] += video[k]->val[i][1];
            value[2][l] += video[k]->val[i][2];
        }
    }
    // printf("nnodes: %d   -  n: %d    -    Imax: %d\n", nnodes, n, Imax);

    for (l = 0; l < nnodes; l++)
    {
        if (size[l] > 0)
        {
            if (randowColorLabels[l][0] == -1)
            {
                randowColorLabels[l][0] = dis(gen);
                randowColorLabels[l][1] = dis(gen);
                randowColorLabels[l][2] = dis(gen);
            }

            if (typeOfColoring == 0)
            {
                value[0][l] /= size[l];
                value[1][l] /= size[l];
                value[2][l] /= size[l];
            }
            else
            {
                value[0][l] = randowColorLabels[l][0];
                value[1][l] = randowColorLabels[l][1];
                value[2][l] = randowColorLabels[l][2];
            }
        } // else {
          // printf("Empty regions in RAG: %d  -  l: %d\n", size[l], l);
        //}
    }
    // gft::FreeIntArray(&size);

    return value;
}

int computedSupervoxelsVideo(int **randowColorLabels, int nSVX)
{
    int nSupervoxelsVideo = 0;

    for (int i = 0; i < nSVX; i++)
    {
        if (randowColorLabels[i][0] != -1)
        {
            nSupervoxelsVideo++;
        }
    }

    return nSupervoxelsVideo;
}

Image **copyTreeIntersection(Image **original, Image **intersecao, int nImagesIntersecao, int num_rows, int num_cols, int num_frames, Intersection *intersection)
{
    int num_pixels = num_cols * num_rows;
    int nframeOriginal = 0;
    int position = -1;

    // Reservo espaco na memoria para vetor com as imagens da intersecao
    intersecao = (Image **)malloc(nImagesIntersecao * sizeof(Image *));

    // Reservar espaço na memoria para cada imagem da intersecao
    for (int f = 0; f < nImagesIntersecao; f++)
    {
        intersecao[f] = createImage(num_rows, num_cols, 1);
    }

    // Copiar label das arvores da intersecao do bloco anterior
    for (int f = 0; f < nImagesIntersecao; f++)
    {
        // Copiar apenas os frames que desejo copiar o label no final do grafo original
        nframeOriginal = num_frames - nImagesIntersecao + f;
        for (int i = 0; i < num_pixels; i++)
        {
            position = original[nframeOriginal]->val[i][0] % intersection->amountLabels;
            intersecao[f]->val[i][0] = intersection->labels[position];
        }
    }
    return intersecao;
}

void Liberar_Intersection(Image **intersecao, int nImagesIntersecao, int num_pixels)
{

    if (intersecao != nullptr)
    {
        for (int j = 0; j < nImagesIntersecao; j++)
        {

            for (int i = 0; i < num_pixels; i++)
            {
                if (intersecao[j]->val[i] != nullptr)
                {
                    free(intersecao[j]->val[i]);
                }
            }
            if (intersecao[j]->val != nullptr)
            {
                free(intersecao[j]->val);
            }
            if (intersecao[j] != nullptr)
            {
                free(intersecao[j]);
            }
        }
        free(intersecao);
    }
}

void printIntersecao(Image **intersecao, int nImagesIntersecao, int num_rows, int num_cols)
{
    int num_pixels = num_cols * num_rows;
    if (intersecao)
    {
        for (int f = 0; f < nImagesIntersecao; f++)
        {
            for (int i = 0; i < num_pixels; i++)
            {
                printf("val = %d \t - ", intersecao[f]->val[i][0]);
            }
        }
    }
}

void gravarArvoresEmArquivo(Image **original, Image **intersecao, int num_rows, int num_cols, int num_frames, int nImagesIntersecao, std::string fileNameBloco, std::string fileNameIntersecao, std::string debug, int bloco, Intersection *intersection)
{
    int num_pixels = num_cols * num_rows;
    int aux = 1;
    int positionLabel = -1;
    std::ofstream intersecaoFile(fileNameIntersecao, std::ios::app);

    intersecaoFile << debug << " - bloco: " << bloco << "\n";
    if (intersecao)
    {
        for (int f = 0; f < nImagesIntersecao; f++)
        {
            aux = 1;
            for (int i = 0; i < num_pixels; i++)
            {
                if (intersecao[f]->val[i][0] < 10)
                {
                    intersecaoFile << intersecao[f]->val[i][0] << "  ";
                }
                else
                {
                    intersecaoFile << intersecao[f]->val[i][0] << " ";
                }

                if ((aux % num_rows) == 0)
                {
                    intersecaoFile << "\n";
                }
                aux++;
            }
            intersecaoFile << "\n";
        }
    }
    intersecaoFile.close();

    std::ofstream originalFile(fileNameBloco, std::ios::app);

    originalFile << debug << " - bloco: " << bloco << "\n";
    for (int f = 0; f < num_frames; f++)
    {
        aux = 1;
        for (int i = 0; i < num_pixels; i++)
        {
            positionLabel = original[f]->val[i][0] % intersection->amountLabels;
            if (intersection->labels[positionLabel] < 10)
            {
                originalFile << intersection->labels[positionLabel] << "  ";
            }
            else
            {
                originalFile << intersection->labels[positionLabel] << " ";
            }

            if ((aux % num_rows) == 0)
            {
                originalFile << "\n";
            }
            aux++;
        }
        originalFile << "\n";
    }
    originalFile.close();
}

// GRAFO TEMPORAL

void writeLabel(int label)
{
    std::string graph_labels_path = "./graph/graph_labels.txt";
    std::ofstream graph_labels(graph_labels_path, std::ios_base::app);

    graph_labels << label << "\n";
    graph_labels.close();
}

int edge_n = 0;

void writeGraph(const std::vector<std::vector<int>> &adjMat, const std::unordered_map<int, int> &vals, int video, int frame)
{
    std::string indicator_path = "./graph/graph_indicator.txt";
    std::string frame_path = "./graph/frame_indicator.txt";
    std::string edges_path = "./graph/A.txt";
    // std::string graph_labels_path = "./graph/graph_labels.txt";
    std::string edge_attributes_path = "./graph/edge_attributes.txt";
    std::string node_path = "./graph/node_labels.txt";

    // std::ofstream graph_labels(graph_labels_path, std::ios_base::app);
    std::ofstream indicator(indicator_path, std::ios_base::app);
    std::ofstream frames(frame_path, std::ios_base::app);
    std::ofstream edges(edges_path, std::ios_base::app);
    std::ofstream edge_attributes(edge_attributes_path, std::ios_base::app);
    std::ofstream node(node_path, std::ios_base::app);

    // graph_labels << label << "\n";
    // graph_labels.close();

    for (std::vector<int>::size_type i = 0; i < adjMat.size(); ++i)
    {
        indicator << video << "\n";
        frames << frame << "\n";
        for (std::vector<int>::size_type j = 0; j < adjMat[i].size(); ++j)
        {
            if (adjMat[i][j] > 0)
            {
                edges << i + 1 + edge_n << ", " << j + 1 + edge_n << "\n";
                edge_attributes << 0 << "\n";
            }
        }
    }
    edge_n += vals.size();
    indicator.close();
    frames.close();
    edges.close();
    edge_attributes.close();
}

std::pair<std::vector<std::vector<int>>, std::unordered_map<int, int>> getAdj(Image *img, int rows, int cols)
{
    std::unordered_map<int, int> vals;
    int size = 0;
    int chn_begin, chn_end;

    chn_begin = 0;
    chn_end = img->num_channels - 1;

    float **matriz;
    matriz = Alocar_matriz(rows, cols);

    for (int i = 0; i < img->num_pixels; i++)
    {
        int pixel = img->val[i][2];
        // printf("pixel value: %d\n", pixel);
        matriz[(int)(i / rows)][(int)(i % rows)] = pixel;
        if (vals.find(pixel) == vals.end())
        {
            vals[pixel] = size++;
        }
    }

    std::vector<std::vector<int>> adjMat(size, std::vector<int>(size, 0));

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            int y = matriz[i][j];
            int reducedValue = vals[y];

            if (i - 1 >= 0)
            { // Left
                int value = matriz[i - 1][j];
                if (y != value)
                {
                    if (vals.find(value) == vals.end())
                    {
                        vals[value] = size++;
                    }
                    adjMat[reducedValue][vals[value]] = 1;
                }
            }
            if (j - 1 >= 0)
            { // Above
                int value = matriz[i][j - 1];
                if (y != value)
                {
                    if (vals.find(value) == vals.end())
                    {
                        vals[value] = size++;
                    }
                    adjMat[reducedValue][vals[value]] = 1;
                }
            }
            if (j + 1 < cols)
            { // Below
                int value = matriz[i][j + 1];
                if (y != value)
                {
                    if (vals.find(value) == vals.end())
                    {
                        vals[value] = size++;
                    }
                    adjMat[reducedValue][vals[value]] = 1;
                }
            }
            if (i + 1 < rows)
            { // Right
                int value = matriz[i + 1][j];
                if (y != value)
                {
                    if (vals.find(value) == vals.end())
                    {
                        vals[value] = size++;
                    }
                    adjMat[reducedValue][vals[value]] = 1;
                }
            }
        }
    }
    Liberar_matriz(rows, cols, matriz);
    return {adjMat, vals};
}

void writeForOneFrame(Image **img, int rows, int cols, int video, int frame)
{
    std::vector<std::vector<int>> adjMat;
    std::unordered_map<int, int> vals;
    std::tie(adjMat, vals) = getAdj(img[0], rows, cols);
    writeGraph(adjMat, vals, video, frame);
}

int main(int argc, char **argv)
{
    gft::sCImage *ctmp;
    gft::sImage32 *img = NULL;
    DIR *pDir;
    struct dirent *pDirent;
    char filename[512], outputPath[512], filepath[1024], filed[2048];
    int propagatedMoreLabels, aux, aux2, sum, nSVX, nSVXcomputed, k, p, j, sp, num_frames, len, frame_id, x, y, z, n, Imax, i, imgPorBloco, imgPorIntersecao, quantDeBlocos, typeOfColoring, nFrame, num_pixels, initial_label, folder;
    clock_t end, start;
    double totaltime, sumTotalTime;
    Image *imgD, *border_img, **label_video, **video, **intersecao;
    Graph *graph, **graphPerBlock;
    struct stat st;
    float **value;
    char concat_ppm[2048];
    gft::sImage32 **scnLabel;
    float rateStopDecrement;

    int num_video = 1;
    int frame = 1;
    int label = 1;

    if (argc < 7)
    {
        fprintf(stdout, "usage:\n");
        fprintf(stdout, "ex ./isf2svx ./datasets/walking/person01_walking_d2_uncomp/frame ./datasets/walking/person01_walking_d2_uncomp/frameResult 1000 50 12 1 1 0 50\n");
        fprintf(stdout, "./isf2svx <videoFolder> <videoOutput>  <p> <nSVX> <imgPorBloco> <imgPorIntersecao> <typeOfColoring> <propagatedMoreLabels> <rateStopDecrement>\n");
        fprintf(stdout, "\t videoFolder: PPM\n");
        fprintf(stdout, "\t p....... Initial number of superpixels DISF\n");
        fprintf(stdout, "\t nSVX....... The desired number of regions\n");
        fprintf(stdout, "\t imgPorBloco....... Images per block\n");
        fprintf(stdout, "\t imgPorIntersecao....... Images by intersection\n");
        fprintf(stdout, "\t typeOfColoring....... Coloring type when saving segmented image. 0: Medium color - 1: Random color\n");
        fprintf(stdout, "\t propagatedMoreLabels....... Do tree merges, 1 if yes, 0 if no\n");
        fprintf(stdout, "\t rateStopDecrement....... Rate to stop supervoxel decrement. Ex.: 10, will be 0.10 of the number of supervoxels\n");
        // ./streamISF ./dataset/testePPM output/testePPM 15 4 3 1 0 1 50
        exit(0);
    }

    // Guardar dados de argv em suas variáveis correspondentes
    strcpy(filename, argv[1]);
    strcpy(outputPath, argv[2]);
    p = atoi(argv[3]);
    nSVX = atoi(argv[4]);
    imgPorBloco = atoi(argv[5]);
    imgPorIntersecao = atoi(argv[6]);
    typeOfColoring = atoi(argv[7]);
    propagatedMoreLabels = atoi(argv[8]);
    rateStopDecrement = (float)atoi(argv[9]) / 100.0;

    nSVXcomputed = 0, k = 0;

    // Verificar se as imagens lidas possuem extensão .ppm e contabilizar o número de frames na pasta
    num_frames = 0;
    pDir = opendir(filename);

    if (pDir != NULL)
    {
        while ((pDirent = readdir(pDir)) != NULL)
        {
            len = strlen(pDirent->d_name);
            if (len >= 4)
            {
                if (strcmp(".ppm", &(pDirent->d_name[len - 4])) == 0)
                    num_frames++;
            }
        }
    }

    if (num_frames == 0)
    {
        fprintf(stderr, "Unable to find video frames at %s\n", filename);
        return 0;
    }

    for (int n = 0; n < num_frames; n++)
    {
        sprintf(filepath, "%s/%05d.ppm", filename, n + 1);
        imgD = loadImage(filepath);

        printf("%05d.ppm", n);

        video[n] = imgD;

        graph = createGraph(video, 1, n);
        label_video = runDISF(graph, p, k, 5, 1, 1, 1, NULL, 0);

        // aqui irei forçar que o resultado do DISF seja o teste que criamos
        label_video = video;
        // conferindo se é do teste
        writeImagePPM(label_video[0], "TESTE.ppm");

        writeForOneFrame(label_video, video[n]->num_rows, video[n]->num_cols, num_video, frame);
        frame++;

        // free(label_video);
    }

    return 0;
}