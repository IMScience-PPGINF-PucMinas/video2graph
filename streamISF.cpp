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
#include <sstream>
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

// GRAFO TEMPORAL

int getLatestEdgeValue() {
    std::string edges_path = "./graph/A.txt";
    std::ifstream edges_file(edges_path);

    int edge_value = 0;

    std::string penultimate_line;
    std::string latest_line;
    while (std::getline(edges_file, latest_line)) {
        if (!latest_line.empty()) {
            penultimate_line = latest_line;
        }
    }

    edges_file.close();

    std::stringstream ss(penultimate_line);
    std::string edge_value_str;
    std::getline(ss, edge_value_str, ',');
    try{
        edge_value = std::stoi(edge_value_str);
    }
    catch (const std::invalid_argument& e){
        edge_value = 0;
    }

    return edge_value;

}

int getLatestFrameIndicator() {
    std::string frame_indicator_path = "./graph/frame_indicator.txt";

    std::ifstream frame_indicator_file(frame_indicator_path);

    if (!frame_indicator_file.is_open()) {
        return 1;
    }

    int latest_frame = 0;
    std::string line;

    while (std::getline(frame_indicator_file, line)) {
        latest_frame = std::stoi(line);
    }

    frame_indicator_file.close();

    return latest_frame + 1;
}

int getLatestGraphIndicator() {
    std::string graph_indicator_path = "./graph/graph_indicator.txt";

    std::ifstream graph_indicator_file(graph_indicator_path);

    if (!graph_indicator_file.is_open()) {
        return 1;
    }

    int latest_graph = 0;
    std::string line;

    while (std::getline(graph_indicator_file, line)) {
        latest_graph = std::stoi(line);
    }

    graph_indicator_file.close();

    return latest_graph + 1;
}

void writeLabel(int label)
{
    std::string graph_labels_path = "./graph/graph_labels.txt";
    std::ofstream graph_labels(graph_labels_path, std::ios_base::app);

    graph_labels << label << "\n";
    graph_labels.close();
}

void writeGraph(const std::vector<std::vector<int>> &adjMat, const std::unordered_map<int, int> &vals, int edge_n, int frame, int video)
{
    std::string indicator_path = "./graph/graph_indicator.txt";
    std::string frame_path = "./graph/frame_indicator.txt";
    std::string edges_path = "./graph/A.txt";
    std::string edge_attributes_path = "./graph/edge_attributes.txt";
    std::string node_path = "./graph/node_labels.txt";

    std::ofstream indicator(indicator_path, std::ios_base::app);
    std::ofstream frames(frame_path, std::ios_base::app);
    std::ofstream edges(edges_path, std::ios_base::app);
    std::ofstream edge_attributes(edge_attributes_path, std::ios_base::app);
    std::ofstream node(node_path, std::ios_base::app);

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

std::pair<std::vector<std::vector<int>>, std::unordered_map<int, int>> getAdj(Image* img, int rows, int cols) {
    std::unordered_map<int, int> vals;
    int size = 0;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int pixel = img->val[i][j];
            if (vals.find(pixel) == vals.end()) {
                vals[pixel] = size++;
            }
        }
    }

    std::vector<std::vector<int>> adjMat(size, std::vector<int>(size, 0));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int y = img->val[i][j];
            int reducedValue = vals[y];
            
            if (i - 1 >= 0) { // Left
                int value = img->val[i - 1][j];
                if (y != value) {
                    if (vals.find(value) == vals.end()) {
                        vals[value] = size++;
                    }
                    adjMat[reducedValue][vals[value]] = 1;
                }
            }
            if (j - 1 >= 0) { // Above
                int value = img->val[i][j - 1];
                if (y != value) {
                    if (vals.find(value) == vals.end()) {
                        vals[value] = size++;
                    }
                    adjMat[reducedValue][vals[value]] = 1;
                }
            }
            if (j + 1 < cols) { // Below
                int value = img->val[i][j + 1];
                if (y != value) {
                    if (vals.find(value) == vals.end()) {
                        vals[value] = size++;
                    }
                    adjMat[reducedValue][vals[value]] = 1;
                }
            }
            if (i + 1 < rows) { // Right
                int value = img->val[i + 1][j];
                if (y != value) {
                    if (vals.find(value) == vals.end()) {
                        vals[value] = size++;
                    }
                    adjMat[reducedValue][vals[value]] = 1;
                }
            }
        }
    }

    return {adjMat, vals};
}

void writeForOneFrame(Image **img, int rows, int cols, int video)
{
    std::vector<std::vector<int>> adjMat;
    std::unordered_map<int, int> vals;
    std::tie(adjMat, vals) = getAdj(img[0], rows, cols);
    int edge_n = getLatestEdgeValue();
    int frame = getLatestFrameIndicator();
    writeGraph(adjMat, vals, edge_n, frame, video);
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
    int label = 0;

    int num_video = getLatestGraphIndicator();

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
    //strcpy(outputPath, argv[2]);
    p = atoi(argv[2]);
    nSVX = atoi(argv[3]);
    imgPorBloco = atoi(argv[4]);
    imgPorIntersecao = atoi(argv[5]);
    typeOfColoring = atoi(argv[6]);
    propagatedMoreLabels = atoi(argv[7]);
    rateStopDecrement = (float)atoi(argv[8]) / 100.0;
    label = atoi(argv[9]);

    
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

    num_frames -= 1;

    if (num_frames == 0)
    {
        fprintf(stderr, "Unable to find video frames at %s\n", filename);
        return 0;
    }
    imgPorBloco = (int) ((((float) imgPorBloco)) * num_frames)/100;

    // printf("Imagens por bloco: %d\n", imgPorBloco);
    // printf("Imagens por Interseção: %d\n", imgPorIntersecao);
    // printf("Total number of frames in fold is %d\n", num_frames);

    // Teto do número de frames dividido pela quantidade de imagens por bloco
    quantDeBlocos = ceil((double)num_frames / imgPorBloco);
    // Copia do número de frames para calcular quantos frames ainda faltam, para verificar se no bloco será o imgPorBloco, ou a quantidade restante que é menor que imgPorBloco
    int numFramesRestantes = num_frames;

    // Número do frame, conforme seu nome 00001, 00002, 00003...
    int imgFrames = 1;
    // Quantidade de frames no bloco, respeitando o número de frames, caso seja inferior a quantidade de imagens por bloco
    int frames = 0;

    sumTotalTime = 0.00;

    // Número do frame atual, pois para gravar, criar na pasta números diferentes para guardar todos os frames, sem gravar um em cima do outro
    nFrame = 0;
    initial_label = 0;

    // Matriz que guardará a respectiva cor de cada label para ao gravar imagem segmentada as cores de um bloco para outro serem as mesmas
    int numMaxTree = quantDeBlocos * nSVX;
    //int numMaxTree = nSVX;
    int** randowColorLabels = new int*[numMaxTree];
    for (int i = 0; i < numMaxTree; i++) {
        randowColorLabels[i] = new int[3];
        randowColorLabels[i][0] = -1;
        randowColorLabels[i][1] = -1;
        randowColorLabels[i][2] = -1;
    }

    // Estratégia para conseguir um número aproximado de supervoxels no vídeo todo de forma que o número de supervoxels
    // por bloco vai decrescendo
    // aux = (int) (((float) nSVX) / quantDeBlocos);
    if (propagatedMoreLabels == 1) {
        aux = (int) (((((float) nSVX)) * quantDeBlocos) / 100.0);
        // aux2 é o limite para não decrementar muito e ficar com muito pouco supervoxel no final
        if (nSVX != 0) {
            aux2 = (int) ( ((float) nSVX * rateStopDecrement));
        }
        sum = 0;
        k = nSVX;
    }

    Intersection* intersection;

    for (int bloco = 0; bloco < quantDeBlocos; bloco++)
    {
        // Calcular quantos frames terá no bloco que será processado
        // printf("Número de frames restantes = %d\n", numFramesRestantes);
        if (numFramesRestantes > imgPorBloco)
        {
            frames = imgPorBloco;
            numFramesRestantes -= imgPorBloco;
        }
        else
        {
            frames = numFramesRestantes;
            numFramesRestantes = 0;
        }
        if (bloco != 0) {
            frames = frames + imgPorIntersecao;
        }
        
        // printf("Quantidade de frames no bloco %d = %d\n", (bloco + 1), frames);

        // Concatenate Images in a video structure
        video = (Image **)malloc(frames * sizeof(Image *));
        Image **video_border = (Image **)malloc(frames * sizeof(Image *));
        sprintf(filepath, "%s/frame1.ppm", filename);
        img = ReadAnyImage(filepath);

        int i = 0;
        // O bloco terá as imagens da interseção juntamente com a quantidade de imagens estabelecidas por bloco
        if (imgFrames != 1) {
            imgFrames = imgFrames - imgPorIntersecao;
        }
        
        // printf("Comeca na imagem: %d   -   Frames no bloco: %d \n", imgFrames, frames);

        // Calcular a quantidade de supervoxels nos blocos
        if (propagatedMoreLabels == 1) {
            sum = aux + sum;
            if (k > aux2) {
                k = nSVX - sum;
            }
        } else {
            k = (nSVX - nSVXcomputed) / (quantDeBlocos - bloco);
        }

        // Escreve o label da classe pra cada video
        writeLabel(label);

        // Ler e armazenar as imagens do bloco atual na variavel video
        while (i < frames && imgFrames <= num_frames)
        {
            frame_id = imgFrames;
            
            sprintf(filepath, "%s/frame%d.ppm", filename, frame_id);
            //sprintf(filepath, "%s/frame%d.ppm", filename, frame_id);
            imgD = loadImage(filepath);

            //printf("%05d.ppm", frame_id);

            border_img = createImage(imgD->num_rows, imgD->num_cols, 1);

            video[i] = imgD;
            video_border[i] = border_img;
        
            graph = createGraph(video, 1, i);
            label_video = runDISF(graph, p, k, 5, 1, 1, 1, NULL, 0);
            
            writeForOneFrame(label_video, video[i]->num_rows, video[i]->num_cols, num_video);
            
            free(label_video);

            i++;
            imgFrames++;

        }

    }

    return 0;
}