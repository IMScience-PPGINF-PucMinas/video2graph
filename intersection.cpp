#include "intersection.hpp"

/**
 * @brief   Alocar espaço em memória para estruturas que serão utilizadas para 
 *          propagar árvores de interseção para o próximo bloco
 * 
 * @param   size Quantidade de labels no bloco, é a quantidade de supervoxels por bloco
 * 
 * @details Alocar espaço em memória para estruturas que serão utilizadas para 
 *          propagar árvores de interseção para o próximo bloco
 * 
 * @return intersection - Estrutura que será alocada
*/
Intersection* newIntersection(int size) {
    Intersection* intersection = new Intersection;

    intersection->amountLabels = size;
    intersection->amountTreeIntersection = 0;
    intersection->labels = new int[size];
    intersection->amountVertices = new int[size];
    intersection->positionLabelsBlock = new int*[size];
    
    for (int i  = 0; i < intersection->amountLabels; i++) {
        // Iniciar posições dos vetores com valores iniciais
        intersection->labels[i] = i;
        intersection->amountVertices[i] = 0;
        
        // Alocar espaço para linhas da matriz com apenas 1 coluna
        intersection->positionLabelsBlock[i] = new int[1];
    }

    return intersection;
}

/**
 * @brief   Alterar valor do label somando o valor máximo de label do bloco anterior
 * 
 * @param   maxLabel Mario label do bloco anterior
 * @param   intersection    Estrutura que guarda todas as árvores e posições de seus vértices
 * 
 * @details Alterar valor do label somando o valor máximo de label do bloco anterior para
 *          não gerar erros ao propagar as árvores da interseção para o próximo bloco
 * 
 * @return intersection - Estrutura que contém os labels que serão alterados
*/
Intersection* changeLabelValue(int maxLabel, Intersection* intersection) {
    for (int i = 0; i < intersection->amountLabels; i++) {
        intersection->labels[i] = intersection->labels[i] + maxLabel;
    }

    return intersection;
}

/**
 * @brief   Calcular a profundidade das árvores em cada bloco
 * 
 * @param   intersection    Estrutura que terá todos os vetores e matrizes preenchidos
 * @param   imgBlock        Bloco que terá a profundidade das árvores calculado
 * @param   num_pixels      Quantidade de pixels em cada frame do bloco
 * @param   numBlock        Qual bloco do vídeo será calculado a profundidade
 * @param   num_frames      Quantidade de frames no bloco
 * 
 * @details A profundidade da árvore é em quantas imagens ela aparece no bloco. Cada linha
 *          da matriz depthTreesBlock é referente a um bloco, e as colunas desta linha são
 *          as profundidades das árvores que estão na última imagem do bloco
*/
void getDepthTreesBlock (Intersection* intersection, Image **imgBlock, int num_pixels, int numBlock, int num_frames) {
    int i = 0, n = 0, pos = -1, size = 0;
    int pos_tree_depth = -1;
    int label_tree = -1;
    int last_img = num_frames - 1;
    int biggest_label_last_block_img = imgBlock[last_img]->val[0][0];
    int smallest_label_last_block_img = imgBlock[last_img]->val[0][0];

    // Pegar o menor e o maior label da última imagem do bloco para saber a quantidade de árvores que serão calculadas a profundidade
    for (i = 1; i < num_pixels; i++) {
        if (imgBlock[last_img]->val[i][0] > biggest_label_last_block_img) {
            biggest_label_last_block_img = imgBlock[last_img]->val[i][0];
        }
        if (imgBlock[last_img]->val[i][0] < smallest_label_last_block_img) {
            smallest_label_last_block_img = imgBlock[last_img]->val[i][0];
        }
    }   
    size = biggest_label_last_block_img - smallest_label_last_block_img + 1;
    intersection->depthTreesBlock[numBlock] = new int[size];
    // Vetor auxiliar para não calcular a profundidade de uma árvore mais de uma vez
    int* isCalculated = new int[size];
    
    // Iniciar todas as profundidades das árvores deste bloco como 0
    for (i = 0; i < size; i++) {
        intersection->depthTreesBlock[numBlock][i] = 0;
        isCalculated[i] = 0;
    }

    // Verificar se existem as árvores nas imagens do bloco e contabilizar a profundidade
    for (n = 0; n < num_pixels; n++) {
        // Pegar label do pixel atual
        //pos = imgBlock[last_img]->val[n][0] % intersection->amountLabels;
        //label_tree = intersection->labels[pos];
        pos_tree_depth = imgBlock[last_img]->val[n][0] - smallest_label_last_block_img;

        // Verificar se a profundidade da árvore já foi calculada
        if (isCalculated[pos_tree_depth] == 0) {
            isCalculated[pos_tree_depth] = 1;

            // Calcular profundidade
            for(i = 0; i < num_frames; i++) {

            }
        }
    }
    
    delete[] isCalculated;
}

/**
 * @brief   Preencher estrutura com interseção e bloco que terá as árvores propagadas
 * 
 * @param   intersectionImage   Imagem do bloco anterior que é a interseção
 * @param   intersectionBlock   Primeira imagem do próximo bloco que terá as árvores propagadas
 * @param   num_pixels          Quantidade de pixels em cada frame do bloco
 * @param   intersection        Estrutura que terá todos os vetores e matrizes preenchidos
 * 
 * @details Preencher todos os vetores e matrizes da estrutura com as imagens de interseção
 * 
 * @return intersection - Estrutura preenchida
*/
void insertIntersection (Image **intersectionImage, Image **intersectionBlock, int num_pixels, Intersection* intersection) {
    int lineIntersection = -1;
    int lineBlock = -1;
    int biggestIntersectionLabel = intersectionImage[0]->val[0][0];
    int smallestIntersectionLabel = intersectionImage[0]->val[0][0];
    int size = 0;

    // Pegar o menor e o maior label da imagem da insterseção
    for (int i = 1; i < num_pixels; i++) {
        if (intersectionImage[0]->val[i][0] > biggestIntersectionLabel) {
            biggestIntersectionLabel = intersectionImage[0]->val[i][0];
        }
        if (intersectionImage[0]->val[i][0] < smallestIntersectionLabel) {
            smallestIntersectionLabel = intersectionImage[0]->val[i][0];
        }
    }   
    // Alocar espaço para posição dos vértices das árvores na imagem de interseção
    size = biggestIntersectionLabel - smallestIntersectionLabel + 1;
    intersection->amountTreeIntersection = size;
    intersection->amountVerticesIntersection = new int[size];
    intersection->positionLabelsIntersection = new int*[size];
    intersection->biggestIntersectionLabel = biggestIntersectionLabel;
    intersection->smallestIntersectionLabel = smallestIntersectionLabel;

    // Inicializar quantidade de vértices das árvores com valor 0 e alocar espaço para linhas da matriz
    for (int i = 0; i < size; i++) {
        intersection->amountVerticesIntersection[i] = 0;
        intersection->positionLabelsIntersection[i] = new int[1];
    }
    
    // Pegar posições dos vértices das árvores da primeira imagem do bloco e da única imagem de interseção
    for (int i = 0; i < num_pixels; i++)
    {
        lineIntersection = intersectionImage[0]->val[i][0] - smallestIntersectionLabel;
        insertVerticePositionIntersection(lineIntersection, i, intersection);

        lineBlock = intersectionBlock[0]->val[i][0] % intersection->amountLabels;
        insertVerticePositionBlock (lineBlock, i, intersection);
    }
}

/**
 * @brief   Escolher label da árvore maior, caso a quantidade máxima de vértices nas posições 
 *          da árvore do próximo bloco entre duas árvores seja igual
 * 
 * @param   amountVerticesTree   Imagem do bloco anterior que é a interseção
 * @param   max   Primeira imagem do próximo bloco que terá as árvores propagadas
 * @param   intersection        Estrutura que terá todos os vetores e matrizes preenchidos
 * 
 * @details Preencher todos os vetores e matrizes da estrutura com as imagens de interseção
 * 
 * @return intersection - Estrutura preenchida
*/
int chooseLabelWithMoreVertices(Intersection* intersection, int *amountVerticesTree, int max) {
    int selectedLabel = -1;
    int maxVerticesTree = 0;

    for (int i = 0; i < intersection->amountTreeIntersection; i++) {
        if ( max == amountVerticesTree[i]) {
            if (intersection->amountVerticesIntersection[i] > maxVerticesTree) {
                maxVerticesTree = intersection->amountVerticesIntersection[i];
                selectedLabel = i + intersection->smallestIntersectionLabel;
            }
        }
    }

    return selectedLabel;
}

/**
 * @brief   Propagar árvores da interseção para o próximo bloco
 * 
 * @param   intersection        Estrutura que terá todos os vetores e matrizes preenchidos
 * 
 * @details Propagar árvores da interseção para o próximo bloco para manter a coerência temporal no vídeo
*/
void propagateIntersectingTrees (Intersection* intersection, Image **intersectionImage, int propagatedMoreLabels) {
    int *amountVerticesTree = new int[intersection->amountTreeIntersection];
    int nVerticesBlockTree = 0;
    int max = 0;
    int labelTree = -1;
    int labelPropagated = -1;
    int amountLabelPropagated = 0;
    int *labelsPropagated = new int[intersection->biggestIntersectionLabel + 1];
    int amountTreeFirstImgBlock = 0;
    int **choiceLabels = new int*[intersection->biggestIntersectionLabel + 1];
    int *amountChoiceLabels = new int[intersection->biggestIntersectionLabel + 1];
    int maxChoiceLabelBlock = 0;
    int choiceLabelBlock = -1;


    for (int i = 0; i <= intersection->biggestIntersectionLabel; i++) {
        labelsPropagated[i] = -1;
        // Estrutura que guarda as árvores que irão receber o mesmo label, para ser escolhida qual delas irá receber o label da interseção
        choiceLabels[i] = new int[intersection->amountLabels * 2];
        amountChoiceLabels[i] = 0;
    }

    for (int i = 0; i < intersection->amountLabels; i++) {
        max = 0;
        labelPropagated = -1;
        // Iniciar quantidade de vértices das árvores da interseção com 0
        for (int x = 0; x < intersection->amountTreeIntersection; x++) {
            amountVerticesTree[x] = 0;
        }

        // Contabilizar quantos vértices de cada árvore estão na mesma posição dos vértices da primeira
        // imagem do próximo bloco
        nVerticesBlockTree = intersection->amountVertices[i];
        for (int j = 0; j < nVerticesBlockTree; j++) {
            // Ver qual label de árvore da interseção está na posição do label da árvore do bloco da linha i
            labelTree = intersectionImage[0]->val[intersection->positionLabelsBlock[i][j]][0];
            // Adicionar +1 na posição referente a este label da interseção
            amountVerticesTree[labelTree - intersection->smallestIntersectionLabel]++;
            // Encontrar a maior quantidade de vértices da árvore da interseção
            if (amountVerticesTree[labelTree - intersection->smallestIntersectionLabel] > max) {
                max = amountVerticesTree[labelTree - intersection->smallestIntersectionLabel];
            }
        }
        if (nVerticesBlockTree != 0 && max != 0) {
            amountTreeFirstImgBlock++;
            // Caso tenha a mesma quantidade de vértices entre labels da interseção, escolher a maior árvore
            labelPropagated = chooseLabelWithMoreVertices(intersection, amountVerticesTree, max);
            labelsPropagated[labelPropagated] = 1;
            if (propagatedMoreLabels != 1) {
                choiceLabels[labelPropagated][amountChoiceLabels[labelPropagated]] = i;
                choiceLabels[labelPropagated][amountChoiceLabels[labelPropagated] + 1] = max;
                amountChoiceLabels[labelPropagated] += 2;
            } else {
                // Alterar label, propagando o label da árvore escolhida da interseção
                intersection->labels[i] = labelPropagated;
            }         
        }
    }

    if (propagatedMoreLabels != 1) {
        // Percorrer as linhas de choiceLabels que cada posição representa um label da interseção que foi selecionado para propagar
        for (int i = 0; i <= intersection->biggestIntersectionLabel; i++) {
            // Se o label da interseção foi escolhido para propagar, ou seja, amountChoiceLabels na posição i recebe um label que receberá a propagação
            if (amountChoiceLabels[i] > 0) {
                // Escolher até então o primeiro label do block como o que possui o maior número de pixels na árvore da interseção
                choiceLabelBlock = choiceLabels[i][0];
                maxChoiceLabelBlock = choiceLabels[i][1];
                
                // Selecionar dentre os labels do bloco qual possui o maior número de pixels na mesma posição do label da árvore da interseção
                for (int j = 2; j < amountChoiceLabels[i]; j += 2) {
                    if (choiceLabels[i][j + 1] > maxChoiceLabelBlock) {
                        choiceLabelBlock = choiceLabels[i][j];
                        maxChoiceLabelBlock = choiceLabels[i][j + 1];
                    }
                }
                // Alterar label, propagando o label da árvore escolhida da interseção
                intersection->labels[choiceLabelBlock] = i;
            }
        }
    }

    for (int i = 0; i <= intersection->biggestIntersectionLabel; i++) {
        if (labelsPropagated[i] != -1) {
            amountLabelPropagated++;
        }
    }
    //printf("\nQuantidade de labels propagados: %d\n", amountLabelPropagated);
    //printf("Quantidade de árvores na primeira imagem do bloco que recebeu os labels: %d\n", amountTreeFirstImgBlock);

    delete[] labelsPropagated;
    delete[] amountVerticesTree;
}

/**
 * @brief   Alocar espaço para uma posição de vértice da árvore do bloco que terá a propagação
 * 
 * @param   line            Qual linha está a árvore que terá a posição do vértice adicionada
 * @param   position        A posição do vértice na árvore
 * @param   intersection    Estrutura que guarda todas as árvores e posições de seus vértices
 * 
 * @details Alocar espaço para uma coluna da linha do label que representa uma árvore
 * 
 * @return insertedColumn - Coluna em que a posição do vértice foi adicionada
*/
int insertVerticePositionBlock (int line, int position, Intersection* intersection) {
    int insertedColumn = -1;
    // Alocar espaço para colunas anteriores + 1 coluna
    int* newLine = new int[intersection->amountVertices[line] + 1];

    if (newLine != nullptr) {
        // Copia os elementos da linha desatualizada para a nova linha
        for (int j = 0; j < intersection->amountVertices[line]; j++) {
            newLine[j] = intersection->positionLabelsBlock[line][j];
        }

        // Insere a posição na nova coluna
        newLine[intersection->amountVertices[line]] = position; 
        // Libera a memória da linha desatualizada
        delete[] intersection->positionLabelsBlock[line]; 
        // Atualiza a linha na matriz
        intersection->positionLabelsBlock[line] = newLine;
        insertedColumn = intersection->amountVertices[line];
        // Atualiza o número de colunas
        intersection->amountVertices[line] += 1;
    }

    return insertedColumn;
}

/**
 * @brief   Alocar espaço para uma posição de vértice da árvore da interseção
 * 
 * @param   line            Qual linha está a árvore que terá a posição do vértice adicionada
 * @param   position        A posição do vértice na árvore
 * @param   intersection    Estrutura que guarda todas as árvores e posições de seus vértices
 * 
 * @details Alocar espaço para uma coluna da linha do label que representa uma árvore
 * 
 * @return insertedColumn - Coluna em que a posição do vértice foi adicionada
*/
int insertVerticePositionIntersection (int line, int position, Intersection* intersection) {
    int insertedColumn = -1;

    // Alocar espaço para colunas anteriores + 1 coluna
    int* newLine = new int[intersection->amountVerticesIntersection[line] + 1];

    if (newLine != nullptr) {
        // Copia os elementos da linha desatualizada para a nova linha
        for (int j = 0; j < intersection->amountVerticesIntersection[line]; j++) {
            newLine[j] = intersection->positionLabelsIntersection[line][j];
        }

        // Insere a posição na nova coluna
        newLine[intersection->amountVerticesIntersection[line]] = position; 
        // Libera a memória da linha desatualizada
        delete[] intersection->positionLabelsIntersection[line]; 
        // Atualiza a linha na matriz
        intersection->positionLabelsIntersection[line] = newLine;
        insertedColumn = intersection->amountVerticesIntersection[line];
        // Atualiza o número de colunas
        intersection->amountVerticesIntersection[line] += 1;
    }

    return insertedColumn;
}

/**
 * @brief   Liberar memória das estruturas que são utilizadas para 
 *          propagar árvores de interseção para o próximo bloco
 * 
 * @param   intersection Estrutura que será liberada
 * 
 * @details Liberar memória das estruturas que são utilizadas para 
 *          propagar árvores de interseção para o próximo bloco
*/
void deleteIntersection(Intersection* intersection) {
    int nVertices = 0;
    int nVerticesIntersection = 0;

    for (int j = 0; j < intersection->amountLabels; j++) {

        delete[] intersection->positionLabelsBlock[j];
        intersection->amountVertices[j] = 0;
    }

    if (intersection->amountTreeIntersection > 0) {
        for (int j = 0; j < intersection->amountTreeIntersection; j++) {
            delete[] intersection->positionLabelsIntersection[j];
            intersection->amountVerticesIntersection[j] = 0;
        }
        delete[] intersection->amountVerticesIntersection;
        delete[] intersection->positionLabelsIntersection;
    }

    delete[] intersection->positionLabelsBlock;
    delete[] intersection->amountVertices;
    delete[] intersection->labels;
    intersection->amountLabels = 0;
    intersection->amountTreeIntersection = 0;
    delete intersection;
}

/**
 * @brief   Printar os valores dos labels atualizados
 * 
 * @param   intersection Estrutura que guarda todas as árvores e posições de seus vértices
 * 
 * @details Printar os valores dos labels atualizados
*/
void printLabels (Intersection* intersection) {

    for (int i = 0; i < intersection->amountLabels; i++) {
        printf("%d -> %d\t\t", i, intersection->labels[i]);
    }
}

/**
 * @brief   Retorna o label de maior valor do bloco
 * 
 * @param   intersection Estrutura que guarda todas as árvores e posições de seus vértices
 * 
 * @details Retorna o label de maior valor do bloco
 * 
 * @return max - Retorna o label de maior valor
*/
int getMaxNumLabel(Intersection* intersection) {
    int max = 0;

    for (int i = 0; i < intersection->amountLabels; i++)
    {
        if(intersection->labels[i] > max) {
            max = intersection->labels[i];
        }
    }
    //printf("maior label: %d\n", max);
    return max;
}

/**
 * @brief   Printar a posição dos vértices nas árvores da interseção
 * 
 * @param   intersection Estrutura que guarda todas as árvores e posições de seus vértices
 * 
 * @details Printar a posição dos vértices nas árvores da interseção
*/
void printPositionVerticesIntersection (Intersection* intersection) {
    int columns = 0;

    printf("\nPosicao dos vertices na intersecao: ");
    for (int i = 0; i < intersection->amountTreeIntersection; i++)
    {   
        printf("\nLabel: %d --> ", (intersection->smallestIntersectionLabel + i));
        columns = intersection->amountVerticesIntersection[i];
        for (int j = 0; j < columns; j++) {
            printf("%d - ", intersection->positionLabelsIntersection[i][j]);
        }
    }
}

/**
 * @brief   Printar a posição dos vértices nas árvores da primeira imagem do bloco
 * 
 * @param   intersection Estrutura que guarda todas as árvores e posições de seus vértices
 * 
 * @details Printar a posição dos vértices nas árvores da primeira imagem do bloco
*/
void printPositionVerticesBlock (Intersection* intersection) {
    int columns = 0;

    printf("\nPosicao dos vertices na primeira imagem do bloco: ");
    for (int i = 0; i < intersection->amountLabels; i++)
    {
        printf("\nLabel: %d --> ", intersection->labels[i]);
        columns = intersection->amountVertices[i];
        for (int j = 0; j < columns; j++) {
            printf("%d - ", intersection->positionLabelsBlock[i][j]);
        }
    }
}

/*
void propagateIntersectingTrees(Image **original, Image **intersecao, int num_rows, int num_cols, int num_frames, int nImagesIntersecao) {
    int num_pixels = num_cols * num_rows;
    // labelDiff contem o label que estava e o proximo no vetor é o que será colocado, então onde tiver o label que estava será alterado para o propagado
    int *labelDiff = (int *)calloc(num_pixels * 2, sizeof(int));
    int pos = 0;

    for (int f = 0; f < nImagesIntersecao - 1; f++)
    {
        for (int i = 0; i < num_pixels; i++)
        {
            original[f]->val[i][0] = intersecao[f]->val[i][0];
        }
    }

    // Armazenar quais labels precisam ser substituídos caso árvores fora da interseção tenham nós dentro do último frame da interseção
    for (int i = 0; i < num_pixels; i++)
    {
        if (original[nImagesIntersecao - 1]->val[i][0] == original[nImagesIntersecao]->val[i][0]) {
            // Valor do label anterior
            labelDiff[pos] = original[nImagesIntersecao - 1]->val[i][0];
            // Valor que veio do label da intersecao
            labelDiff[pos+1] = intersecao[nImagesIntersecao - 1]->val[i][0];
            pos = pos + 2;
        }
        // Propagar os labels da interseção do último frame da interseção para o frame do bloco
        original[nImagesIntersecao - 1]->val[i][0] = intersecao[nImagesIntersecao - 1]->val[i][0];
    }

    // Substituir quais labels de árvores precisam ser atualizados para propagar para árvores fora da interseção
    for (int k = 0; k < pos; k = k + 2) {
        for (int f = nImagesIntersecao; f < num_frames; f++) {
            for (int i = 0; i < num_pixels; i++) {

                if (original[f]->val[i][0] == labelDiff[k]) {
                    original[f]->val[i][0] = labelDiff[k+1];
                }
            }
        }
    }
}*/

/*
int getMaxNumLabel(Image **label_video, int num_pixels, int frames) {
    int max = 0;

    for (int f = 0; f < frames; f++)
    {
        for (int i = 0; i < num_pixels; i++)
        {
            if(label_video[f]->val[i][0] > max) {
                max = label_video[f]->val[i][0];
            }
        }
    }
    printf("maior label: %d\n", max);
    return max;
}
*/
/*
void changeLabelValue(Image **label_video, int num_pixels, int frames, int maxLabel) {

    for (int f = 0; f < frames; f++)
    {
        for (int i = 0; i < num_pixels; i++)
        {
            label_video[f]->val[i][0] = label_video[f]->val[i][0] + maxLabel;   
        }
    }
}*/