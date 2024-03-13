#include "Utils.h"

//=============================================================================
// Int*
//=============================================================================
int* randomSampling(const int lb, const int hb, const int nsampl)
{
    int num_qtd, nsel;
    bool *was_selected;
    int *indexes, *selected;

    srand(time(NULL));

    num_qtd = hb - lb + 1;
    indexes = calloc(num_qtd, sizeof(int));
    was_selected = calloc(num_qtd, sizeof(bool));
    selected = calloc(nsampl, sizeof(int));

    for(int i = 0; i < num_qtd; ++i)
    {
        indexes[i] = lb + i;
        was_selected[i] = false;
    }

    nsel = 0;
    while(nsel < nsampl)
    {
        int index = rand() % num_qtd;
        if(!was_selected[index])
        {
            selected[nsel] = indexes[index];
            was_selected[index] = true;
            nsel++;
        }
    }

    free(indexes);
    free(was_selected);

    return selected;
}

//=============================================================================
// Void
//=============================================================================
void printError(const char* function_name, const char* message, ...)
{
    va_list args;
    char full_msg[4096];

    va_start(args, message);
    vsprintf(full_msg, message, args);
    va_end(args);

    fprintf(stderr, "\nError in %s:\n%s!\n", function_name, full_msg);
    fflush(stdout);
    exit(-1);
}

void printWarning(const char *function_name, const char *message, ...)
{
    va_list args;
    char full_msg[4096];

    va_start(args, message);
    vsprintf(full_msg, message, args);
    va_end(args);

    fprintf(stdout, "\nWarning in %s:\n%s!\n", function_name, full_msg);
}