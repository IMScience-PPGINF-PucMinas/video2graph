#include "PrioQueue.h"

//=============================================================================
// Bool functions
//=============================================================================
bool insertPrioQueue(PrioQueue **queue, int index)
{
    bool success;

    if(isPrioQueueFull(*queue)) 
    {
        printWarning("insertPrioQueue", "The queue is full");
        success = false;
    }
    else
    {
        PrioQueue *tmp;

        tmp = *queue;

        // Number of inserted elements + 1
        tmp->last_elem_pos++;
        
        // Put the new index in the last position
        tmp->node[tmp->last_elem_pos] = index;
        tmp->state[index] = GRAY_STATE; // Newly inserted
        tmp->pos[index] = tmp->last_elem_pos;

        // Update its position to its correct location
        moveIndexUpPrioQueue(queue, index);

        success = true;
    }

    return success;
}

inline bool isPrioQueueEmpty(PrioQueue *queue)
{
    return (queue->last_elem_pos == -1);
}

inline bool isPrioQueueFull(PrioQueue *queue)
{
    return (queue->last_elem_pos == queue->size - 1);
}

//=============================================================================
// Int functions
//=============================================================================
inline int getFatherPos(int pos)
{
    return (pos - 1)/2; 
}

inline int getLeftSonPos(int pos)
{
    return (2 * pos + 1); 
}

inline int getRightSonPos(int pos)
{
    return (2 * pos + 2); 
}

int popPrioQueue(PrioQueue **queue)
{    
    int index;
    
    if(isPrioQueueEmpty(*queue)) 
    {
        printWarning("popPrioQueue", "The queue is empty");
        index = -1;
    }
    else
    {
        PrioQueue *tmp;

        tmp = *queue;

        // Auxiliary var
        index = tmp->node[0];

        // Frees the current position
        tmp->pos[index] = -1;
        tmp->state[index] = BLACK_STATE; // Ordely removed

        // Puts the last index in the first position
        tmp->node[0] = tmp->node[tmp->last_elem_pos];
        tmp->pos[tmp->node[0]] = 0;

        // Frees the last position
        tmp->node[tmp->last_elem_pos] = -1;
        
        // Number of inserted elements - 1
        tmp->last_elem_pos--;

        // Updates the moved node to its correct position
        moveIndexDownPrioQueue(queue, tmp->node[0]);
    }

    return index;
}

//=============================================================================
// PrioQueue* functions
//=============================================================================
PrioQueue* createPrioQueue(int size, double *prio, RemPolicy rem_policy)
{
    PrioQueue *queue;

    queue = (PrioQueue*)calloc(1,sizeof(PrioQueue));

    queue->size = size;
    queue->prio = prio;
    queue->state = (ElemState*)calloc(size, sizeof(ElemState));
    queue->node = (int*)calloc(size, sizeof(int));
    queue->pos = (int*)calloc(size, sizeof(int));
    queue->last_elem_pos = -1; // It is empty
    queue->rem_policy = rem_policy;

    // Set heap to its initial values
    for(int i = 0; i < queue->size; i++)
    {
        queue->state[i] = WHITE_STATE;
        queue->pos[i] = -1;
        queue->node[i] = -1;
    }

    return queue;
}

//=============================================================================
// Void functions
//=============================================================================
void freePrioQueue(PrioQueue **queue)
{
    if(*queue != NULL)
    {
        PrioQueue *tmp;

        tmp = *queue;

        free(tmp->state);
        free(tmp->node);
        free(tmp->pos);
        free(*queue);
    }
}

void moveIndexDownPrioQueue(PrioQueue **queue, int index)
{    
    int pos, curr_pos, left_pos, right_pos, curr_index;
    float curr_prio;
    PrioQueue *tmp;

    tmp = *queue;

    if(index < tmp->size && index >= 0) 
    {
        // Auxiliary vars of the current index
        pos = tmp->pos[index];

        // Positions of the respective sons' nodes
        left_pos = getLeftSonPos(pos);
        right_pos = getRightSonPos(pos);

        // Evaluation vars
        curr_pos = pos;
        curr_index = index;
        curr_prio = tmp->prio[curr_index];

        // Check the left son
        if(left_pos <= tmp->last_elem_pos)
        {
            int left_index;
            float left_prio;

            left_index = tmp->node[left_pos];
            left_prio = tmp->prio[left_index];

            if((tmp->rem_policy == MINVAL_POLICY && left_prio < curr_prio) ||
               (tmp->rem_policy == MAXVAL_POLICY && left_prio > curr_prio))
            {
                curr_pos = left_pos;
                curr_index = left_index;
                curr_prio = left_prio;
            }
            else if(left_prio == curr_prio)
            {
                curr_pos = left_pos;
                curr_index = left_index;
                curr_prio = left_prio;
            }
        }

        // Check the right son
        if(right_pos <= tmp->last_elem_pos)
        {
            int right_index;
            float right_prio;

            right_index = tmp->node[right_pos];
            right_prio = tmp->prio[right_index];

            if((tmp->rem_policy == MINVAL_POLICY && right_prio < curr_prio) ||
               (tmp->rem_policy == MAXVAL_POLICY && right_prio > curr_prio))
            {
                curr_pos = right_pos;
                curr_index = right_index;
                curr_prio = right_prio;
            }
            else if(right_prio == curr_prio)
            {
                curr_pos = right_pos;
                curr_index = right_index;
                curr_prio = right_prio;
            }
        }

        // Swap
        if(curr_pos != pos)
        {
            tmp->node[curr_pos] = index;
            tmp->pos[curr_index] = pos;

            tmp->node[pos] = curr_index;
            tmp->pos[index] = curr_pos;

            // Order through recursion
            moveIndexDownPrioQueue(queue, index);
        }
    }
}

void moveIndexUpPrioQueue(PrioQueue **queue, int index)
{
    int pos, curr_pos, curr_index, father_pos;
    float prio, curr_prio;
    PrioQueue *tmp;

    tmp = *queue;

    if(index < tmp->size && index >= 0) 
    {
        // Auxiliary vars of the current index
        pos = tmp->pos[index];
        prio = tmp->prio[index];

        // Position of the father's node
        father_pos = getFatherPos(pos);

        // Evaluation vars
        curr_pos = pos;
        curr_index = index;
        curr_prio = prio;

        // Check the father node
        if(father_pos >= 0)
        {
            int father_index;
            float father_prio;

            father_index = tmp->node[father_pos];
            father_prio = tmp->prio[father_index];

            if((tmp->rem_policy == MINVAL_POLICY && father_prio > curr_prio) ||
               (tmp->rem_policy == MAXVAL_POLICY && father_prio < curr_prio))
            {
                curr_pos = father_pos;
                curr_index = father_index;
                curr_prio = father_prio;
            }
        }

        // Swap
        if(curr_pos != pos)
        {
            tmp->node[curr_pos] = index;
            tmp->pos[curr_index] = pos;

            tmp->node[pos] = curr_index;
            tmp->pos[index] = curr_pos;

            // Order through recursion
            moveIndexUpPrioQueue(queue, index);
        }
    }
}
void removePrioQueueElem(PrioQueue **queue, int index)
{    
    double prio;
    PrioQueue *tmp;

    tmp = *queue;

    prio = tmp->prio[index];

    // Guaranteeing that the desired index will be returned in POP
    if(tmp->rem_policy == MINVAL_POLICY)
        tmp->prio[index] = -INFINITY;
    else
        tmp->prio[index] = INFINITY;

    // Get the desired index
    moveIndexUpPrioQueue(queue, index);
    popPrioQueue(queue); // No desire for returning the index

    // Change back the priority value
    tmp->prio[index] = prio;
    // It was not orderly removed
    tmp->state[index] = WHITE_STATE;
}

void resetPrioQueue(PrioQueue **queue)
{
    PrioQueue *tmp;

    tmp = *queue;

    // Clears all the structures of the heap
    for(int i = 0; i < tmp->size; i++)
    {
        tmp->state[i] = WHITE_STATE;
        tmp->pos[i] = -1;
        tmp->node[i] = -1;
    }
    // No element inserted
    tmp->last_elem_pos = -1;
}