#ifndef PARALLEL
#define PARALLEL

void get_work(state_stack *stack, int amount);
void get_mpi_globals();
int is_work_needed();
void fold_all_parallel(state_stack *stack);
int *pack_state(state * s, int length);
state *unpack_state(int *ints);
void send_work(state_stack *stack, int to);

//------------------- MPI jazz for parallelization  ------------------//
// So, there may be some mystery about why this communication scheme is
// as complicated as it is. Mostly, there are some terrifying race conditions.
#include "mpi.h"

#define MPI_REQUEST_WORK (1)
#define MPI_WORK (2)
#define MPI_DIE (4)

int mpi_rank, mpi_from, mpi_to, mpi_size;


void get_mpi_globals()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    mpi_from = mpi_rank - 1;
    if (mpi_from < 0) {
        mpi_from = mpi_size - 1;
    }
    mpi_to = mpi_rank + 1;
    if (mpi_to == mpi_size) {
        mpi_to = 0;
    }
}

void free_state(state *s){
    free(s->intervals);
    free(s->structure);
    free(s);
}

void fold_parallel(state_stack *stack)
{
    if (mpi_rank == 0) {
        new_state(stack);
        fold_all_parallel(stack);
        get_work(stack, 0);
    } else {
        int message[4] = {mpi_rank, 0, 0, 0};
        MPI_Ssend(message, 4, MPI_INT, mpi_to, MPI_REQUEST_WORK,
                  MPI_COMM_WORLD);
        get_work(stack, -1);
    }

    while (stack->current >= 0) {
        fold_all_parallel(stack);
        get_work(stack, 0);
    }
}

void fold_all_parallel(state_stack *stack)
{
    long counter = 0;
    while(stack->current >= 0){        
        if(counter % 64 == 63){
            int to = is_work_needed(stack);
            if(to != -1){
                send_work(stack, to);
            }
        }
        counter++;
        fold1(stack);
    }
}

// message format is: 
// { message source, states sent left, states sent, states recieved }

int parallel_sent_work_left = 0;
int parallel_work_received = 0;
int parallel_work_sent = 0;

int is_work_needed(state_stack *stack)
{
    MPI_Status ms;
    int flag = 0;
    MPI_Iprobe(mpi_from, MPI_REQUEST_WORK, MPI_COMM_WORLD, &flag, &ms);

    if (flag) { 
        int message[4];
        MPI_Recv(message, 4, MPI_INT, ms.MPI_SOURCE,
                 MPI_REQUEST_WORK, MPI_COMM_WORLD, &ms);

        //should I send work?
        if (stack->states[0].unpair == 1 ||
            stack->current > 0)
            return message[0];

        // or should I pass the message along?
        if (message[0] == 0){
            message[2] += 10000;
        }
        // This send could fail for process 0 on the first send. Hmm. 
        MPI_Send(message, 4, MPI_INT, mpi_to, MPI_REQUEST_WORK, MPI_COMM_WORLD);
    }
    return -1;
}


void prepare_to_work(MPI_Status ms, state_stack *stack)
{
    parallel_work_received++;

    int num_elements;
    MPI_Get_count(&ms, MPI_INT, &num_elements);
    int *message = malloc(sizeof(int) * (num_elements+1));
    MPI_Recv(message, num_elements, MPI_INT, ms.MPI_SOURCE,
             MPI_WORK, MPI_COMM_WORLD, &ms);

    state *b = unpack_state(message);
    duplicate_state(b, stack); // push it onto the stack.

    free(message);
    free_state(b);
}


void prepare_to_die()
{
    int message = 0;
    MPI_Send(&message, 1, MPI_INT, mpi_to, MPI_DIE, MPI_COMM_WORLD);
}


void get_work(state_stack *stack, int amount)
{
    //print("I RAN OUT OF STATES!\n");
    int message[4] = { mpi_rank, 0, 0, 0 };

    if (amount != -1) {
        if (mpi_rank == 0) {
            message[1] = 0;     // 0 can't send to the left.
            message[2] = parallel_work_sent;
            message[3] = parallel_work_received;
            MPI_Send(message, 4, MPI_INT, mpi_to, MPI_REQUEST_WORK,
                     MPI_COMM_WORLD);
        } else {
            MPI_Send(message, 4, MPI_INT, mpi_to, MPI_REQUEST_WORK,
                     MPI_COMM_WORLD);
        }
    }

    MPI_Status ms;
    while (1) {
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ms);
        //print("%d is listening.\n", mpi_rank);
        // possible tags are: WORK, REQUEST_WORK, and DIE.
        if (ms.MPI_TAG == MPI_WORK) {
            prepare_to_work(ms, stack);
            return;

        } else if (ms.MPI_TAG == MPI_REQUEST_WORK) {
            MPI_Recv(message, 4, MPI_INT, mpi_from, MPI_REQUEST_WORK,
                     MPI_COMM_WORLD, &ms);

            if (message[0] == 0 && mpi_rank == 0) {
                // If zero gets back it's own request, and
                // the flags say no work was en route, 
                // the DIE.
                //print("zero got back the flag.\n");
                if (message[1] == 0 && message[2] == message[3]) {
                    prepare_to_die();
                    return;
                }
                message[1] = 0;
                message[2] = 0;
                message[3] = 0;
            }

            if (message[0] == 0) {
                // If it came from 0, he wants extra info.
                //print("it's a flag.\n");
                message[1] += parallel_sent_work_left;
                message[2] += parallel_work_sent;
                message[3] += parallel_work_received;
                parallel_sent_work_left = 0;
                parallel_work_sent = 0;
                parallel_work_received = 0;
            }

            MPI_Send(message, 4, MPI_INT, mpi_to, MPI_REQUEST_WORK,
                     MPI_COMM_WORLD);

        } else if (ms.MPI_TAG == MPI_DIE) {
            prepare_to_die();
            return;
        }
    }
}


//-------------------------- SPLITTING AND SENDING ---------------------------//

state * split_stack(state_stack *stack)
{
    // This is the state closest to the root of the tree.
    // It will be the most fruitful to split.
    state *s = &stack->states[0];


    // Copy the topmost state.
    // This is a bit of a mess because we're copying OUTSIDE the stack.
    state *news = malloc(sizeof(state));

    news->intervals = malloc(sizeof(interval)*stack->length);
    news->structure = malloc(sizeof(int)*stack->length);
    int i;
    for (i=0; i<stack->length; i++){
        news->structure[i] = s->structure[i];
    }
    for (i=0; i<=s->current_interval; i++){
        news->intervals[i].start = s->intervals[i].start;
        news->intervals[i].end = s->intervals[i].end;
    }
    news->current_interval = s->current_interval;
    news->current_base = s->current_base;
    news->max_base = s->max_base;
    news->unpair = s->unpair;

    // OK, now there are three different things we may have to split.
    // First, a state that is currently pairing and will eventually not-pair:
    if(s->unpair == 1){
        s->unpair = 0;
        news->current_base = news->intervals[news->current_interval].end; 
    // Next, a state that's pairing but will not not-pair:
    } else if (s->unpair == 0 && s->max_base > s->current_base){
        int curmax = s->max_base;
        s->max_base = curmax-1;
        news->current_base = curmax-1;
    } else {
        shift_stack_up(stack);
        free_state(news);
        return split_stack(stack);
    }
    return news;
}


void send_work(state_stack *stack, int to)
{
    //print(" I am sending\n");
    parallel_work_sent++;

    if (to < mpi_rank) {
        parallel_sent_work_left = 1;
    }

    state *s = split_stack(stack);
    int *p = pack_state(s, stack->length);
    // magic numbers. Fix.
    MPI_Send(p, (5+stack->length+stack->length*2), MPI_INT, to, MPI_WORK,
             MPI_COMM_WORLD);
    free(p);
    free_state(s);
    //print_stack(stack);
}

//---------------------------- SERIALIZATION OF STATES ------------------------// 

int *pack_state(state * s, int length)
{
    int *m = malloc(sizeof(int)*(6+length+length*2));
    // { length, current_interval, current_base, unpair, max_base, flag,
    // structure[length], [i, j]*length 
    m[0] = length;
    m[1] = s->current_interval;
    m[2] = s->current_base;
    m[3] = s->unpair;
    m[4] = s->max_base;

    int i;
    for(i=0; i<length; i++){
        m[i +6] = s->structure[i];
    }

    for(i=0; i<length; i++){
        m[2*i   +length+6] = s->intervals[i].start;
        m[2*i+1 +length+6] = s->intervals[i].end;
    }

    return m;
}


state *unpack_state(int *m)
{
    state *s = malloc(sizeof(state));
    int length           = m[0];
    s->current_interval  = m[1];
    s->current_base      = m[2];
    s->unpair            = m[3];
    s->max_base          = m[4];

    s->structure = malloc(sizeof(int)*length);
    int i;
    for(i=0; i<length; i++){
        s->structure[i] = m[i +6];
    }

    s->intervals = malloc(sizeof(interval)*length);
    for(i=0; i<length; i++){
        s->intervals[i].start = m[2*i   +length+6];
        s->intervals[i].end   = m[2*i+1 +length+6];
    }

    return s;
}

 #endif 
