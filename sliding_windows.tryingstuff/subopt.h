#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>


#if defined(_MPI)

#include "mpi.h"
#define MPI_RANK_MASTER (0)
#define MPI_FIVE (6)
int mpi_rank, mpi_size;

#endif // _MPI

typedef struct {
    int noLP;
    int noGU;
    int max_dist;
    int count;
    int constraints;
    int noML;
    int allpair;
    FILE * infile;
    FILE * outfile;
    int helix_badnesses;
    int helix_number;
    int helix_size;
    long statecount;
} options;

//options OPTIONS;

#define print(...) fprintf(OPTIONS.outfile, __VA_ARGS__)

typedef struct interval_struct interval;

int MODIFIED_BASE = -3;
int UNPAIRED_BASE = -2;

typedef enum {
  TIP, MULTI, LOOP,
} interval_flag;

typedef enum {
    IDLE, HALT, STATE
} mpi_flag;

typedef int bool;

#define TRUE	  1
#define FALSE	  0
#define VERYBIG   1000000000

struct interval_struct {
    int i;
    int j;
    interval_flag flag;
    interval *next;
};


typedef struct state_struct {
    interval *intervals;
    int *structure;
    int length; 
    char *sequence;
    int* constraints;
} state;

#define MIN_PAIR_DIST 3

void start(char *seq, int *constraints);
void start_with_constraints(state *s, int i);
void refine_state(state *s);
void refine_state_locally(state *s);

interval *make_interval(state *s, int ai, int aj, interval_flag aflag)  ;
void unmake_interval(state *s);

void print_soln(state *s);
void print_state(state *s);

int allow_pair(state *s, int i, int j); 
int can_pair(char a, char b);
int constrained_to_not_pair(state *s, int i, int j);
int chemical_modification_ok(state *s, int i, int j);

int allow_state(state *s); 
int room_for_helices(state *s);
int helix_count(state *s);
int lonely_pair(state *s);

void set_args(int argc, char *argv[]);
void print_usage();
int *interpreted_constraints(char *constraints);

void pstart(char *seq, int *constraints);
int * pack_state(state * S);
int *wait_for_state();
state *interpret_message(int *ints);
/* int* package_state(state *s); // MALLOC */
/* state *unpack_state(int *s); */
/* int send_state(state *s, int destination); */
/* state *receive_state(int from); */
/* int check_for_idle(int *idleness); */
/* int next_idle_node(); */


int label_struct(char* structure, int start, char* outprefix);
