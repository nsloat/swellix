#ifdef _MPI
#ifndef __MPI_JUMP_TREE
#define __MPI_JUMP_TREE

#include "main.h"

//*
//* This is nearly an exact algorithm transportation from Parallel Crumple. The main/only differences
//* will be in data structures transferred between processing elements and the requisite code to facilitate
//* these transferral.
//*

int commit_datatypes();
void get_work(/*params*/, int amount);
void get_mpi_globals();
int is_work_needed();
void make_jump_tree_parallel(config* seq, global* crik, knob** cmpnts, int* cmpntTracker);
void send_work(/*params*/, int to);
int* pack_state(/*params*/);
void unpack_state(/*params*/);

#define MPI_REQUEST_WORK (1)
#define MPI_WORK (2)
#define MPI_DIE (4)
//#define MPI_WORK_LEFT (1)

#define KNOB_SIZE_INT (19)

#endif // __MPI_JUMP_TREE
#endif
