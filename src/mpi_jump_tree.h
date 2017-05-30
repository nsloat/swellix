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
void get_work(int amount, config* seq, global* crik, local* todd);
void get_mpi_globals();
int is_work_needed();
void make_jump_tree_parallel(config* seq, global* crik, local* todd, knob** cmpnts, int* cmpntTracker);
void send_work(global* crik, local* todd, int to);
void unpack_todd(global* crik, local* todd, int** msg);
void unpack_structure(global* crik, int* msg, int count);
void unpack_crikInfo(global* crik, int** msg);
void unpack_swellix_structure(global* crik, int* msg, int count);

#define MPI_REQUEST_WORK_1 (1)
#define MPI_REQUEST_WORK_2 (2)
#define MPI_WORK (3)
#define MPI_DIE (4)
//#define MPI_WORK_LEFT (1)

#define KNOB_SIZE_INT (18)

#define CHUNK_SIZE 64

#define mpi_stage_1 5
#define mpi_stage_2 6
#define mpi_change_stage 7

#endif // __MPI_JUMP_TREE
#endif
