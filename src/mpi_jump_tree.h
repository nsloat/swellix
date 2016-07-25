#ifndef __MPI_JUMP_TREE
#define __MPI_JUMP_TREE

//*
//* This is nearly an exact algorithm transportation from Parallel Crumple. The main/only differences
//* will be in data structures transferred between processing elements and the requisite code to facilitate
//* these transferral.
//*

int commit_datatypes();
void get_work(/*params*/, int amount);
void get_mpi_globals();
int is_work_needed();
void make_jump_tree_parallel(/*params*/);
void send_work(/*params*/, int to);

#include "mpi.h"
#include "main.h"
#include "jump_tree.h"

#define MPI_REQUEST_WORK (1)
#define MPI_WORK (2)
#define MPI_DIE (4)
#define MPI_WORK_LEFT (1)


#endif // __MPI_JUMP_TREE
