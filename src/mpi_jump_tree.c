#include "mpi_jump_tree.h"

extern int rank, wsize;
int mpi_from, mpi_to;

void get_mpi_globals() {
  mpi_from = rank-1;
  if(mpi_from < 0) {
    mpi_from = wsize-1;
  }
  mpi_to = rank+1;
  if(mpi_to == wsize) {
    mpi_to = 0;
  }
}

void make_jump_tree_parallel(config* seq, global* crik) {
}
