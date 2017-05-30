#CFLAGS=-Wall
#DEBUG_FLAGS=-g3
#PROFILE_FLAGS= -fprofile-arcs -ftest-coverage -pg
MPI_FLAGS=-D_MPI

all:
	gcc $(CFLAGS) -o crumple crumple.c 

debug:
	gcc $(CFLAGS) -o crumple crumple.c $(DEBUG_FLAGS)

profile:
	gcc $(CFLAGS) -o crumple crumple.c $(PROFILE_FLAGS) 

cc:
	gcc $(CFLAGS) -o constraintcheck constraintcheck.c

ccdbg:
	gcc $(CFLAGS) $(DEBUG_FLAGS) -o constraintcheck constraintcheck.c

ccprofile:
	gcc $(CFLAGS) $(DEBUG_FLAGS) $(PROFILE_FLAGS) -o constraintcheck constraintcheck.c

mpi:
	cc $(CFLAGS) -o crumple crumple.c $(MPI_FLAGS)

mpi-debug:
	export MPI_VENDOR=openmpi && export MPI_INTERCONNECT=ib && export MPI_COMPILER=gnu && mpicc $(CFLAGS) -o crumple crumple.c $(MPI_FLAGS) $(DEBUG_FLAGS)

mpi-profile:
	export MPI_VENDOR=openmpi && export MPI_INTERCONNECT=ib && export MPI_COMPILER=gnu && mpicc -mpilog $(CFLAGS) -o crumple crumple.c $(MPI_FLAGS)

test:
	@gcc $(CFLAGS) -o testes test.c $(DEBUG_FLAGS)
	@ ./testes 
	@rm testes

mpi-test:
	@mpicc $(CFLAGS) -o testes test.c $(MPI_FLAGS) $(DEBUG_FLAGS) # -D_GDB_ATTACH
	@ ./testes

14mer:
	@gcc $(CFLAGS) -o subopt subopt.c 
	echo "GCUCUAAAAGAGAG" | ./subopt  -noLP 
	./subopt -h14mer:

