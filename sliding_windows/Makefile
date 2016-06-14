CFLAGS=-Wall \
	-g
DEBUG_FLAGS=-g3 -pg
MPI_FLAGS=-march=core2 -mtune=core2 -D_MPI

all:
	gcc $(CFLAGS) -o subopt_and_label subopt.c 

debug:
	gcc $(CFLAGS) -o subopt_and_label subopt.c $(DEBUG_FLAGS)

profile:
	gcc $(CFLAGS) -o subopt_and_label subopt.c -pg

mpi:
	export MPI_VENDOR=openmpi
	export MPI_INTERCONNECT=ib
	export MPI_COMPILER=gcc
	mpicc $(CFLAGS) -o subopt_and_label subopt.c $(MPI_FLAGS) $(DEBUG_FLAGS) # -D_GDB_ATTACH

test:
	@gcc -o subopt_and_label subopt.c -lRNA -lm -I../pvienna/H/
	@echo "GCUCUAAAAGAGAG" | ./subopt_and_label > test
	@wc -l test | awk '{ if ($$1 == 120) print "PASSED";}'
	@rm test

14mer:
	@gcc -o subopt_and_label subopt.c -lRNA -lm -I../pvienna/H/
	echo "GCUCUAAAAGAGAG" | ./subopt_and_label  -noLP 
	./subopt -h14mer:

gauntlet:
	@gcc -o subopt_and_label subopt.c -lRNA -lm -I../pvienna/H/
	echo "gcgcuuugcgc" | ./subopt_and_label  -e 5  > ours  
	echo "gcgcuuugcgc" | RNAsubopt -e 5  > thiers
	echo "auacccuau"   | ./subopt_and_label  -e 5 >> ours
	echo "auacccuau"   | RNAsubopt -e 5 >> thiers
	echo "cauacccuau"  | ./subopt_and_label  -e 5 >> ours
	echo "cauacccuau"  | RNAsubopt -e 5 >> thiers
	diff -y ours thiers

