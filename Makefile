SHELL := /bin/bash

CC = mpicc

EXECUTABLE=swellix.exe

EXECDIR = $(PWD)/
LIBDIR = $(EXECDIR)viennabuild/lib
INCLUDEDIR = $(EXECDIR)viennabuild/include/
BINDIR = $(EXECDIR)viennabuild/bin/
PARAMFILE = $(EXECDIR)viennabuild/rna_turner2004.par

VPATH=$(EXECDIR)src/:$(INCLUDEDIR):$(LIBDIR)

# Change this directory according to your filesystem. It doesn't necessarily need to have a large storage
# capacity. This is used to facilitate the bundling mechanism within Swellix if you choose to activate it.
BUNDLINGDIRECTORY=/scratch/nsloat/

FILES = $(wildcard $(EXECDIR)src/*.c)

OBJS = $(FILES:.c=.o)

OBJDIR=$(EXECDIR)src

export MPICC  = mpicc
MPIFLAGS = -D_MPI
#MPIRUN = /usr/mpi/gcc/mvapich2-2.1/bin/mpirun_rsh
#CC = GNU
export MPI_HARDWARE=ib
export MPI_SOFTWARE=openmpi
export MPI_COMPILER=intel
CFLAGS=-g -fPIC

RNACONF = --prefix=$(EXECDIR)viennabuild --without-kinfold --without-forester --without-kinwalker \
 	   --without-perl --without-python --without-doc --without-doc-html --without-doc-pdf


RNAFLAGS=-L$(LIBDIR) -lRNA -lgomp -lm -I$(INCLUDEDIR) -I$(INCLUDEDIR)ViennaRNA -D_PARAM='"$(PARAMFILE)"'

BUNDLINGFLAG=-D_BUNDLE='"$(BUNDLINGDIRECTORY)"'

VIENNADIRS=`ls -I rna\_turner2004\.par $(EXECDIR)viennabuild`

all:
	cd ViennaRNA-2.2.5; \
	./configure $(RNACONF); \
	make; \
	make install
	$(CC) -g -o exe--swellix $(CFILES) $(RNAFLAGS) -W -Wall -g3 -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG)

swellix.exe: $(OBJS)
	echo objdir = $(OBJDIR)
	$(CC) $(CFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -o $@ $^ 

$(OBJDIR)/main.o: main.c main.h init_general.h init_constraint.h paren_lookup_table.h component_list.h interval_lookup_table.h bundle_list.h jump_tree.h unit_tests.h close_up.h statistics.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<

$(OBJDIR)/bundle_list.o: bundle_list.c bundle_list.h main.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<

$(OBJDIR)/close_up.o: close_up.c close_up.h main.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<

$(OBJDIR)/component_list.o: component_list.c component_list.h main.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<

$(OBJDIR)/init_constraint.o: init_constraint.c init_constraint.h main.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<

$(OBJDIR)/init_general.o: init_general.c init_general.h main.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<

$(OBJDIR)/interval_lookup_table.o: interval_lookup_table.c interval_lookup_table.h main.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<

$(OBJDIR)/jump_tree.o: jump_tree.c jump_tree.h main.h close_up.h unit_tests.h statistics.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<

$(OBJDIR)/main_common.o: main_common.c main.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<

$(OBJDIR)/mpi_bundle_list.o: mpi_bundle_list.c mpi_bundle_list.h main.h bundle_list.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<

$(OBJDIR)/paren_lookup_table.o: paren_lookup_table.c paren_lookup_table.h main.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<

$(OBJDIR)/statistics.o: statistics.c statistics.h main.h
	$(CC) $(CFLAGS) $(MPIFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $< $(RNAFLAGS)

$(OBJDIR)/unit_tests.o: unit_tests.c unit_tests.h main.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<

$(OBJDIR)/subopt.o: subopt.c subopt.h
	$(CC) $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -c -o $@ $<


swellix-mpi.exe: $(OBJS)
	$(CC) $(CFLAGS) $(MPIFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -o $@ $^ $(RNAFLAGS)

vienna:
#	cd viennabuild; \
#	rm -r $(VIENNADIRS)
#	cd ViennaRNA-2.2.5; \
#	make clean; \
#	make distclean
	cd ViennaRNA-2.2.5; \
	./configure $(RNACONF); \
	make; \
	make install

clean:
	rm -f exe--swellix exe--swellix.mpi exe--swellix-tests
	cd ViennaRNA-2.2.5; \
	make clean; \
	make distclean
	cd viennabuild; \
	rm -r $(VIENNADIRS)

dispOn:
	gcc -g -o exe--swellix $(CFILES) $(RNAFLAGS) -W -Wall -g3 -D_EXECDIR='"$(EXECDIR)"' \
	-D _display -D _test $(BUNDLINGFLAG)

mpi_disp:
	$(MPICC) $(CFLAGS) -o exe--swellix.mpi $(CFILES) $(RNAFLAGS) -W -Wall -g3 \
	-D_MPI -D_display -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) -lX11 -lm

#mpi_prof:
#	mpicc ${CFLAGS} -o exe--swellix $(CFILES) -W -Wall -g3 -pg -lX11 -lm

#mpi_prof2:
#	mpicc-vt ${CFLAGS} -o exe--swellix $(CFILES) -W -Wall -g3 -pg -lX11 -lm

test:
	gcc -g -o exe--swellixTest $(CFILES) -W -Wall -g3 -D _display $(BUNDLINGFLAG)

profile0: 
	cc -h profile_generate -lm -o exe--swellix_craypath $(CFILES) $(RNAFLAGS) \
	-D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG)

profile1:
	cc -O3 -h pl=exe--swellix.pl -h wp -o exe--swellix_craypath $(CFILES) $(RNAFLAGS) \
	-D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG)

report:
	gcc -g -o exe--swellix $(CFILES) -W -Wall -g3 -pg -D _report $(BUNDLINGFLAG)

unittests:
	gcc -g -o exe--swellix-tests $(TFILES) -Wall -g3
