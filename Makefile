SHELL := /bin/bash

myCC =gcc

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

CFILES = $(wildcard $(EXECDIR)src/*.c)

OBJS = $(CFILES:.c=.o)

OBJDIR=$(EXECDIR)src

#export MPICC  = mpic
MPIFLAGS =
ifeq (mpicc, $(myCC))
MPIFLAGS +=-D_MPI
endif

#CC = GNU
#export MPI_HARDWARE=ib
#export MPI_SOFTWARE=openmpi
#export MPI_COMPILER=intel
CFLAGS=-g -fPIC# -fopenmp

RNACONF = --prefix=$(EXECDIR)viennabuild --without-kinfold --without-forester --without-kinwalker \
 	   --without-perl --without-python --without-doc --without-doc-html --without-doc-pdf


RNAFLAGS=-L$(LIBDIR) -lRNA -lm -I$(INCLUDEDIR) -I$(INCLUDEDIR)ViennaRNA -D_PARAM='"$(PARAMFILE)"'

BUNDLINGFLAG=-D_BUNDLE='"$(BUNDLINGDIRECTORY)"'

VIENNADIRS=`ls -I rna\_turner2004\.par $(EXECDIR)viennabuild`

.PHONY: serial disp mpi mpi-disp pre-build-mpi pre-build

serial: pre-build
serial: swellix.exe

swellix.exe: $(OBJS)
	@echo myCC = $(shell which $(myCC))
	$(myCC) -o $@ $^ $(CFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG)

$(OBJDIR)/main.o: main.c main.h init_general.h init_constraint.h paren_lookup_table.h component_list.h interval_lookup_table.h bundle_list.h jump_tree.h unit_tests.h close_up.h statistics.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) 

$(OBJDIR)/bundle_list.o: bundle_list.c bundle_list.h main.h subopt.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) 

$(OBJDIR)/close_up.o: close_up.c close_up.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) 

$(OBJDIR)/component_list.o: component_list.c component_list.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) 

$(OBJDIR)/init_constraint.o: init_constraint.c init_constraint.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) 

$(OBJDIR)/init_general.o: init_general.c init_general.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) 

$(OBJDIR)/interval_lookup_table.o: interval_lookup_table.c interval_lookup_table.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) 

$(OBJDIR)/jump_tree.o: jump_tree.c jump_tree.h main.h close_up.h unit_tests.h statistics.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) 

$(OBJDIR)/main_common.o: main_common.c main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) 

$(OBJDIR)/mpi_bundle_list.o: mpi_bundle_list.c mpi_bundle_list.h main.h bundle_list.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) 

$(OBJDIR)/paren_lookup_table.o: paren_lookup_table.c paren_lookup_table.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG)

$(OBJDIR)/statistics.o: statistics.c statistics.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG) 

$(OBJDIR)/unit_tests.o: unit_tests.c unit_tests.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG)

$(OBJDIR)/subopt.o: subopt.c subopt.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG)

swellix-mpi.exe: $(OBJS)
	@echo myCC = $(shell which $(myCC))
	$(myCC) -o $@ $^ $(CFLAGS) $(RNAFLAGS) $(MPIFLAGS) -D_EXECDIR='"$(EXECDIR)"' $(BUNDLINGFLAG)

pre-build-mpi:

pre-build:

clean:
	rm $(OBJS)

disp: CFLAGS+=-D_display
disp: myCC=gcc
disp: pre-build
disp: swellix.exe

mpi: myCC=mpicc
mpi: MPIFLAGS+=-D_MPI
mpi: pre-build-mpi
mpi: swellix-mpi.exe

mpi-disp: myCC=mpicc
mpi-disp: MPIFLAGS+=-D_MPI
mpi-disp: CFLAGS+=-D_display
mpi-disp: pre-build-mpi
mpi-disp: swellix-mpi.exe

vienna:
	cd ViennaRNA-2.2.5; \
	./configure $(RNACONF); \
	make; \
	make install

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

