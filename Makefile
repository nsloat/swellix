SHELL := /bin/bash

myCC=gcc

EXECDIR = $(PWD)/
LIBDIR = $(EXECDIR)viennabuild/lib
INCLUDEDIR = $(EXECDIR)viennabuild/include/
BINDIR = $(EXECDIR)viennabuild/bin/
PARAMFILE = $(EXECDIR)viennabuild/rna_turner2004.par

VPATH=$(EXECDIR)src/:$(INCLUDEDIR):$(LIBDIR)

CFILES = $(wildcard $(EXECDIR)src/*.c)

OBJS = $(CFILES:.c=.o)

OBJDIR=$(EXECDIR)src

CFLAGS=-g

RNACONF = --prefix=$(EXECDIR)viennabuild --without-kinfold --without-forester --without-kinwalker \
 	   --without-perl --without-python --without-doc --without-doc-html --without-doc-pdf --disable-openmp

RNAFLAGS=-L$(LIBDIR) -lRNA -lm -I$(INCLUDEDIR) -I$(INCLUDEDIR)ViennaRNA -D_PARAM='"$(PARAMFILE)"'

.PHONY: serial ddisp cdisp disp mpi mpi-disp pre-build-mpi pre-build

serial: pre-build
serial: swellix.exe

swellix.exe: $(OBJS)
	@echo myCC = $(shell which $(myCC))
	$(myCC) -o $@ $^ $(CFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' 

$(OBJDIR)/main.o: main.c main.h init_general.h init_constraint.h paren_lookup_table.h component_list.h interval_lookup_table.h bundle_list.h jump_tree.h unit_tests.h close_up.h statistics.h mpi_jump_tree.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' 

$(OBJDIR)/bundle_list.o: bundle_list.c bundle_list.h main.h subopt.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"'  

$(OBJDIR)/close_up.o: close_up.c close_up.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"'  

$(OBJDIR)/component_list.o: component_list.c component_list.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"'  

$(OBJDIR)/init_constraint.o: init_constraint.c init_constraint.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"'  

$(OBJDIR)/init_general.o: init_general.c init_general.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"'  

$(OBJDIR)/interval_lookup_table.o: interval_lookup_table.c interval_lookup_table.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"'  

$(OBJDIR)/jump_tree.o: jump_tree.c jump_tree.h main.h close_up.h unit_tests.h statistics.h mpi_jump_tree.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"'  

$(OBJDIR)/main_common.o: main_common.c main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"'  

$(OBJDIR)/mpi_jump_tree.o: mpi_jump_tree.c mpi_jump_tree.h main.h jump_tree.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"'  

$(OBJDIR)/paren_lookup_table.o: paren_lookup_table.c paren_lookup_table.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' 

$(OBJDIR)/statistics.o: statistics.c statistics.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"'  

$(OBJDIR)/unit_tests.o: unit_tests.c unit_tests.h main.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' 

$(OBJDIR)/subopt.o: subopt.c subopt.h
	@echo myCC = $(shell which $(myCC))
	$(myCC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) $(RNAFLAGS) -D_EXECDIR='"$(EXECDIR)"' 

swellix-mpi.exe: $(OBJS)
	@echo myCC = $(shell which $(myCC))
	$(myCC) -o $@ $^ $(CFLAGS) $(RNAFLAGS) $(MPIFLAGS) -D_EXECDIR='"$(EXECDIR)"' 

pre-build-mpi:

pre-build:

clean:
	rm $(OBJS)

cdisp: clean
cdisp: disp

ddisp: CFLAGS+=-g 
ddisp: myCC=cc
ddisp: disp

disp: CFLAGS+=-D_display
disp: myCC=cc
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
	make install;

test:
	cc -g -o exe--swellixTest $(CFILES) -W -Wall -g3 -D _display 

profile0: 
	cc -h profile_generate -lm -o exe--swellix_craypath $(CFILES) $(RNAFLAGS) \
	-D_EXECDIR='"$(EXECDIR)"' 

profile1:
	cc -O3 -h pl=exe--swellix.pl -h wp -o exe--swellix_craypath $(CFILES) $(RNAFLAGS) \
	-D_EXECDIR='"$(EXECDIR)"' 

report:
	cc -g -o exe--swellix $(CFILES) -W -Wall -g3 -pg -D _report 

unittests:
	cc -g -o exe--swellix-tests $(TFILES) -Wall -g3

