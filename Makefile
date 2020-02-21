FC=/usr/bin/gfortran
#FFLAGS=-O2 -Wall -std=f95 -pedantic -frecord-marker=4

LIBS= -llapack -lblas -lfftw3

##########################################################

PREFIX=$(PWD)
MODDIR=$(PREFIX)/modules
SRCDIR=$(PREFIX)/src
OBJDIR=$(PREFIX)/objects
BINDIR=$(PREFIX)/bin

vpath %.x $(BINDIR)
vpath %.f90 $(SRCDIR)
vpath %.o $(OBJDIR)

INC=-I$(MODDIR) -J$(MODDIR)

##########################################################

SANDBOX = calcul_deplacement.x
MODULES = distance_pbc.o num_lines_file.o constants.o

all: $(SANDBOX)

.SUFFIXES:
.SUFFIXES: .f90 .x .o

%.x : %.o $(OBJECTS)
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(INC) $(OBJDIR)/$< $(addprefix $(OBJDIR)/,$(MODULES)) $(LIBS)

%.o : %.f90
	$(FC) $(FFLAGS) -c $(INC) $< -o $(OBJDIR)/$@

clean:
	rm -f $(SANDBOX) $(OBJDIR)/*.o $(MODDIR)/*.mod *~

calcul_deplacement.o: calcul_deplacement.f90 $(MODULES)
