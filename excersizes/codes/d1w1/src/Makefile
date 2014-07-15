/SHELL = /bin/tcsh

# The C compiler
CC = gcc

# The Fortran compiler
FC = gfortran

# The executable goes here
BINDIR = ~/Talent_ECT2014_repo.git/excersizes/codes/d1w1/bin

# Name of executable
EXENAME = matel_test

LIBS = -L$/. -llapack -lblas

OPTS = -O3  -funroll-loops -ffast-math -g -fcheck=all -fbounds-check -ffree-form -ffree-line-length-1024 -fbacktrace -ffpe-trap=invalid -fdefault-double-8 -fdefault-real-8

# Notice that use statements may only invoke modules lying further upp (to the left) 
# in the list
OBJECTS = types.o symmat_mod.o matrix_elements_test.o

#.f90.o :
#	$(FC) -c $(OPTS) $<

$(BINDIR)/$(EXENAME) : $(OBJECTS)
	$(FC) $(OBJECTS) $(LIBS) -o $@

types.o : types.f90
	$(FC) -c $(OPTS) $<

matrix_elements_test.o : matrix_elements_test.f90
	$(FC) -c $(OPTS) $<

symmat_mod.o : symmat_mod.f90
	$(FC) -c $(OPTS) $<



clean:
	rm -f $(OBJECTS)
	rm -f *.mod