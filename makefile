# Source Directory
#SRCDIR = src/

# Put compiler command in FC

#Intel MPI fortran: MPIIFORT compiler
# FC=mpiifort
#OpenMPI fortran compiler
FC = mpifort
# Put compile time options in FFLAGS


#MPIIFORT compiler flags
# FFLAGS = -Ofast -align -real-size 64 -mcmodel=large
#OpenMPI ocmpiler flags
FFLAGS = -mcmodel=large -fdefault-real-8 -Ofast

OBJ = main.o \
      derivatives4thOrder.o \
      energy.o \
      evolveeuler.o \
      fluxesNumRel.o \
      covariantDerivs.o \
      fieldStrengths.o \
      icMMbar.o \
      ghostsend.o \
      ghostrecvl.o \
      ghostrecvu.o \
      printenergy.o \
      evolveloop.o \
      processCoord.o \
      icgauge.o	\
      boundaryIndices.o \
      electromagnetics.o\

# Program name
PROGRAM = EW_dumbbell_relax.out

# Default make command is all
all: $(PROGRAM) 

# Create the program executable
$(PROGRAM): $(OBJ) 
	$(FC) -o  $@ $^ 
	rm *.o

# Object files (.o) are generated when .f files are compiled
%.o: %.f
	$(FC) $(FFLAGS) -c $*.f

# Clean all the junk
clean:
	rm: -f *.o *.dvi *.aux *.log core.* *~

# In this file, FC, FFLAGS, and PROGRAM are macros (variables) that can be used 
# throughout the file. If we just type make in the terminal, the command all in 
# this file be executed. $(PROGRAM) in all command is it's dependency and this 
# will cause $(PROGRAM) command to execute which, in turn, depends on the 
# object files $(OBJ). $(PROGRAM) commands is then executed. But $(OBJ) again 
# have dependencies on %.f files to be generated. So then the %.o command is 
# executed. Flag -c says to generate the object file, the -o $@ says to put 
# output of the compilation in the file named on the left side of the :, $< is 
# the first term in the dependency list. $@ and $^ are left and right sides of 
# :, respectively.
