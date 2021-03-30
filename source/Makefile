#! /bin/csh -f

# fortran pour compiler.

MYFOR = gfortran

# linker also fortran.
LD = ${MYFOR}

# flags for compilation
CFLAGS = -O3 -C -fdollar-ok -fbounds-check -fsignaling-nans -ffpe-trap=invalid,zero,overflow
#CFLAGS = -O0 -C -fdollar-ok -fbounds-check -fsignaling-nans -ffpe-trap=invalid,zero,overflow
#CFLAGS = -O0 -C -Wall -fdollar-ok -fbounds-check -fsignaling-nans -ffpe-trap=invalid,zero,overflow

# debug flag
DFLAGS = -g

# Linker flags 
LDFLAGS =


#suffixes
.SUFFIXES: .o .f90

# ----------------------------
# sources
SF90 =  \
        module_pdma.f90 \
        module_plotflex.f90\
        flex-hs.f90


# -------------------------------
# objets

OF90 = $(SF90:.f90=.o)


SOURCESFILES = $(SF90)
OBJECTSFILES = $(OF90)

# executable name
PROG = flex-hs.exe

# linking 
$(PROG): $(OBJECTSFILES)
	$(LD) $(LDFLAGS) -o $@ $(OBJECTSFILES) $(LIBS)


# compiling sources

.f90.o:
	$(MYFOR) -c $(CFLAGS) $(DFLAGS) $(INCL) $<


