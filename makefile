F90 = gfortran
LLAPACK  = -lblas -llapack
DEFS = -DDP 

OPTFLAGS = -fall-intrinsics #-O3 #-ffpe-trap=inexact
FFLAGS  = -cpp $(DEFS) $(OPTFLAGS) $(MPFUN) $(MKLINCLUDE)
LDFLAGS = $(OPTFLAGS) $(LLAPACK)

###### Section for developers

.SUFFIXES: 
.SUFFIXES: .o .f90 .mod

## main modules and objects ##
EXE  = mpol
OBJC = typedefs.o consts.o lapackmodule.o inverse.o systm.o output.o main.o 
MODF = typedefs.mod consts.mod lapackmodule.mod inverse.mod systm.mod output.mod 


.f90.mod:
	$(F90) -c $(FFLAGS) $<

.f90.o:
	$(F90) -c $(FFLAGS) $<

mpol: $(MODC) $(OBJC)
	$(F90) -o $@ $(OBJC) $(LDFLAGS)

clean:
	rm -f *.o *.mod *~
 
