## Main config file
include ../../../config/make.config

## Objects and sources
FOBJS	= $(FSRCS:.f90=.o)
FSRCS	= $(wildcard *.f90)

## Fortran module

all:	$(FOBJS) 
	$(FC) $(LDFLAGS) $(FOBJS) -o $(EXEC_TOOLS_CUBIT_REFINE)  

clean:
	$(RM) $(FOBJS) $(EXEC_TOOLS_CUBIT_REFINE) *.mod *.o
