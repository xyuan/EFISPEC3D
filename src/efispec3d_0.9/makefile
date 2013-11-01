## Main config file
include ../../config/make.config

## Objects and sources
FOBJS	= $(FSRCS:.f90=.o)
FSRCS	= $(wildcard *.f90)

## Fortran module
MODSRCS	= module_global_variables.o module_init_memory.o module_lagrange.o module_jacobian.o module_init_efi.o module_init_mesh.o module_source_function.o module_vtk_io.o module_init_medium.o\
          module_coordinate.o module_gll_value.o module_init_mpibuffer.o module_solver.o module_receiver.o module_source.o module_snapshot.o module_mpi_assemble.o

all:	$(MODSRCS) $(FOBJS) 
	$(FC) $(LDFLAGS) $(FOBJS) $(LIBS)  -o $(EXEC)  

clean:
	$(RM) $(FOBJS) $(EXEC) *.mod *.o
