########################################
#    Master Makefile for EFISPEC3D     #
#                                      #
#  Please Edit ./config/make.config    #
#                                      #
#    Philippe THIERRY july 2013        #
########################################
export ROOTEFI=$${PWD}

include ./config/make.config

notarget:
	@echo ""
	@echo ""
	@echo "================================================================================================================================="
	@echo "                                               EFISPEC3D master Makefile : WARNING"
	@echo "================================================================================================================================="
	@echo ""
	@echo "      You cannot run \'make' without an explicit argument."
	@echo "      Read the installation instructions : ./docs/html/index.html"
	@echo "      Here is a summary :"
	@echo ""
	@echo ""
	@echo "      If you just extracted files from the archive EFISPEC3D.tgz"
	@echo "      To make EFISPEC3D :"
	@echo "      ================================================================================================="
	@echo "             run \'make all' to make all EFISPEC3D versions with Intel SSE, AVX and AVX2 instructions"
	@echo ""
	@echo "             or select one of the following builds:"
	@echo ""
	@echo "             make efispec version=version simd=intel_instruction"
	@echo "             with"
	@echo "             efi_version       = 0.9 or 1.0"
	@echo "             intel_instruction = sse or avx or avx2"
	@echo ""
	@echo "             ex: to make EFISPEC3D version 0.9 with AVX instruction"
	@echo "             run \'make efispec version=0.9 simd=avx' "
	@echo ""
	@echo ""
	@echo "      To make mesh partitioning tool :"
	@echo "      ================================================================================================="
	@echo "             run \'make mesh_partitioning'" to make mesh_partitioning.exe
	@echo ""
	@echo ""
	@echo "      To make EFISPEC3D tools :"
	@echo "      ================================================================================================="
	@echo "             run \'make efispec_tools'" to make some useful tools
	@echo ""
	@echo ""
	@echo "       To clean :"
	@echo "      ================================================================================================="
	@echo "             run \'make clean_all'                     to clean all libraries and executables"
	@echo "             run \'make clean_efispec version=  simd=' to clean selected EFISPEC3D version"
	@echo "             run \'make clean_mesh'                    to clean mesh partitioning tool"
	@echo "      ================================================================================================="
	@echo ""
	@echo ""

#################
clean_all:
	cd $(SDIR);  $(RM) */*.o */*/*.o  */a.out */core  */*.mod 
	cd $(BINDIR) ; rm -f *.exe
#######################
clean_efispec: 
	cd $(EFI_SRC);  $(RM) *.o  a.out core  *.mod ${EXEC}
clean_mesh: 
	cd $(EFI_MESH);  $(RM) *.o  a.out core  *.mod
##################################################
#### just do it 
efispec:
	@echo ""
	@cd $(EFI_SRC); $(MAKE) clean ; $(MAKE) all
	@echo ""
	@echo ""
	@echo ""
	@echo "         ********************************   "
	@echo "         *                              *   "
	@echo "         *  EFISPEC3D compilation done  *   "
	@echo "         *                              *   "
	@echo "         ********************************   "
	@echo ""
	@echo ""
	@echo ""
	@cd $(BINDIR); ls -lhtra
##################################################
#### just do it 
mesh_partitioning:
	@echo ""
	@cd $(EFI_MESH); $(MAKE) clean ; $(MAKE) all
	@echo ""
	@echo ""
	@echo ""
	@echo "         *********************************************   "
	@echo "         *                                           *   "
	@echo "         *  Mesh partitioning tool compilation done  *   "
	@echo "         *                                           *   "
	@echo "         *********************************************   "
	@echo ""
	@echo ""
	@echo ""
	@cd $(BINDIR); ls -ltra
##################################################

##################################################
#### just do it 
efispec_tools:
	@echo ""
	@cd $(EFI_TOOLS_CUBIT_TOPO)  ; $(MAKE) clean ; $(MAKE) all
	@cd $(EFI_TOOLS_CUBIT_REFINE); $(MAKE) clean ; $(MAKE) all
	@echo ""
	@echo ""
	@echo ""
	@echo "         **************************************   "
	@echo "         *                                    *   "
	@echo "         *  EFISPEC3D tools compilation done  *   "
	@echo "         *                                    *   "
	@echo "         **************************************   "
	@echo ""
	@echo ""
	@echo ""
	@cd $(BINDIR); ls -ltra
##################################################



##################################################
all: 
	for vect in sse avx avx2; do \
	for vers in efi_1.0 efi_0.9 ;  do \
		make efispec version=$$vers simd=$$vect ;\
	done	;\
	done
