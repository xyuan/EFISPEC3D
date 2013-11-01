!=====================================================================================================================================
!                                      EFISPEC3D                                             !
!                            (Elements FInis SPECtraux 3D)                                   !
!                                                                                            !
!                               http://efispec.free.fr                                       !
!                                                                                            !
!                                                                                            !
!                            This file is part of EFISPEC3D                                  !
!              Please refer to http://efispec.free.fr if you use it or part of it            !
!                                                                                            !
!                                                                                            !
!1 ---> French License: CeCILL V2                                                            !
!                                                                                            !
!         Copyright BRGM 2009 "contributeurs : Florent  DE MARTIN (BRGM)                     !
!                                              David    MICHEA    (BRGM)                     !
!                                              Philippe THIERRY   (Intel)"                   !
!                                                                                            !
!         Contact: f.demartin at brgm.fr                                                     !
!                                                                                            !
!         Ce logiciel est un programme informatique servant a resoudre l'equation du         !
!         mouvement en trois dimensions via une methode des elements finis spectraux.        !
!                                                                                            !
!         Ce logiciel est regi par la licence CeCILL soumise au droit francais et            !
!         respectant les principes de diffusion des logiciels libres. Vous pouvez            !
!         utiliser, modifier et/ou redistribuer ce programme sous les conditions de la       !
!         licence CeCILL telle que diffusee par le CEA, le CNRS et l'INRIA sur le site       !
!         "http://www.cecill.info".                                                          !
!                                                                                            !
!         En contrepartie de l'accessibilite au code source et des droits de copie, de       !
!         modification et de redistribution accordes par cette licence, il n'est offert      !
!         aux utilisateurs qu'une garantie limitee. Pour les memes raisons, seule une        !
!         responsabilite restreinte pese sur l'auteur du programme, le titulaire des         !
!         droits patrimoniaux et les concedants successifs.                                  !
!                                                                                            !
!         A cet egard l'attention de l'utilisateur est attiree sur les risques associes      !
!         au chargement, a l'utilisation, a la modification et/ou au developpement et a      !
!         la reproduction du logiciel par l'utilisateur etant donne sa specificite de        !
!         logiciel libre, qui peut le rendre complexe a manipuler et qui le reserve donc     !
!         a des developpeurs et des professionnels avertis possedant des connaissances       !
!         informatiques approfondies. Les utilisateurs sont donc invites a charger et        !
!         tester l'adequation du logiciel a leurs besoins dans des conditions permettant     !
!         d'assurer la securite de leurs systemes et ou de leurs donnees et, plus            !
!         generalement, a l'utiliser et l'exploiter dans les memes conditions de             !
!         securite.                                                                          !
!                                                                                            !
!         Le fait que vous puissiez acceder a cet en-tete signifie que vous avez pris        !
!         connaissance de la licence CeCILL et que vous en avez accepte les termes.          !
!                                                                                            !
!                                                                                            !
!2 ---> International license: GNU GPL V3                                                    !
!                                                                                            !
!         EFISPEC3D is a computer program that solves the three-dimensional equations of     !
!         motion using a finite spectral-element method.                                     !
!                                                                                            !
!         Copyright (C) 2009 Florent DE MARTIN                                               !
!                                                                                            !
!         Contact: f.demartin at brgm.fr                                                     !
!                                                                                            !
!         This program is free software: you can redistribute it and/or modify it under      !
!         the terms of the GNU General Public License as published by the Free Software      !
!         Foundation, either version 3 of the License, or (at your option) any later         !
!         version.                                                                           !
!                                                                                            !
!         This program is distributed in the hope that it will be useful, but WITHOUT ANY    !
!         WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    !
!         PARTICULAR PURPOSE. See the GNU General Public License for more details.           !
!                                                                                            !
!         You should have received a copy of the GNU General Public License along with       !
!         this program. If not, see http://www.gnu.org/licenses/.                            !
!                                                                                            !
!                                                                                            !
!3 ---> Third party libraries                                                                !
!                                                                                            !
!         EFISPEC3D uses the following of third party libraries:                             !
!                                                                                            !
!           --> METIS 5.1.0 under the Apache License Version 2.0 (compatible with GNU GPL V3)!
!               see http://glaros.dtc.umn.edu/gkhome/metis/metis/overview                    !
!                                                                                            !
!           --> Lib_VTK_IO under GNU GPL V3 License                                          !
!               see S. Zaghi's website: https://github.com/szaghi/Lib_VTK_IO                 !
!                                                                                            !
!           --> INTERP_LINEAR under GNU GPL License                                          !
!               see J. Burkardt website: http://people.sc.fsu.edu/~jburkardt/                !
!                                                                                            !
!                                                                                            !
!4 ---> Related Articles                                                                     !
!                                                                                            !
!         De Martin, F., Matsushima, M., Kawase, H. (BSSA, 2013)                             !
!            Impact of geometric effects on near-surface Green's functions                   !
!            doi:10.1785/0120130039                                                          !
!                                                                                            !
!         Aochi, H., Ducellier, A., Dupros, F., Delatre, M., Ulrich, T., De Martin, F.,      !
!         Yoshimi, M., (Pure Appl. Geophys. 2013)                                            !
!            Finite Difference Simulations of Seismic Wave Propagation for the 2007 Mw 6.6   !
!            Niigata-ken Chuetsu-Oki Earthquake: Validity of Models and Reliable             !
!            Input Ground Motion in the Near-Field                                           !
!            doi:10.1007/s00024-011-0429-5                                                   !
!                                                                                            !
!         De Martin, F. (BSSA, 2011)                                                         !
!            Verification of a Spectral-Element Method Code for the Southern California      !
!            Earthquake Center LOH.3 Viscoelastic Case                                       !
!            doi:10.1785/0120100305                                                          !
!                                                                                            !
!=====================================================================================================================================

!>@file
!!This file contains the module to define all global variables of EFISPEC3D.

!>@brief
!!This module defines all global variables of EFISPEC3D.
!!Scalar variables are initialized directly in this module.
!!This module also contains some general purpose subroutines.
module mod_global_variables
   
   use mpi
   
   implicit none

!
!
!***********************************************************************************************************************************************************************************
!parameters
!***********************************************************************************************************************************************************************************

!>if .true. --> activate viscoelasticity in the solver module: see mod_solver
   logical, parameter :: LG_VISCO = .false.

!>if .true. --> write snapshots in VTK format
   logical, parameter :: LG_SNAPSHOT_VTK = .true.

!>if .true. --> write snapshots in gmt format
   logical, parameter :: LG_SNAPSHOT_GMT = .false.

!>if .true. --> write medium in VTK format
   logical, parameter :: LG_OUTPUT_MEDIUM_VTK = .true.

!>order of Lagrange's polynomials. Should be set from 4 to 6.
   integer, parameter :: IG_LAGRANGE_ORDER = 4

!>number of GLL nodes in the reference domain [-1:1]
   integer, parameter :: IG_NGLL = IG_LAGRANGE_ORDER + 1

!>unit for listing file *.lst
   integer, parameter :: IG_LST_UNIT = 10

!>number of degrees of freedom of the problem (i.e., 3)
   integer, parameter :: IG_NDOF = 3

!>number of relaxation times and weights for viscoelastic stress-strain law (following liu and archuleta, 2006)
   integer, parameter :: IG_NRELAX = 8

!>RG_NEWMARK_GAMMA : gamma coefficient of Newmark method                                                                                                         
   real   , parameter :: RG_NEWMARK_GAMMA = 0.5

!>RG_PI : approximation of pi
   real   , parameter :: RG_PI = 3.141592654

!>if .true. --> activate asynchrone MPI communications between cpu
   logical, parameter :: LG_ASYNC_MPI_COMM = .false.

!>if .true. --> output cpu_time information in file *.time.cpu.*
   logical, parameter :: LG_OUTPUT_CPUTIME = .false.

!>if .true. --> output debug.* files
   logical, parameter :: LG_OUTPUT_DEBUG_FILE = .false.

!
!
!***********************************************************************************************************************************************************************************
!variables (scalars, or n-order tensors)
!***********************************************************************************************************************************************************************************

!>distributed forces at GLL nodes to generate double couple point source in cpu myrank
   real, dimension(:,:,:,:,:), allocatable :: rg_dcsource_gll_force

!>jacobian determinant at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_jacobian_det

!>derivative of local coordinate (@f$\xi@f$) with respect to global coordinate @f$x@f$ at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_dxidx

!>derivative of local coordinate (@f$\xi@f$) with respect to global coordinate @f$y@f$ at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_dxidy

!>derivative of local coordinate (@f$\xi@f$) with respect to global coordinate @f$z@f$ at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_dxidz

!>derivative of local coordinate (@f$\eta@f$) with respect to global coordinate @f$x@f$ at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_detdx

!>derivative of local coordinate (@f$\eta@f$) with respect to global coordinate @f$y@f$ at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_detdy

!>derivative of local coordinate (@f$\eta@f$) with respect to global coordinate @f$z@f$ at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_detdz

!>derivative of local coordinate (@f$\zeta@f$) with respect to global coordinate @f$x@f$ at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_dzedx

!>derivative of local coordinate (@f$\zeta@f$) with respect to global coordinate @f$y@f$ at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_dzedy

!>derivative of local coordinate (@f$\zeta@f$) with respect to global coordinate @f$z@f$ at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_dzedz

!>density (@f$ \rho @f$) at GLL nodes for hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_rho

!>density times squared shear-wave velocity (@f$\rho \beta^{2}@f$) at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_rhovs2

!>density times squared pressure-wave velocity (@f$\rho \alpha^{2}@f$) at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_hexa_gll_rhovp2

!>density times shear-wave velocity (@f$\rho \beta@f$) at GLL nodes of absorbing quadrangle elements in cpu myrank
   real, dimension(:,:,:), allocatable :: rg_quadp_gll_rhovs

!>density times shear-wave velocity (@f$\rho \alpha@f$) at GLL nodes of absorbing quadrangle elements in cpu myrank
   real, dimension(:,:,:), allocatable :: rg_quadp_gll_rhovp

!>jacobian determinant at GLL nodes of quadrangle elements in cpu myrank
   real, dimension(:,:,:), allocatable :: rg_quadp_gll_jaco_det

!>vector normal to GLL nodes of quadrangle elements in x,y,z coordinate system in cpu myrank
   real, dimension(:,:,:,:), allocatable :: rg_quadp_gll_normal

!>weights @f$ w^{S}_{k} @f$ for a given quality factor @f$ Q_{S} @f$ at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:,:), allocatable :: rg_hexa_gll_wkqs

!>weights @f$ w^{P}_{k} @f$ for a given quality factor @f$ Q_{P} @f$ at GLL nodes of hexahedron elements in cpu myrank
   real, dimension(:,:,:,:,:), allocatable :: rg_hexa_gll_wkqp

!>memory variables @f$ ksi_{xx} @f$ at GLL nodes of hexahedron elements in cpu myrank used for viscoelastic simulation 
   real, dimension(:,:,:,:,:), allocatable :: rg_hexa_gll_ksixx

!>memory variables @f$ ksi_{yy} @f$ at GLL nodes of hexahedron elements in cpu myrank used for viscoelastic simulation 
   real, dimension(:,:,:,:,:), allocatable :: rg_hexa_gll_ksiyy

!>memory variables @f$ ksi_{zz} @f$ at GLL nodes of hexahedron elements in cpu myrank used for viscoelastic simulation 
   real, dimension(:,:,:,:,:), allocatable :: rg_hexa_gll_ksizz

!>memory variables @f$ ksi_{xy} @f$ at GLL nodes of hexahedron elements in cpu myrank used for viscoelastic simulation 
   real, dimension(:,:,:,:,:), allocatable :: rg_hexa_gll_ksixy

!>memory variables @f$ ksi_{xz} @f$ at GLL nodes of hexahedron elements in cpu myrank used for viscoelastic simulation 
   real, dimension(:,:,:,:,:), allocatable :: rg_hexa_gll_ksixz

!>memory variables @f$ ksi_{yz} @f$ at GLL nodes of hexahedron elements in cpu myrank used for viscoelastic simulation 
   real, dimension(:,:,:,:,:), allocatable :: rg_hexa_gll_ksiyz

!>coefficients for viscoelastic simulation
   real, dimension(IG_NRELAX,3), parameter :: RG_RELAX_COEFF = [             &
                                                                 1.72333e-3, & 
                                                                 1.80701e-3, &
                                                                 5.38887e-3, &
                                                                 1.99322e-2, &
                                                                 8.49833e-2, &
                                                                 4.09335e-1, &
                                                                 2.05951e+0, &
                                                                13.26290e+0, &
                                                                 1.66958e-2, &
                                                                 3.81644e-2, &
                                                                 9.84666e-3, &
                                                                -1.36803e-2, &
                                                                -2.85125e-2, &
                                                                -5.37309e-2, &
                                                                -6.65035e-2, &
                                                                -1.33696e-1, &
                                                                 8.98758e-2, &
                                                                 6.84635e-2, &
                                                                 9.67052e-2, &
                                                                 1.20172e-1, &
                                                                 1.30728e-1, &
                                                                 1.38746e-1, &
                                                                 1.40705e-1, &
                                                                 2.14647e-1  &
                                                                ]

!>exponential used for viscoelastic simulation
   real, dimension(IG_NRELAX) :: rg_mem_var_exp

!>@f$x,y,z@f$ displacement at step n or n+1 at the GLL nodes (global numbering) in cpu myrank
   real, dimension(:,:), allocatable :: rg_gll_displacement

!>@f$x,y,z@f$ velocity at step n or n+1 at the GLL nodes (global numbering) in cpu myrank
   real, dimension(:,:), allocatable :: rg_gll_velocity

!>@f$x,y,z@f$ acceleration or external forces at step n+1 at the GLL nodes (global numbering) in cpu myrank
   real, dimension(:,:), allocatable :: rg_gll_acceleration

!>@f$x,y,z@f$ acceleration at step n at the GLL nodes (global numbering) in cpu myrank
   real, dimension(:,:), allocatable :: rg_gll_acctmp

!>local coordinate of geometric nodes in the reference domain [-1:1] 
   real, dimension(:), allocatable :: rg_gnode_abscissa

!>denominator for computing lagrange polynomial for geometric nodes in the reference domain [-1:1]
   real, dimension(:,:), allocatable :: rg_gnode_abscissa_dist

!>diagonal mass matrix at GLL nodes (global numbering) in cpu myrank
   real, dimension(:), allocatable :: rg_gll_mass_matrix

!>@f$x@f$ coordinate of the geometric nodes in cpu myrank
   real, dimension(:), allocatable :: rg_gnode_x

!>@f$y@f$ coordinate of the geometric nodes in cpu myrank
   real, dimension(:), allocatable :: rg_gnode_y

!>@f$z@f$ coordinate of the geometric nodes in cpu myrank
   real, dimension(:), allocatable :: rg_gnode_z

!>@f$z@f$ coordinate of all receivers used for snapshot (only in cpu0)
   real, dimension(:), allocatable :: rg_receiver_snapshot_z

!>local number of receiver used for snapshot
   integer, dimension(:), allocatable :: ig_receiver_snapshot_locnum

!>user defined source function for double couple source
   real, dimension(:), allocatable :: rg_dcsource_user_func

!>user defined source function for single force point source
   real, dimension(:), allocatable :: rg_sfsource_user_func

!>send buffer for asynchrone communications
   real, dimension(:), allocatable :: rg_mpi_buffer_send

!>receive buffer for asynchrone communications
   real, dimension(:), allocatable :: rg_mpi_buffer_recv

!>time step increment @f$ \Delta t @f$
   real :: rg_dt = 0.0

!>squared time step increment @f$ \Delta t^{2} @f$
   real :: rg_dt2 = 0.0

!>current time of simulation @f$ t_{simu} @f$ (i.e. rg_dt*ig_idt)
   real :: rg_simu_current_time = 0.0

!>total time of simulation @f$ tt_{simu} @f$
   real :: rg_simu_total_time = 0.0

!>maximum @f$x@f$ coordinate of the physical domain of study (over all cpus)
   real :: rg_mesh_xmax = 0.0

!>minimum @f$x@f$ coordinate of the physical domain of study (over all cpus) 
   real :: rg_mesh_xmin = 0.0

!>maximum @f$y@f$ coordinate of the physical domain of study (over all cpus)
   real :: rg_mesh_ymax = 0.0

!>minimum @f$y@f$ coordinate of the physical domain of study (over all cpus)
   real :: rg_mesh_ymin = 0.0

!>maximum @f$z@f$ coordinate of the physical domain of study (over all cpus)
   real :: rg_mesh_zmax = 0.0

!>minimum @f$z@f$ coordinate of the physical domain of study (over all cpus)
   real :: rg_mesh_zmin = 0.0

!>@f$dx@f$ and @f$dy@f$ of receivers for snapshot
   real :: rg_receiver_snapshot_dxdy = 0.0

!>number of receivers for snapshot in the @f$x@f$-direction
   integer :: rg_receiver_snapshot_nx = 0

!>number of receivers for snapshot in the @f$y@f$-direction
   integer :: rg_receiver_snapshot_ny = 0

!>global GLL nodes of hexahedron elements from local gll nodes
   integer, target, dimension(:,:,:,:), allocatable :: ig_hexa_gll_glonum

!>global GLL nodes of paraxial quadrangle elements from local gll nodes
   integer, dimension(:,:,:), allocatable :: ig_quadp_gll_glonum

!>global GLL nodes of free surface quadrangle elements from local gll nodes
   integer, dimension(:,:,:), allocatable :: ig_quadf_gll_glonum

!>global geometric nodes of hexahedron elements in cpu myrank 
   integer, dimension(:,:), allocatable :: ig_hexa_gnode_glonum

!>global geometric nodes of paraxial quadrangle elements in cpu myrank
   integer, dimension(:,:), allocatable :: ig_quadp_gnode_glonum

!>global geometric nodes of free surface quadrangle elements in cpu myrank
   integer, dimension(:,:), allocatable :: ig_quadf_gnode_glonum

!>local position @f$ \xi @f$ of geometric nodes in hexahedron elements
   integer, dimension(:), allocatable :: ig_hexa_gnode_xiloc

!>local position @f$ \eta @f$ of geometric nodes in hexahedron elements
   integer, dimension(:), allocatable :: ig_hexa_gnode_etloc

!>local position @f$ \zeta @f$ of geometric nodes in hexahedron elements
   integer, dimension(:), allocatable :: ig_hexa_gnode_zeloc

!>local position @f$ \xi @f$ of geometric nodes in quadrangle elements
   integer, dimension(:), allocatable :: ig_quad_gnode_xiloc

!>local position @f$ \eta @f$ of geometric nodes in quadrangle elements
   integer, dimension(:), allocatable :: ig_quad_gnode_etloc

!>fortran unit number of receivers in hexahedron elements in cpu myrank
   integer, dimension(:), allocatable :: ig_hexa_receiver_unit

!>fortran unit number of receivers in quadrangle elements in cpu myrank
   integer, dimension(:), allocatable :: ig_quad_receiver_unit

!>maximum number of GLL nodes in all mpi buffers for cpu myrank
   integer :: ig_mpi_buffer_sizemax = 0

!>mpi send request for asynchrone communications
   integer, dimension(:), allocatable :: ig_mpi_request_send

!>mpi receive request for asynchrone communications
   integer, dimension(:), allocatable :: ig_mpi_request_recv

!>offset buffer for asynchrone communications
   integer, dimension(:), allocatable :: ig_mpi_buffer_offset 

!>material number of hexahedron elements
   integer, dimension(:), allocatable :: ig_hexa_material_number

!>material type: 1=elastic, 2=viscoelastic
   integer(kind=1) , dimension(:), allocatable :: ig_material_type

!>global numbering of all receivers of all cpus to write *.grd files (gathered by cpu0 only)
   integer, dimension(:), allocatable :: ig_receiver_snapshot_glonum

!>shift for gatherv for snapshot receivers
   integer, dimension(:), allocatable :: ig_receiver_snapshot_mpi_shift

!>number of receivers in each cpu (gathered by cpu0 only)
   integer, dimension(:), allocatable :: ig_receiver_snapshot_total_number

!>total number of cpu used for computation
   integer :: ig_ncpu = 0

!>rank of the cpu myrank 
   integer :: ig_myrank = 0

!>number of cpus neighbor to cpu myrank
   integer :: ig_ncpu_neighbor = 0

!>info about cpus neighbor to cpu myrank: ig_cpu_neighbor_info(1,:) --> element number, ig_cpu_neighbor_info(2,:) --> face number, ig_cpu_neighbor_info(3,:) --> cpu number
   integer, dimension(:,:), allocatable :: ig_cpu_neighbor_info

!>number of hexahedron elements in cpu myrank
   integer :: ig_nhexa = 0

!>number of outer hexahedron elements in cpu myrank
   integer :: ig_nhexa_outer = 0

!>number of inner hexahedron elements in cpu myrank
   integer :: ig_nhexa_inner = 0

!>number of paraxial quadrangle elements in cpu myrank
   integer :: ig_nquad_parax = 0

!>number of free surface quadrangle elements in cpu myrank
   integer :: ig_nquad_fsurf = 0

!>total number of geometric nodes in cpu myrank
   integer :: ig_mesh_nnode = 0

!>number of geometric nodes to define a hexahedron element
   integer :: ig_hexa_nnode = 0

!>number of geometric nodes to define a quadrangle element
   integer :: ig_quad_nnode = 0

!>number of geometric nodes to define a line element
   integer :: ig_line_nnode = 0

!>indirection from hexahedron geometrical nodes to local GLL nodes (only for vertex nodes)
   integer, dimension(3,8) :: ig_hexa_node2gll

!>number of time step 
   integer :: ig_ndt = 0

!>current time step of analysis
   integer :: ig_idt = 0

!>saving increment for receivers time history
   integer :: ig_receiver_saving_incr = 0

!>saving increment for free surface snapshot
   integer :: ig_snapshot_saving_incr = 0 

!>if .true. --> boundary absorption is activated
   logical :: lg_boundary_absorption = .true.

!>number of double couple point sources in cpu myrank 
   integer :: ig_ndcsource = 0

!>number of single force point sources in cpu myrank
   integer :: ig_nsfsource = 0

!>number of receivers in hexahedron elements in cpu myrank
   integer :: ig_nreceiver_hexa = 0

!>number of receivers in quadrangle elements in cpu myrank
   integer :: ig_nreceiver_quad = 0

!>number of receivers in quadrangle elements for free surface snapshot in cpu myrank
   integer :: ig_nreceiver_snapshot = 0

!>media type: if = 0 --> material properties are identical for all GLL nodes of a hexahedron
   integer :: ig_medium_type = 0

!>number of material
   integer :: ig_nmaterial = 0

!>hexahedron element number attached to paraxial quadrangle elements
   integer, dimension(:), allocatable :: ig_quadp_neighbor_hexa

!>face number of hexahedron element attached to paraxial quadrangle elements
   integer, dimension(:), allocatable :: ig_quadp_neighbor_hexaface

!>number of neighbors (either face, edge and node) between cpu myrank and its neighbors
   integer :: ig_nneighbor_all_kind = 0

!>length of the name of cpu myrank
   integer :: ig_cpu_name_len = 0

!>name of cpu myrank
   character(len=MPI_MAX_PROCESSOR_NAME) :: cg_cpu_name

!>prefix of all input files
   character(len=92) :: cg_prefix

!>myrank number in character format
   character(len= 6) :: cg_myrank

!>total number of cpu in character format
   character(len= 6) :: cg_ncpu

!>if .true. --> snapshot of free surface is output either in gmt or VTK format
   logical :: lg_snapshot = .false.

!>if .true. --> activation of displacement free surface snapshot
   logical :: lg_snapshot_displacement = .false.

!>if .true. --> activation of velocity free surface snapshot
   logical :: lg_snapshot_velocity = .false.

!>if .true. --> activation of acceleration free surface snapshot
   logical :: lg_snapshot_acceleration = .false.

!>GLL nodes @f$x,y,z@f$coordinates in cpu myrank
   real, dimension(:,:), allocatable :: rg_gll_coordinate

!>weight of GLL nodes for numerical integration
   real, dimension(IG_NGLL) :: rg_gll_weight

!>abscissa of GLL nodes in reference domain [-1:1]
   real, dimension(IG_NGLL) :: rg_gll_abscissa

!>derivative of lagrange polynomial of node p at GLL node i in [-1:1] --> rg_gll_lagrange_deriv(p,i)
   real, dimension(IG_NGLL,IG_NGLL) :: rg_gll_lagrange_deriv

!>denominator for computing lagrange polynomial for GLL nodes in the reference domain [-1:1]
   real, dimension(IG_NGLL,IG_NGLL) :: rg_gll_abscissa_dist

!>total number of unique global GLL nodes in cpu myrank
   integer :: ig_ngll_total = 0

!
!
!***********************************************************************************************************************************************************************************
!single force point sources
!***********************************************************************************************************************************************************************************

!>type for single force point sources
   type type_single_force_source 
      real    :: x          !<@f$x@f$ coordinate of single force point source
      real    :: y          !<@f$y@f$ coordinate of single force point source
      real    :: z          !<@f$z@f$ coordinate of single force point source
      real    :: fac        !<factor on single force point source function
      real    :: rise_time  !<rise time of source function
      real    :: shift_time !<time shift of source function
      real    :: var1       !<other variable needed by some functions
      real    :: dmin       !<minimal distance of point source to a GLL node
      integer :: icur       !<function's number associated to single force source. See mod_source_function to see the definition of functions' number.
      integer :: idir       !<direction (@f$x@f$=1, @f$y@f$=2, @f$z@f$=3) in which the force is applied
      integer :: cpu        !<cpu that computes the single force point source
      integer :: iel        !<hexahedron element that contains the single force point source
      integer :: iequ       !<GLL node (global numbering) on which the force is applied
      integer :: kgll       !<closest local GLL node in the @f$\zeta@f$-direction
      integer :: lgll       !<closest local GLL node in the @f$\eta @f$-direction
      integer :: mgll       !<closest local GLL node in the @f$\xi  @f$-direction
   end type
!>data structure for single force point sources in cpu myrank 
   type(type_single_force_source), allocatable, dimension(:) :: tg_sfsource

!
!                                                      
!***********************************************************************************************************************************************************************************
!type for double couple sources
!***********************************************************************************************************************************************************************************

!>type for double couple point sources
   type type_double_couple_source                                            
      real    :: x          !<@f$x@f$ coordinate of double couple point source
      real    :: y          !<@f$y@f$ coordinate of double couple point source
      real    :: z          !<@f$z@f$ coordinate of double couple point source
      real    :: mxx        !<@f$xx@f$ component of the source tensor
      real    :: myy        !<@f$yy@f$ component of the source tensor
      real    :: mzz        !<@f$zz@f$ component of the source tensor
      real    :: mxy        !<@f$xy@f$ component of the source tensor
      real    :: mxz        !<@f$xz@f$ component of the source tensor
      real    :: myz        !<@f$yz@f$ component of the source tensor
      real    :: shift_time !<time shift of double couple source function
      real    :: rise_time  !<rise time of double couple source function                            
      real    :: xi         !<local coordinate @f$\xi@f$ of double couple source
      real    :: et         !<local coordinate @f$\eta@f$ of double couple source
      real    :: ze         !<local coordinate @f$\zeta@f$ of double couple source
      real    :: dmin       !<minimal distance of point source to a GLL node
      real    :: str        !<strike angle in degree of double couple source (starting clockwise from the @f$y@f$-axis (=north)). @f$ 0 \le str \le 360 @f$
      real    :: dip        !<dip angle in degree of double couple source @f$ 0 \le dip \le 90 @f$
      real    :: rak        !<rake angle in degree of double couple source @f$ -180 \le rak \le -180 @f$
      real    :: mw         !<moment magnitude @f$M_{\mathrm{w}}@f$ of double couple source
      integer :: icur       !<function's number associated to double couple source. See mod_source_function to see the definition of functions' number.
      integer :: cpu        !<cpu that computes the double couple source
      integer :: iel        !<hexahedron element that contains the double couple source
      integer :: kgll       !<closest local GLL node in @f$\zeta@f$-direction
      integer :: lgll       !<closest local GLL node in @f$\eta@f$-direction
      integer :: mgll       !<closest local GLL node in @f$\xi@f$-direction
   end type                   
!>data structure for double couple sources in cpu myrank                            
   type(type_double_couple_source), allocatable, dimension(:) :: tg_dcsource

!
!
!***********************************************************************************************************************************************************************************
!type for receivers inside hexahedron
!***********************************************************************************************************************************************************************************

!>type for receivers (i.e., stations) located inside hexahedron elements
   type type_receiver_hexa
      real                                        :: x       !<@f$x@f$ coordinate of receiver
      real                                        :: y       !<@f$y@f$ coordinate of receiver
      real                                        :: z       !<@f$z@f$ coordinate of receiver
      real                                        :: xi      !<local coordinate @f$\xi@f$ of receiver
      real                                        :: eta     !<local coordinate @f$\eta@f$ of receiver
      real                                        :: zeta    !<local coordinate @f$\zeta@f$ of receiver
      real                                        :: dmin    !<minimal distance of receiver to a GLL node
      real   , dimension(IG_NGLL,IG_NGLL,IG_NGLL) :: lag     !<Lagrange polynomial at receiver location (local coordinates) in hexahedron element iel
      integer, dimension(IG_NGLL,IG_NGLL,IG_NGLL) :: gll     !<global GLL nodes number used for Lagrange interpolation at the receiver location
      integer                                     :: cpu     !<cpu that contains the receiver
      integer                                     :: iel     !<hexahedron element that contains the receiver
      integer                                     :: kgll    !<closest local GLL node in @f$\zeta@f$-direction
      integer                                     :: lgll    !<closest local GLL node in @f$\eta@f$-direction
      integer                                     :: mgll    !<closest local GLL node in @f$\xi@f$-direction
      integer                                     :: rglo    !<global number of receiver among all cpus
   end type
!>data structure for receivers located inside hexahedron elements
   type(type_receiver_hexa), allocatable, dimension(:) :: tg_receiver_hexa

!
!
!***********************************************************************************************************************************************************************************
!type for receivers inside quad
!***********************************************************************************************************************************************************************************

!>type for receivers (i.e., stations) located inside quadrangle elements
   type type_receiver_quad
      real                                :: x       !<@f$x@f$ coordinate of receiver
      real                                :: y       !<@f$y@f$ coordinate of receiver
      real                                :: z       !<@f$z@f$ coordinate of receiver
      real                                :: xi      !<local coordinate @f$\xi@f$ of receiver
      real                                :: eta     !<local coordinate @f$\eta@f$ of receiver
      real                                :: dmin    !<minimal distance of receiver to GLL node
      real   , dimension(IG_NGLL,IG_NGLL) :: lag     !<Lagrange polynomial at receiver location (local coordinates) in quadrangle element iel
      real                                :: pgd_x   !<Peak Ground Displacement (PGD) in @f$x@f$-direction
      real                                :: pgd_y   !<Peak Ground Displacement (PGD) in @f$y@f$-direction
      real                                :: pgd_z   !<Peak Ground Displacement (PGD) in @f$z@f$-direction
      real                                :: pgv_x   !<Peak Ground Velocity (PGV) in @f$x@f$ direction
      real                                :: pgv_y   !<Peak Ground Velocity (PGV) in @f$y@f$ direction
      real                                :: pgv_z   !<Peak Ground Velocity (PGV) in @f$z@f$ direction
      real                                :: pgv_xyz !<PGV module = @f$max_{0 \le t_{simu} \le tt_{simu}}\left[ \sqrt{v_x(t)^2 + v_y(t)^2 + v_z(t)^2} \right] @f$
      real                                :: pga_x   !<Peak Ground Acceleration (PGA) in @f$x@f$ direction
      real                                :: pga_y   !<Peak Ground Acceleration (PGA) in @f$y@f$ direction
      real                                :: pga_z   !<Peak Ground Acceleration (PGA) in @f$z@f$ direction
      integer, dimension(IG_NGLL,IG_NGLL) :: gll     !<global GLL nodes number used for Lagrange interpolation at the receiver location
      integer                             :: cpu     !<cpu that contains the receiver
      integer                             :: iel     !<quadrangle element that contains the receiver
      integer                             :: lgll    !<closest local GLL node in @f$\eta@f$-direction
      integer                             :: mgll    !<closest local GLL node in @f$\xi@f$-direction
      integer                             :: rglo    !<global number of receiver among all cpus
   end type
!>data structure for receivers used for snapshot
   type(type_receiver_quad), allocatable, dimension(:) :: tg_receiver_snapshot_quad
!>data structure for receivers located on free surface quadrangle elements
   type(type_receiver_quad), allocatable, dimension(:) :: tg_receiver_quad

!
!
!***********************************************************************************************************************************************************************************
!type for cpus connected to cpu myrank
!***********************************************************************************************************************************************************************************

!>type that gathers information about cpus connected to cpu myrank (cpus not connected to cpu myrank are not stored)
   type type_cpu_neighbor
      integer                            :: ngll     !<total number of unique common GLL nodes between cpu myrank and cpu icon
      integer                            :: icpu     !<global cpu number (from call mpi_comm_rank) connected to cpu myrank
      integer, allocatable, dimension(:) :: gll_send !<GLL nodes number to be sent to cpu icon
      integer, allocatable, dimension(:) :: gll_recv !<GLL nodes number to be received from cpu icon
   endtype
!>data structure of cpu myrank that contains information about connected cpu icon.
   type(type_cpu_neighbor), allocatable, dimension(:) :: tg_cpu_neighbor

!
!
!***********************************************************************************************************************************************************************************
!type for linear elastic properties of materials
!***********************************************************************************************************************************************************************************

!>type for linear elastic properties of materials
   type type_elastic_material
      real :: dens !<density @f$\rho@f$
      real :: svel !<S-wave velocity @f$\beta@f$
      real :: pvel !<P-wave veocity @f$\alpha@f$
      real :: rvel !<Rayleigh wave velocity (approximated)
      real :: pois !<Poisson coefficient @f$\nu@f$
      real :: lam1 !<Lame parameter @f$\lambda@f$
      real :: lam2 !<Lame parameter @f$\mu = \rho \beta^2@f$ (S-wave modulus)
      real :: youn !<Young modulus
      real :: bulk !<Bulk modulus
      real :: pwmo !<P-wave modulus @f$\rho \alpha^2@f$
   end type
!>data structure for linear elastic properties of materials 
   type(type_elastic_material), allocatable, dimension(:) :: tg_elastic_material

!
!
!***********************************************************************************************************************************************************************************
!type for viscoelastic properties of materials
!***********************************************************************************************************************************************************************************

!>type for viscoelastic properties of materials
   type type_viscoelastic_material
      real                            :: freq  !<reference frequency @f$ f_r @f$ for which unrelaxed s-wave and p-wave modulus are defined
      real                            :: qs    !<s-wave quality factor @f$ Q_s @f$
      real                            :: qp    !<p-wave quality factor @f$ Q_p @f$
      real                            :: uswm  !<unrelaxed s-wave modulus @f$ M^{S}_{u} @f$
      real                            :: upwm  !<unrelaxed p-wave modulus @f$ M^{P}_{u} @f$
      real, dimension(:), allocatable :: wkqs  !<weight coefficients @f$ w^{S}_{k} @f$ associated to @f$ Q_s @f$
      real, dimension(:), allocatable :: wkqp  !<weight coefficients @f$ w^{P}_{k} @f$ associated to @f$ Q_p @f$
   end type
!>data structure for viscoelastic properties of materials
   type(type_viscoelastic_material), allocatable, dimension(:) :: tg_viscoelastic_material


!***********************************************************************************************************************************************************************************
!***********************************************************************************************************************************************************************************
!***********************************************************************************************************************************************************************************
!***********************************************************************************************************************************************************************************


   public  :: get_newunit
   public  :: error_stop
   public  :: info_all_cpu
   public  :: sweep_blanks
   public  :: strupcase
   public  :: strlowcase

   contains
!>@cond
!
!
!>@brief function to search for an available unit.
!!lun_min and lun_max define the range of possible luns to check.
!!the unit value is returned by the function, and also by the optional
!!argument. this allows the function to be used directly in an open
!!statement, and optionally save the result in a local variable.
!!if no units are available, -1 is returned.
!>@param myunit : free unit 
!***********************************************************************************************************************************************************************************
   integer function get_newunit(myunit)
!***********************************************************************************************************************************************************************************

      implicit none

      integer, intent(out), optional :: myunit
      integer, parameter             :: lun_min=101
      integer, parameter             :: lun_max=20101
      logical                        :: is_opened
      integer                        :: lun

      get_newunit=-1
      do lun=lun_min,lun_max
         inquire(unit=lun,opened=is_opened)
         if (.not. is_opened) then
            get_newunit=lun
            exit
         endif
      enddo
      if (present(myunit)) myunit=get_newunit

!***********************************************************************************************************************************************************************************
   end function get_newunit
!***********************************************************************************************************************************************************************************

!
!
!>@brief subroutine to abort and stop simulation with an error message
!>@param info : error message printed to standard output
!***********************************************************************************************************************************************************************************
   subroutine error_stop(info)
!***********************************************************************************************************************************************************************************

      use mpi

      implicit none

      character(len=*), intent(in) :: info
      integer                      :: ios

      write(*,'(a,i5,1x,a)') "cpu ",ig_myrank,trim(adjustl(info))

      call mpi_abort(mpi_comm_world,100,ios)
      stop

!***********************************************************************************************************************************************************************************
   end subroutine error_stop
!***********************************************************************************************************************************************************************************

!
!
!>@brief subroutine to gather and sum integer variables (done by cpu0 only)
!>@param i    : variable gathered by cpu0
!>@param info : name of gathered variables
!***********************************************************************************************************************************************************************************
   subroutine info_all_cpu(i,info)
!***********************************************************************************************************************************************************************************

      use mpi

      implicit none

      integer           , intent(in)       :: i
      character(len=255), intent(in)       :: info

      integer, dimension(ig_ncpu) :: buffer
      integer                              :: ios
      integer                              :: icpu
      integer(kind=8)                      :: isum

      isum = 0

      call mpi_gather(i,1,mpi_integer,buffer,1,mpi_integer,0,mpi_comm_world,ios)

      if (ig_myrank == 0) then

         write(IG_LST_UNIT,'(a)') " "

         do icpu = 1,ig_ncpu
            write(IG_LST_UNIT,'(a,i6,a,i12,1x,a)') "cpu ",icpu-1," computes      ",buffer(icpu),trim(info)
            isum = isum + buffer(icpu)
         enddo

         write(IG_LST_UNIT,'("              -----------------------",/,a,i29)') "total = ",isum

      endif

      return
!***********************************************************************************************************************************************************************************
   end subroutine info_all_cpu
!***********************************************************************************************************************************************************************************

!
!
!>@brief function to sweep blank characters in string
!>@param in_str       : input string to be swept
!***********************************************************************************************************************************************************************************
   character(len=100) function sweep_blanks(in_str)
!***********************************************************************************************************************************************************************************

     implicit none
     
     character(len=100), intent(in) :: in_str
     character(len=100)             :: out_str
     character(len=1)               :: ch
     integer                        :: j
     
     out_str = " "
     do j=1, len_trim(in_str)
!    
!----->get j-th char
       ch = in_str(j:j)
     
       if (ch .ne. " ") then
         out_str = trim(out_str) // ch
       endif
       sweep_blanks = out_str 
     enddo
!***********************************************************************************************************************************************************************************
  end function sweep_blanks
!***********************************************************************************************************************************************************************************

!
!
!>@brief function to convert lowercase characters to uppercase characters
!>@param input_string : lowercase input string 
!***********************************************************************************************************************************************************************************
   function strupcase (input_string) result (output_string)
!***********************************************************************************************************************************************************************************

      implicit none

      character(len=*), parameter      :: lower_case = 'abcdefghijklmnopqrstuvwxyz'
      character(len=*), parameter      :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

      character(len=*), intent( in )   :: input_string
      character(len=len(input_string)) :: output_string

      integer                          :: i
      integer                          :: n
      
      output_string = input_string

      do i = 1,len(output_string)
         n = index(lower_case,output_string(i:i))
         if (n /= 0) output_string(i:i) = upper_case(n:n)
      enddo

!***********************************************************************************************************************************************************************************
   end function strupcase
!***********************************************************************************************************************************************************************************

!
!
!>@brief function to convert uppercase characters to lowercase characters
!>@param input_string : uppercase input string 
!***********************************************************************************************************************************************************************************
   function strlowcase (input_string) result (output_string)
!***********************************************************************************************************************************************************************************

      implicit none

      character(len=*), parameter      :: lower_case = 'abcdefghijklmnopqrstuvwxyz'
      character(len=*), parameter      :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

      character(len=*), intent(in)     :: input_string
      character(len=len(input_string)) :: output_string

      integer                          :: i
      integer                          :: n
      
      output_string = input_string

      do i = 1,len(output_string)
         n = index(upper_case,output_string(i:i))
         if (n /= 0) output_string(i:i) = lower_case(n:n)
      enddo

!***********************************************************************************************************************************************************************************
   end function strlowcase
!***********************************************************************************************************************************************************************************
!>@endcond
   
end module mod_global_variables
